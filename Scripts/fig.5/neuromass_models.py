from tvb.simulator.models.base import Model, ModelNumbaDfun
from tvb.basic.neotraits.api import NArray, List, Range, Final
from tvb.simulator.backend.ref import RefBase
import numpy
import math
from numba import guvectorize, float64
def sigmoid(x, sig_a=1, sig_b=0, sig_c=1):
    return sig_c / (1 + numpy.exp(-sig_a*(x-sig_b))) - sig_c/2

def sigmoid1 (x, sig_a=1, sig_b=0, sig_c=1):
    return sig_c / (1 + numpy.exp(-sig_a*(x-sig_b)))


class MontbrioPazoRoxin(Model):
    r"""
    2D model describing the Ott-Antonsen reduction of infinite all-to-all
    coupled QIF neurons (Theta-neurons) as in [Montbrio_Pazo_Roxin_2015]_.

    The two state variables :math:`r` and :math:`V` represent the average
    firing rate and the average membrane potential of our QIF neurons.

    The equations of the infinite QIF 2D population model read

    .. math::
            \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
            \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
    
    Input from the network enters in the :math:`V` variable as 
    :math:`1/\tau(c_r C_r + c_v C_V)` where C is the incomming coupling. In 
    other words, depending on the parameters :math:`c_r`, :math:`c_v` we couple
    the neural masses via the firing rate and/or the membrane potential.
    
    .. [Montbrio_Pazo_Roxin_2015] Montbrió, E., Pazó, D., & Roxin, A. (2015). Macroscopic description for networks of spiking neurons. *Physical Review X*, 5(2), 021028.
    """

    # Define traited attributes for this model, these represent possible kwargs.

    tau = NArray(
        label=r":math:`\tau`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.001, hi=15.0, step=0.01),
        doc="""Characteristic time""",
    )

    I = NArray(
        label=":math:`I_{ext}`",
        default=numpy.array([0.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.01),
        doc="""External Current""",
    )

    Delta = NArray(
        label=r":math:`\Delta`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.0, hi=10.0, step=0.01),
        doc="""Mean heterogeneous noise""",
    )

    J = NArray(
        label=":math:`J`",
        default=numpy.array([15.0]),
        domain=Range(lo=-25.0, hi=25.0, step=0.0001),
        doc="""Mean Synaptic weight.""",
    )

    eta = NArray(
        label=r":math:`\eta`",
        default=numpy.array([-5.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.0001),
        doc="""Constant parameter to scale the rate of feedback from the
            firing rate variable to itself""",
    )

    Gamma = NArray(
        label=r":math:`\Gamma`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=10.0, step=0.01),
        doc="""Half-width of synaptic weight distribution""",
    )

    cr = NArray(
        label=":math:`cr`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    cv = NArray(
        label=":math:`cv`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable V.""",
    )

    

    # Informational attribute, used for phase-plane and initial()
    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"r": numpy.array([0., 2.0]),
                 "V": numpy.array([-2.0, 1.5])},
        doc="""Expected ranges of the state variables for initial condition generation and phase plane setup.""",
    )

    state_variable_boundaries = Final(
        label="State Variable boundaries [lo, hi]",
        default={
            "r": numpy.array([0.0, numpy.inf])
        },
    )

    # TODO should match cvars below..
    coupling_terms = Final(
        label="Coupling terms",
        # how to unpack coupling array
        default=["Coupling_Term_r", "Coupling_Term_V"]
    )

    state_variable_dfuns = Final(
        label="Drift functions",
        default={
            "r": "1/tau * ( Delta / (pi * tau) + 2 * V * r)",
            "V": "1/tau * ( V*V - pi*pi*tau*tau*r*r + eta + J * tau * r + I + cr * Coupling_Term_r + cv * Coupling_Term_V)"
        }
    )

    variables_of_interest = List(
        of=str,
        label="Variables or quantities available to Monitors",
        choices=("r", "V"),
        default=("r", "V"),
        doc="The quantities of interest for monitoring for the Infinite QIF 2D oscillator.",
    )

    parameter_names = List(
        of=str,
        label="List of parameters for this model",
        default='tau Delta eta J I cr cv'.split())

    state_variables = ('r', 'V')
    _nvar = 2
    # Cvar is the coupling variable. 
    cvar = numpy.array([0, 1], dtype=numpy.int32)
    # Stvar is the variable where stimulus is applied.
    stvar = numpy.array([1], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.0):
        r"""
            2D model describing the Ott-Antonsen reduction of infinite all-to-all
            coupled QIF neurons (Theta-neurons) as in [Montbrio_Pazo_Roxin_2015]_.

            The two state variables :math:`r` and :math:`V` represent the average
            firing rate and the average membrane potential of our QIF neurons.

            The equations of the infinite QIF 2D population model read

            .. math::
                    \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
                    \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
        """

        r, V = state_variables

        # [State_variables, nodes]
        I = self.I
        Delta = self.Delta
        Gamma = self.Gamma
        eta = self.eta
        tau = self.tau
        J = self.J
        cr = self.cr
        cv = self.cv


        Coupling_Term_r = coupling[0, :]  # This zero refers to the first element of cvar (r in this case)
        Coupling_Term_V = coupling[1, :]  # This zero refers to the second element of cvar (V in this case)
        derivative = numpy.empty_like(state_variables)

        derivative[0] = 1 / tau * (Delta / (numpy.pi * tau) + 2 * V * r)
        derivative[1] = 1 / tau * (
                    V ** 2 - numpy.pi ** 2 * tau ** 2 * r ** 2 + eta + J * tau * r + I + cr * Coupling_Term_r + cv * Coupling_Term_V)

        return derivative



class MontbrioPazoRoxin_m(Model):
    r"""
    2D model describing the Ott-Antonsen reduction of infinite all-to-all
    coupled QIF neurons (Theta-neurons) as in [Montbrio_Pazo_Roxin_2015]_.

    The two state variables :math:`r` and :math:`V` represent the average
    firing rate and the average membrane potential of our QIF neurons.

    The equations of the infinite QIF 2D population model read

    .. math::
            \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
            \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
    
    Input from the network enters in the :math:`V` variable as 
    :math:`1/\tau(c_r C_r + c_v C_V)` where C is the incomming coupling. In 
    other words, depending on the parameters :math:`c_r`, :math:`c_v` we couple
    the neural masses via the firing rate and/or the membrane potential.
    
    .. [Montbrio_Pazo_Roxin_2015] Montbrió, E., Pazó, D., & Roxin, A. (2015). Macroscopic description for networks of spiking neurons. *Physical Review X*, 5(2), 021028.
    """

    # Define traited attributes for this model, these represent possible kwargs.

    tau = NArray(
        label=r":math:`\tau`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.001, hi=15.0, step=0.01),
        doc="""Characteristic time""",
    )

    I = NArray(
        label=":math:`I_{ext}`",
        default=numpy.array([0.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.01),
        doc="""External Current""",
    )

    Delta = NArray(
        label=r":math:`\Delta`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.0, hi=10.0, step=0.01),
        doc="""Mean heterogeneous noise""",
    )

    J = NArray(
        label=":math:`J`",
        default=numpy.array([15.0]),
        domain=Range(lo=-25.0, hi=25.0, step=0.0001),
        doc="""Mean Synaptic weight.""",
    )

    eta = NArray(
        label=r":math:`\eta`",
        default=numpy.array([-5.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.0001),
        doc="""Constant parameter to scale the rate of feedback from the
            firing rate variable to itself""",
    )

    Gamma = NArray(
        label=r":math:`\Gamma`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=10.0, step=0.01),
        doc="""Half-width of synaptic weight distribution""",
    )

    cr = NArray(
        label=":math:`cr`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    cv = NArray(
        label=":math:`cv`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable V.""",
    )

    sig_a=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([1.0]),
        domain=Range(lo=-3.0, hi=3.0, step=0.1),
        doc="""scale inputs step.""") 
    sig_b=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([0.0]),
        domain=Range(lo=-3.0, hi=3.0, step=0.1),
        doc="""scale inputs (x left or right).""") 
    sig_c=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([2.0]),
        domain=Range(lo=0.0, hi=3.0, step=0.1),
        doc="""max couping input.""") 
    

    # Informational attribute, used for phase-plane and initial()
    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"r": numpy.array([0., 2.0]),
                 "V": numpy.array([-2.0, 1.5])},
        doc="""Expected ranges of the state variables for initial condition generation and phase plane setup.""",
    )

    state_variable_boundaries = Final(
        label="State Variable boundaries [lo, hi]",
        default={
            "r": numpy.array([0.0, numpy.inf])
        },
    )

    # TODO should match cvars below..
    coupling_terms = Final(
        label="Coupling terms",
        # how to unpack coupling array
        default=["Coupling_Term_r", "Coupling_Term_V"]
    )

    state_variable_dfuns = Final(
        label="Drift functions",
        default={
            "r": "1/tau * ( Delta / (pi * tau) + 2 * V * r)",
            "V": "1/tau * ( V*V - pi*pi*tau*tau*r*r + eta + J * tau * r + I + cr * Coupling_Term_r + cv * Coupling_Term_V)"
        }
    )

    variables_of_interest = List(
        of=str,
        label="Variables or quantities available to Monitors",
        choices=("r", "V"),
        default=("r", "V"),
        doc="The quantities of interest for monitoring for the Infinite QIF 2D oscillator.",
    )

    parameter_names = List(
        of=str,
        label="List of parameters for this model",
        default='tau Delta eta J I cr cv'.split())

    state_variables = ('r', 'V')
    _nvar = 2
    # Cvar is the coupling variable. 
    cvar = numpy.array([0, 1], dtype=numpy.int32)
    # Stvar is the variable where stimulus is applied.
    stvar = numpy.array([1], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.0):
        r"""
            2D model describing the Ott-Antonsen reduction of infinite all-to-all
            coupled QIF neurons (Theta-neurons) as in [Montbrio_Pazo_Roxin_2015]_.

            The two state variables :math:`r` and :math:`V` represent the average
            firing rate and the average membrane potential of our QIF neurons.

            The equations of the infinite QIF 2D population model read

            .. math::
                    \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
                    \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
        """

        r, V = state_variables

        # [State_variables, nodes]
        I = self.I
        Delta = self.Delta
        Gamma = self.Gamma
        eta = self.eta
        tau = self.tau
        J = self.J
        cr = self.cr
        cv = self.cv
        sig_a = self.sig_a  
        sig_b = self.sig_b
        sig_c = self.sig_c

        Coupling_Term_r = coupling[0, :]  # This zero refers to the first element of cvar (r in this case)
        Coupling_Term_V = coupling[1, :]  # This zero refers to the second element of cvar (V in this case)
        Coupling_Term_r = sigmoid(Coupling_Term_r, sig_a=sig_a, sig_b=sig_b, sig_c=sig_c)

        derivative = numpy.empty_like(state_variables)

        derivative[0] = 1 / tau * (Delta / (numpy.pi * tau) + 2 * V * r)
        derivative[1] = 1 / tau * (
                    V ** 2 - numpy.pi ** 2 * tau ** 2 * r ** 2 + eta + J * tau * r + I + cr * Coupling_Term_r + cv * Coupling_Term_V)

        return derivative



class Montbrio_Receptors(Model):
    r"""
    This model add transmitter feature on montbrio model

    .. math::
            \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
            \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
    
    Input from the network enters in the :math:`V` variable as 
    :math:`1/\tau(c_r C_r + c_v C_V)` where C is the incomming coupling. In 
    other words, depending on the parameters :math:`c_r`, :math:`c_v` we couple
    the neural masses via the firing rate and/or the membrane potential.
    
    .. [Montbrio_Pazo_Roxin_2015] Montbrió, E., Pazó, D., & Roxin, A. (2015). Macroscopic description for networks of spiking neurons. *Physical Review X*, 5(2), 021028.
    """

    # Define traited attributes for this model, these represent possible kwargs.

    tau = NArray(
        label=r":math:`\tau`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.001, hi=15.0, step=0.01),
        doc="""Characteristic time""",
    )

    tau_Sa = NArray(
        label=r":math:`\tau_Sa`",
        default=numpy.array([5.0]),
        domain=Range(lo=0.001, hi=15.0, step=0.01),
        doc="""Characteristic time""",
    )

    tau_Sg = NArray(
        label=r":math:`\tau_Sg`",
        default=numpy.array([5.0]),
        domain=Range(lo=0.001, hi=15.0, step=0.01),
        doc="""Characteristic time""",
    )

    I = NArray(
        label=":math:`I_{ext}`",
        default=numpy.array([0.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.01),
        doc="""External Current""",
    )

    Delta = NArray(
        label=r":math:`\Delta`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.0, hi=10.0, step=0.01),
        doc="""Mean heterogeneous noise""",
    )

    J = NArray(
        label=":math:`J`",
        default=numpy.array([15.0]),
        domain=Range(lo=-50.0, hi=50.0, step=0.0001),
        doc="""Mean Synaptic weight.""",
    )

    eta = NArray(
        label=r":math:`\eta`",
        default=numpy.array([-5.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.0001),
        doc="""Constant parameter to scale the rate of feedback from the
            firing rate variable to itself""",
    )

    Gamma = NArray(
        label=r":math:`\Gamma`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=10.0, step=0.01),
        doc="""Half-width of synaptic weight distribution""",
    )

    cr = NArray(
        label=":math:`cr`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    cv = NArray(
        label=":math:`cv`",
        default=numpy.array([0.]),
        domain=Range(lo=0., hi=1, step=0.1),
        doc="""It is the weight on Coupling through variable V.""",
    )

    c_glu = NArray(
        label=":math:`c_glu`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    c_gaba = NArray(
        label=":math:`c_gaba`",
        default=numpy.array([1.]),
        domain=Range(lo=-2., hi=0., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    r_glu = NArray(
        label=":math:`r_glu`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    r_gaba = NArray(
        label=":math:`r_gaba`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    E_ampa = NArray(
        label=":math:`E_ampa`",
        default=numpy.array([10.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc=""
        "It is the weight on Coupling through variable r.""",
    )

    E_gabaa = NArray(
        label=":math:`E_gabaa`",
        default=numpy.array([-10.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    g_ampa = NArray(
        label=":math:`g_ampa`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc=""
        "It is the weight on Coupling through variable r.""",
    )

    g_gabaa = NArray(
        label=":math:`g_gabaa`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""It is the weight on Coupling through variable r.""",
    )

    gama_ampa = NArray(
        label=":math:`gama_ampa`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""control the open rate of a AMPA recptor""",
    )

    gama_gabaa = NArray(
        label=":math:`gama_ampa`",
        default=numpy.array([1.]),
        domain=Range(lo=0., hi=1., step=0.1),
        doc="""control the open rate of a gabaa recptor""",
    )

    # Informational attribute, used for phase-plane and initial()
    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"r": numpy.array([0., 2.0]),
                 "V": numpy.array([-2.0, 1.5]),
                 "rel_glu": numpy.array([0., 1.0]),
                 "rel_gaba": numpy.array([0., 1.0]),
                 "Sa": numpy.array([0., 1.0]),
                 "Sg": numpy.array([0., 1.0]),},
                 
        doc="""Expected ranges of the state variables for initial condition generation and phase plane setup.""",
    )

    state_variable_boundaries = Final(
        label="State Variable boundaries [lo, hi]",
        default={
            "r": numpy.array([0.0, numpy.inf]),
            "rel_glu": numpy.array([0., numpy.inf]),
            "rel_gaba": numpy.array([0., numpy.inf])
        },
    )

    # TODO should match cvars below..
    coupling_terms = Final(
        label="Coupling terms",
        # how to unpack coupling array
        default=["Coupling_Term_r", "Coupling_Term_V", "Coupling_Term_glu", "Coupling_Term_gaba"]
    )

    state_variable_dfuns = Final(
        label="Drift functions",
        default={
            "r": "1/tau * ( Delta / (pi * tau) + 2 * V * r)",
            "V": "1/tau * ( V*V - pi*pi*tau*tau*r*r + eta + J * tau * r + I +  + )",
            "rel_glu": "r_glu * r",
            "rel_gaba": "r_gaba * r",
            "Sa": "-Sa + c_glu * Coupling_Term_glu",
            "Sg": "-Sg + c_gaba * Coupling_Term_gaba"
        }
    )

    variables_of_interest = List(
        of=str,
        label="Variables or quantities available to Monitors",
        choices=("r", "V", "rel_glu", "rel_gaba", "Sa", "Sg"),
        default=("r", "V", "rel_glu", "rel_gaba", "Sa", "Sg"),
        doc="The quantities of interest for monitoring for the Infinite QIF 2D oscillator.",
    )

    parameter_names = List(
        of=str,
        label="List of parameters for this model",
        default='tau tau_Sa tau_Sg Delta eta J I cr cv c_glu c_gaba r_glu r_gaba E_ampa E_gaba g_ampa g_gabaa gama_ampa gama_gabaa'.split())

    state_variables = ('r', 'V', 'rel_glu', 'rel_gaba', "Sa", "Sg")
    _nvar = 6
    # Cvar is the coupling variable. 
    cvar = numpy.array([2 ,3], dtype=numpy.int32)
    # Stvar is the variable where stimulus is applied.
    stvar = numpy.array([1], dtype=numpy.int32)

    def dfun(self, state_variables, coupling, local_coupling=0.0):
        r"""
            2D model describing the Ott-Antonsen reduction of infinite all-to-all
            coupled QIF neurons (Theta-neurons) as in [Montbrio_Pazo_Roxin_2015]_.

            The two state variables :math:`r` and :math:`V` represent the average
            firing rate and the average membrane potential of our QIF neurons.

            The equations of the infinite QIF 2D population model read

            .. math::
                    \dot{r} &= 1/\tau (\Delta/(\pi \tau) + 2 V r)\\
                    \dot{V} &= 1/\tau (V^2 - \tau^2 \pi^2 r^2 + \eta + J \tau r + I)
        """

        r, V, rel_glu, rel_gaba, Sa, Sg = state_variables

        # [State_variables, nodes]
        I = self.I
        Delta = self.Delta
        Gamma = self.Gamma
        eta = self.eta
        tau = self.tau
        J = self.J
        c_glu = self.c_glu
        c_gaba = self.c_gaba
        r_glu = self.r_glu
        r_gaba = self.r_gaba
        tau_Sa = self.tau_Sa
        tau_Sg = self.tau_Sg
        E_ampa = self.E_ampa
        E_gaba = self.E_gabaa
        g_ampa = self.g_ampa
        g_gabaa = self.g_gabaa
        gama_ampa = self.gama_ampa
        gama_gabaa = self.gama_gabaa


        Coupling_Term_glu = coupling[0, :]
        Coupling_Term_gaba = coupling[1, :]

        Input_Term = I + Sa * E_ampa * g_ampa + Sg * E_gaba * g_gabaa
        response_Trem = Sa * g_ampa * rel_glu + Sg * g_gabaa * rel_gaba 

        dr = 1 / tau * (Delta / (numpy.pi * tau) + 2 * V * r - response_Trem)
        dV = 1 / tau * (
                    V ** 2 - numpy.pi ** 2 * tau ** 2 * r ** 2 + eta + J * tau * r + Input_Term)
        dr_glu = r_glu * dr
        dr_gaba = r_gaba * dr
        dSa = -Sa/tau_Sa + (c_glu * Coupling_Term_glu + J * rel_glu) * (1-Sa) * gama_ampa
        dSg = -Sg/tau_Sg + (c_gaba * Coupling_Term_gaba + J * rel_gaba) * (1-Sg) * gama_gabaa

        derivative = numpy.empty_like(state_variables)
        derivative[0] = dr
        derivative[1] = dV
        derivative[2] = dr_glu
        derivative[3] = dr_gaba
        derivative[4] = dSa
        derivative[5] = dSg

        return derivative
    
# modifed Generic2dOscillator to add restrict of input (scale)
class Generic2dOscillator_m(ModelNumbaDfun):
    r"""
    The Generic2dOscillator model is a generic dynamic system with two state
    variables. The dynamic equations of this model are composed of two ordinary
    differential equations comprising two nullclines. The first nullcline is a
    cubic function as it is found in most neuron and population models; the
    second nullcline is arbitrarily configurable as a polynomial function up to
    second order. The manipulation of the latter nullcline's parameters allows
    to generate a wide range of different behaviours.

    Equations:

    .. math::
            \dot{V} &= d \, \tau (-f V^3 + e V^2 + g V + \alpha W + \gamma I) \\
            \dot{W} &= \dfrac{d}{\tau}\,\,(c V^2 + b V - \beta W + a)

    See:


        .. [FH_1961] FitzHugh, R., *Impulses and physiological states in theoretical
            models of nerve membrane*, Biophysical Journal 1: 445, 1961.

        .. [Nagumo_1962] Nagumo et.al, *An Active Pulse Transmission Line Simulating
            Nerve Axon*, Proceedings of the IRE 50: 2061, 1962.

        .. [SJ_2011] Stefanescu, R., Jirsa, V.K. *Reduced representations of
            heterogeneous mixed neural networks with synaptic coupling*.
            Physical Review E, 83, 2011.

        .. [SJ_2010]	Jirsa VK, Stefanescu R.  *Neural population modes capture
            biologically realistic large-scale network dynamics*. Bulletin of
            Mathematical Biology, 2010.

        .. [SJ_2008_a] Stefanescu, R., Jirsa, V.K. *A low dimensional description
            of globally coupled heterogeneous neural networks of excitatory and
            inhibitory neurons*. PLoS Computational Biology, 4(11), 2008).


    The model's (:math:`V`, :math:`W`) time series and phase-plane its nullclines
    can be seen in the figure below.

    The model with its default parameters exhibits FitzHugh-Nagumo like dynamics.

    +---------------------------+
    |  Table 1                  |
    +--------------+------------+
    |  EXCITABLE CONFIGURATION  |
    +--------------+------------+
    |Parameter     |  Value     |
    +==============+============+
    | a            |     -2.0   |
    +--------------+------------+
    | b            |    -10.0   |
    +--------------+------------+
    | c            |      0.0   |
    +--------------+------------+
    | d            |      0.02  |
    +--------------+------------+
    | I            |      0.0   |
    +--------------+------------+
    |  limit cycle if a is 2.0  |
    +---------------------------+


    +---------------------------+
    |   Table 2                 |
    +--------------+------------+
    |   BISTABLE CONFIGURATION  |
    +--------------+------------+
    |Parameter     |  Value     |
    +==============+============+
    | a            |      1.0   |
    +--------------+------------+
    | b            |      0.0   |
    +--------------+------------+
    | c            |     -5.0   |
    +--------------+------------+
    | d            |      0.02  |
    +--------------+------------+
    | I            |      0.0   |
    +--------------+------------+
    | monostable regime:        |
    | fixed point if Iext=-2.0  |
    | limit cycle if Iext=-1.0  |
    +---------------------------+


    +---------------------------+
    |  Table 3                  |
    +--------------+------------+
    |  EXCITABLE CONFIGURATION  |
    +--------------+------------+
    |  (similar to Morris-Lecar)|
    +--------------+------------+
    |Parameter     |  Value     |
    +==============+============+
    | a            |      0.5   |
    +--------------+------------+
    | b            |      0.6   |
    +--------------+------------+
    | c            |     -4.0   |
    +--------------+------------+
    | d            |      0.02  |
    +--------------+------------+
    | I            |      0.0   |
    +--------------+------------+
    | excitable regime if b=0.6 |
    | oscillatory if b=0.4      |
    +---------------------------+


    +---------------------------+
    |  Table 4                  |
    +--------------+------------+
    |  GhoshetAl,  2008         |
    |  KnocketAl,  2009         |
    +--------------+------------+
    |Parameter     |  Value     |
    +==============+============+
    | a            |    1.05    |
    +--------------+------------+
    | b            |   -1.00    |
    +--------------+------------+
    | c            |    0.0     |
    +--------------+------------+
    | d            |    0.1     |
    +--------------+------------+
    | I            |    0.0     |
    +--------------+------------+
    | alpha        |    1.0     |
    +--------------+------------+
    | beta         |    0.2     |
    +--------------+------------+
    | gamma        |    -1.0    |
    +--------------+------------+
    | e            |    0.0     |
    +--------------+------------+
    | g            |    1.0     |
    +--------------+------------+
    | f            |    1/3     |
    +--------------+------------+
    | tau          |    1.25    |
    +--------------+------------+
    |                           |
    |  frequency peak at 10Hz   |
    |                           |
    +---------------------------+


    +---------------------------+
    |  Table 5                  |
    +--------------+------------+
    |  SanzLeonetAl  2013       |
    +--------------+------------+
    |Parameter     |  Value     |
    +==============+============+
    | a            |    - 0.5   |
    +--------------+------------+
    | b            |    -10.0   |
    +--------------+------------+
    | c            |      0.0   |
    +--------------+------------+
    | d            |      0.02  |
    +--------------+------------+
    | I            |      0.0   |
    +--------------+------------+
    |                           |
    |  intrinsic frequency is   |
    |  approx 10 Hz             |
    |                           |
    +---------------------------+

    NOTE: This regime, if I = 2.1, is called subthreshold regime.
    Unstable oscillations appear through a subcritical Hopf bifurcation.


    .. figure :: img/Generic2dOscillator_01_mode_0_pplane.svg
    .. _phase-plane-Generic2D:
        :alt: Phase plane of the generic 2D population model with (V, W)

        The (:math:`V`, :math:`W`) phase-plane for the generic 2D population
        model for default parameters. The dynamical system has an equilibrium
        point.

    .. automethod:: Generic2dOscillator.dfun

    """

    # Define traited attributes for this model, these represent possible kwargs.
    tau = NArray(
        label=r":math:`\tau`",
        default=numpy.array([1.0]),
        domain=Range(lo=1.0, hi=5.0, step=0.01),
        doc="""A time-scale hierarchy can be introduced for the state
        variables :math:`V` and :math:`W`. Default parameter is 1, which means
        no time-scale hierarchy.""")

    I = NArray(
        label=":math:`I_{ext}`",
        default=numpy.array([0.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.01),
        doc="""Baseline shift of the cubic nullcline""")

    a = NArray(
        label=":math:`a`",
        default=numpy.array([-2.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.01),
        doc="""Vertical shift of the configurable nullcline""")

    b = NArray(
        label=":math:`b`",
        default=numpy.array([-10.0]),
        domain=Range(lo=-20.0, hi=15.0, step=0.01),
        doc="""Linear slope of the configurable nullcline""")

    c = NArray(
        label=":math:`c`",
        default=numpy.array([0.0]),
        domain=Range(lo=-10.0, hi=10.0, step=0.01),
        doc="""Parabolic term of the configurable nullcline""")

    d = NArray(
        label=":math:`d`",
        default=numpy.array([0.02]),
        domain=Range(lo=0.0001, hi=1.0, step=0.0001),
        doc="""Temporal scale factor. Warning: do not use it unless
        you know what you are doing and know about time tides.""")

    e = NArray(
        label=":math:`e`",
        default=numpy.array([3.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.0001),
        doc="""Coefficient of the quadratic term of the cubic nullcline.""")

    f = NArray(
        label=":math:`f`",
        default=numpy.array([1.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.0001),
        doc="""Coefficient of the cubic term of the cubic nullcline.""")

    g = NArray(
        label=":math:`g`",
        default=numpy.array([0.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.5),
        doc="""Coefficient of the linear term of the cubic nullcline.""")

    alpha = NArray(
        label=r":math:`\alpha`",
        default=numpy.array([1.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.0001),
        doc="""Constant parameter to scale the rate of feedback from the
            slow variable to the fast variable.""")

    beta = NArray(
        label=r":math:`\beta`",
        default=numpy.array([1.0]),
        domain=Range(lo=-5.0, hi=5.0, step=0.0001),
        doc="""Constant parameter to scale the rate of feedback from the
            slow variable to itself""")

    # This parameter is basically a hack to avoid having a negative lower boundary in the global coupling strength.
    gamma = NArray(
        label=r":math:`\gamma`",
        default=numpy.array([1.0]),
        domain=Range(lo=-1.0, hi=1.0, step=0.1),
        doc="""Constant parameter to reproduce FHN dynamics where
               excitatory input currents are negative.
               It scales both I and the long range coupling term.""")
    
    sig_a=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([1.0]),
        domain=Range(lo=-3.0, hi=3.0, step=0.1),
        doc="""scale inputs step.""") 
    sig_b=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([0.0]),
        domain=Range(lo=-3.0, hi=3.0, step=0.1),
        doc="""scale inputs (x left or right).""") 
    sig_c=NArray(
        label=r":math:`\sig_a`",
        default=numpy.array([2.0]),
        domain=Range(lo=0.0, hi=3.0, step=0.1),
        doc="""max couping input.""") 

    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"V": numpy.array([-2.0, 4.0]),
                 "W": numpy.array([-6.0, 6.0])},
        doc="""The values for each state-variable should be set to encompass
            the expected dynamic range of that state-variable for the current
            parameters, it is used as a mechanism for bounding random initial
            conditions when the simulation isn't started from an explicit
            history, it is also provides the default range of phase-plane plots.""")

    variables_of_interest = List(
        of=str,
        label="Variables or quantities available to Monitors",
        choices=("V", "W", "V + W", "V - W"),
        default=("V",),
        doc="The quantities of interest for monitoring for the generic 2D oscillator.")

    state_variables = ('V', 'W')
    _nvar = 2
    cvar = numpy.array([0], dtype=numpy.int32)

    def _numpy_dfun(self, state_variables, coupling, local_coupling=0.0):
        V = state_variables[0, :]
        W = state_variables[1, :]

        # [State_variables, nodes]
        c_0 = coupling[0, :]

        tau = self.tau
        I = self.I
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        f = self.f
        g = self.g
        beta = self.beta
        alpha = self.alpha
        gamma = self.gamma
        sig_c = self.sig_c
        sig_a = self.sig_a  
        sig_b = self.sig_b

        lc_0 = local_coupling * V

        # Pre-allocate the result array then instruct numexpr to use it as output.
        # This avoids an expensive array concatenation
        derivative = numpy.empty_like(state_variables)

        ev = RefBase.evaluate
        ev('d * tau * (alpha * W - f * V**3 + e * V**2 + g * V + gamma * I + gamma *c_0 + lc_0)', out=derivative[0])
        ev('d * (a + b * V + c * V**2 - beta * W) / tau', out=derivative[1])

        return derivative

    def dfun(self, vw, c, local_coupling=0.0):
        r"""
        The two state variables :math:`V` and :math:`W` are typically considered
        to represent a function of the neuron's membrane potential, such as the
        firing rate or dendritic currents, and a recovery variable, respectively.
        If there is a time scale hierarchy, then typically :math:`V` is faster
        than :math:`W` corresponding to a value of :math:`\tau` greater than 1.

        The equations of the generic 2D population model read

        .. math::
                \dot{V} &= d \, \tau (-f V^3 + e V^2 + g V + \alpha W + \gamma I) \\
                \dot{W} &= \dfrac{d}{\tau}\,\,(c V^2 + b V - \beta W + a)

        where external currents :math:`I` provide the entry point for local,
        long-range connectivity and stimulation.

        """

        sig_a = self.sig_a  
        sig_b = self.sig_b
        sig_c = self.sig_c
        
        lc_0 = local_coupling * vw[0, :, 0]
        vw_ = vw.reshape(vw.shape[:-1]).T
        c_ = c.reshape(c.shape[:-1]).T
        c_ = sigmoid(c_, sig_a=sig_a, sig_b=sig_b, sig_c=sig_c)
        deriv = _numba_dfun_g2d(vw_, c_, self.tau, self.I, self.a, self.b, self.c, self.d, self.e, self.f, self.g,
                                self.beta, self.alpha, self.gamma, lc_0)
        return deriv.T[..., numpy.newaxis]


@guvectorize([(float64[:],) * 16], '(n),(m)' + ',()' * 13 + '->(n)', nopython=True)
def _numba_dfun_g2d(vw, c_0, tau, I, a, b, c, d, e, f, g, beta, alpha, gamma, lc_0, dx):
    "Gufunc for Generic2dOscillator model equations."
    V = vw[0]
    V2 = V * V
    W = vw[1]
    dx[0] = d[0] * tau[0] * (
                alpha[0] * W - f[0] * V2 * V + e[0] * V2 + g[0] * V + gamma[0] * I[0] + gamma[0] * c_0[0] + lc_0[0])
    dx[1] = d[0] * (a[0] + b[0] * V + c[0] * V2 - beta[0] * W) / tau[0]



class SupHopf(ModelNumbaDfun):
    r"""
    The supHopf model describes the normal form of a supercritical Hopf bifurcation in Cartesian coordinates.
    This normal form has a supercritical bifurcation at a=0 with a the bifurcation parameter in the model. So 
    for a < 0, the local dynamics has a stable fixed point and the system corresponds to a damped oscillatory 
    state, whereas for a > 0, the local dynamics enters in a stable limit cycle and the system switches to an 
    oscillatory state.

    See for examples:

    .. [Kuznetsov_2013] Kuznetsov, Y.A. *Elements of applied bifurcation theory.* Springer Sci & Business
        Media, 2013, vol. 112.

    .. [Deco_2017a] Deco, G., Kringelbach, M.L., Jirsa, V.K., Ritter, P. *The dynamics of resting fluctuations
       in the brain: metastability and its dynamical cortical core* Sci Reports, 2017, 7: 3095.

    The equations of the supHopf equations read as follows:

    .. math::
        \dot{x}_{i} &= (a_{i} - x_{i}^{2} - y_{i}^{2})x_{i} - {\omega}{i}y_{i} \\
        \dot{y}_{i} &= (a_{i} - x_{i}^{2} - y_{i}^{2})y_{i} + {\omega}{i}x_{i}

    where a is the local bifurcation parameter and omega the angular frequency.
    """

    a = NArray(
        label=r":math:`a`",
        default=numpy.array([-0.5]),
        domain=Range(lo=-10.0, hi=10.0, step=0.01),
        doc="""Local bifurcation parameter.""")

    omega = NArray(
        label=r":math:`\omega`",
        default=numpy.array([1.]),
        domain=Range(lo=0.05, hi=630.0, step=0.01),
        doc="""Angular frequency.""")

    # Initialization.
    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"x": numpy.array([-5.0, 5.0]),
                 "y": numpy.array([-5.0, 5.0])},
        doc="""The values for each state-variable should be set to encompass
               the expected dynamic range of that state-variable for the current
               parameters, it is used as a mechanism for bounding random initial
               conditions when the simulation isn't started from an explicit
               history, it is also provides the default range of phase-plane plots.""")

    variables_of_interest = List(
        of=str,
        label="Variables watched by Monitors",
        choices=("x", "y"),
        default=("x",),
        doc="Quantities of supHopf available to monitor.")

    state_variables = ["x", "y"]

    _nvar = 2  # number of state-variables
    cvar = numpy.array([0, 1], dtype=numpy.int32)  # coupling variables

    def _numpy_dfun(self, state_variables, coupling, local_coupling=0.0,
                    array=numpy.array, where=numpy.where, concat=numpy.concatenate):
        y = state_variables
        ydot = numpy.empty_like(state_variables)

        # long-range coupling
        c_0 = coupling[0]
        c_1 = coupling[1]

        # short-range (local) coupling
        lc_0 = local_coupling * y[0]

        # supHopf's equations in Cartesian coordinates:
        ydot[0] = (self.a - y[0] ** 2 - y[1] ** 2) * y[0] - self.omega * y[1] + c_0 + lc_0
        ydot[1] = (self.a - y[0] ** 2 - y[1] ** 2) * y[1] + self.omega * y[0] + c_1

        return ydot

    def dfun(self, x, c, local_coupling=0.0):
        r"""
        Computes the derivatives of the state-variables of supHopf
        with respect to time.

        The equations of the supHopf equations read as follows:

        .. math::
            \dot{x}_{i} &= (a_{i} - x_{i}^{2} - y_{i}^{2})x_{i} - {\omega}{i}y_{i} \\
            \dot{y}_{i} &= (a_{i} - x_{i}^{2} - y_{i}^{2})y_{i} + {\omega}{i}x_{i}

        where a is the local bifurcation parameter and omega the angular frequency.
        """
        x_ = x.reshape(x.shape[:-1]).T
        c_ = c.reshape(c.shape[:-1]).T
        lc_0 = local_coupling * x[0, :, 0]
        deriv = _numba_dfun_supHopf(x_, c_, self.a, self.omega, lc_0)

        return deriv.T[..., numpy.newaxis]


@guvectorize([(float64[:],) * 6], '(n),(m)' + ',()' * 3 + '->(n)', nopython=True)
def _numba_dfun_supHopf(y, c, a, omega, lc_0, ydot):
    "Gufunc for supHopf model equations."

    # long-range coupling
    c_0 = c[0]
    c_1 = c[1]

    # supHopf equations
    ydot[0] = (a[0] - y[0] ** 2 - y[1] ** 2) * y[0] - omega[0] * y[1] + c_0 + lc_0[0]
    ydot[1] = (a[0] - y[0] ** 2 - y[1] ** 2) * y[1] + omega[0] * y[0] + c_1



class JansenRit_m(ModelNumbaDfun):
    r"""
    The Jansen and Rit is a biologically inspired mathematical framework
    originally conceived to simulate the spontaneous electrical activity of
    neuronal assemblies, with a particular focus on alpha activity, for instance,
    as measured by EEG. Later on, it was discovered that in addition to alpha
    activity, this model was also able to simulate evoked potentials.

    .. [JR_1995]  Jansen, B., H. and Rit V., G., *Electroencephalogram and
        visual evoked potential generation in a mathematical model of
        coupled cortical columns*, Biological Cybernetics (73) 357:366, 1995.

    .. [J_1993] Jansen, B., Zouridakis, G. and Brandt, M., *A
        neurophysiologically-based mathematical model of flash visual evoked
        potentials*

    .. figure :: img/JansenRit_45_mode_0_pplane.svg
        :alt: Jansen and Rit phase plane (y4, y5)

        The (:math:`y_4`, :math:`y_5`) phase-plane for the Jansen and Rit model.

    The dynamic equations were taken from [JR_1995]_

    .. math::
        \dot{y_0} &= y_3 \\
        \dot{y_3} &= A a\,S[y_1 - y_2] - 2a\,y_3 - a^2\, y_0 \\
        \dot{y_1} &= y_4\\
        \dot{y_4} &= A a \,[p(t) + \alpha_2 J + S[\alpha_1 J\,y_0]+ c_0]
                    -2a\,y - a^2\,y_1 \\
        \dot{y_2} &= y_5 \\
        \dot{y_5} &= B b (\alpha_4 J\, S[\alpha_3 J \,y_0]) - 2 b\, y_5
                    - b^2\,y_2 \\
        S[v] &= \frac{2\, \nu_{max}}{1 + \exp^{r(v_0 - v)}}

    """

    # Define traited attributes for this model, these represent possible kwargs.
    A = NArray(
        label=":math:`A`",
        default=numpy.array([3.25]),
        domain=Range(lo=2.6, hi=9.75, step=0.05),
        doc="""Maximum amplitude of EPSP [mV]. Also called average synaptic gain.""")

    B = NArray(
        label=":math:`B`",
        default=numpy.array([22.0]),
        domain=Range(lo=17.6, hi=110.0, step=0.2),
        doc="""Maximum amplitude of IPSP [mV]. Also called average synaptic gain.""")

    a = NArray(
        label=":math:`a`",
        default=numpy.array([0.1]),
        domain=Range(lo=0.05, hi=0.15, step=0.01),
        doc="""Reciprocal of the time constant of passive membrane and all
        other spatially distributed delays in the dendritic network [ms^-1].
        Also called average synaptic time constant.""")

    b = NArray(
        label=":math:`b`",
        default=numpy.array([0.05]),
        domain=Range(lo=0.025, hi=0.075, step=0.005),
        doc="""Reciprocal of the time constant of passive membrane and all
        other spatially distributed delays in the dendritic network [ms^-1].
        Also called average synaptic time constant.""")

    v0 = NArray(
        label=":math:`v_0`",
        default=numpy.array([5.52]),
        domain=Range(lo=3.12, hi=6.0, step=0.02),
        doc="""Firing threshold (PSP) for which a 50% firing rate is achieved.
        In other words, it is the value of the average membrane potential
        corresponding to the inflection point of the sigmoid [mV].

        The usual value for this parameter is 6.0.""")

    nu_max = NArray(
        label=r":math:`\nu_{max}`",
        default=numpy.array([0.0025]),
        domain=Range(lo=0.00125, hi=0.00375, step=0.00001),
        doc="""Determines the maximum firing rate of the neural population
        [ms^-1].""")

    r = NArray(
        label=":math:`r`",
        default=numpy.array([0.56]),
        domain=Range(lo=0.28, hi=0.84, step=0.01),
        doc="""Steepness of the sigmoidal transformation [mV^-1].""")

    J = NArray(
        label=":math:`J`",
        default=numpy.array([135.0]),
        domain=Range(lo=65.0, hi=1350.0, step=1.),
        doc="""Average number of synapses between populations.""")

    a_1 = NArray(
        label=r":math:`\alpha_1`",
        default=numpy.array([1.0]),
        domain=Range(lo=0.5, hi=1.5, step=0.1),
        doc="""Average probability of synaptic contacts in the feedback excitatory loop.""")

    a_2 = NArray(
        label=r":math:`\alpha_2`",
        default=numpy.array([0.8]),
        domain=Range(lo=0.4, hi=1.2, step=0.1),
        doc="""Average probability of synaptic contacts in the slow feedback excitatory loop.""")

    a_3 = NArray(
        label=r":math:`\alpha_3`",
        default=numpy.array([0.25]),
        domain=Range(lo=0.125, hi=0.375, step=0.005),
        doc="""Average probability of synaptic contacts in the feedback inhibitory loop.""")

    a_4 = NArray(
        label=r":math:`\alpha_4`",
        default=numpy.array([0.25]),
        domain=Range(lo=0.125, hi=0.375, step=0.005),
        doc="""Average probability of synaptic contacts in the slow feedback inhibitory loop.""")

    p_min = NArray(
        label=":math:`p_{min}`",
        default=numpy.array([0.12]),
        domain=Range(lo=0.0, hi=0.12, step=0.01),
        doc="""Minimum input firing rate.""")

    p_max = NArray(
        label=":math:`p_{max}`",
        default=numpy.array([0.32]),
        domain=Range(lo=0.0, hi=0.32, step=0.01),
        doc="""Maximum input firing rate.""")

    mu = NArray(
        label=r":math:`\mu_{max}`",
        default=numpy.array([0.22]),
        domain=Range(lo=0.0, hi=0.22, step=0.01),
        doc="""Mean input firing rate""")

    # Used for phase-plane axis ranges and to bound random initial() conditions.
    state_variable_range = Final(
        label="State Variable ranges [lo, hi]",
        default={"y0": numpy.array([-1.0, 1.0]),
                 "y1": numpy.array([-500.0, 500.0]),
                 "y2": numpy.array([-50.0, 50.0]),
                 "y3": numpy.array([-6.0, 6.0]),
                 "y4": numpy.array([-20.0, 20.0]),
                 "y5": numpy.array([-500.0, 500.0])},
        doc="""The values for each state-variable should be set to encompass
        the expected dynamic range of that state-variable for the current
        parameters, it is used as a mechanism for bounding random inital
        conditions when the simulation isn't started from an explicit history,
        it is also provides the default range of phase-plane plots.""")

    variables_of_interest = List(
        of=str,
        label="Variables watched by Monitors",
        choices=("y0", "y1", "y2", "y3", "y4", "y5"),
        default=("y0", "y1", "y2", "y3"),
        doc="""This represents the default state-variables of this Model to be
                                    monitored. It can be overridden for each Monitor if desired. The
                                    corresponding state-variable indices for this model are :math:`y0 = 0`,
                                    :math:`y1 = 1`, :math:`y2 = 2`, :math:`y3 = 3`, :math:`y4 = 4`, and
                                    :math:`y5 = 5`""")

    state_variables = tuple('y0 y1 y2 y3 y4 y5'.split())
    _nvar = 6
    cvar = numpy.array([1, 2], dtype=numpy.int32)

    def _numpy_dfun(self, state_variables, coupling, local_coupling=0.0):
        y0, y1, y2, y3, y4, y5 = state_variables

        # NOTE: This is assumed to be \sum_j u_kj * S[y_{1_j} - y_{2_j}]
        lrc = coupling[0, :]
        short_range_coupling = local_coupling*(y1 - y2)

        # NOTE: for local couplings
        # 0: pyramidal cells
        # 1: excitatory interneurons
        # 2: inhibitory interneurons
        # 0 -> 1,
        # 0 -> 2,
        # 1 -> 0,
        # 2 -> 0,

        exp = numpy.exp
        sigm_y1_y2 = 2.0 * self.nu_max / (1.0 + exp(self.r * (self.v0 - (y1 - y2))))
        sigm_y0_1  = 2.0 * self.nu_max / (1.0 + exp(self.r * (self.v0 - (self.a_1 * self.J * y0))))
        sigm_y0_3  = 2.0 * self.nu_max / (1.0 + exp(self.r * (self.v0 - (self.a_3 * self.J * y0))))

        return numpy.array([
            y3,
            y4,
            y5,
            self.A * self.a * sigm_y1_y2 - 2.0 * self.a * y3 - self.a ** 2 * y0,
            self.A * self.a * (self.mu + self.a_2 * self.J * sigm_y0_1 + lrc + short_range_coupling)
                - 2.0 * self.a * y4 - self.a ** 2 * y1,
            self.B * self.b * (self.a_4 * self.J * sigm_y0_3) - 2.0 * self.b * y5 - self.b ** 2 * y2,
        ])

    def dfun(self, y, c, local_coupling=0.0):
        r"""
        The dynamic equations were taken from [JR_1995]_

        .. math::
            \dot{y_0} &= y_3 \\
            \dot{y_3} &= A a\,S[y_1 - y_2] - 2a\,y_3 - 2a^2\, y_0 \\
            \dot{y_1} &= y_4\\
            \dot{y_4} &= A a \,[p(t) + \alpha_2 J S[\alpha_1 J\,y_0]+ c_0]
                        -2a\,y - a^2\,y_1 \\
            \dot{y_2} &= y_5 \\
            \dot{y_5} &= B b (\alpha_4 J\, S[\alpha_3 J \,y_0]) - 2 b\, y_5
                        - b^2\,y_2 \\
            S[v] &= \frac{2\, \nu_{max}}{1 + \exp^{r(v_0 - v)}}


        :math:`p(t)` can be any arbitrary function, including white noise or
        random numbers taken from a uniform distribution, representing a pulse
        density with an amplitude varying between 120 and 320

        For Evoked Potentials, a transient component of the input,
        representing the impulse density attribuable to a brief visual input is
        applied. Time should be in seconds.

        .. math::
            p(t) = q\,(\frac{t}{w})^n \, \exp{-\frac{t}{w}} \\
            q = 0.5 \\
            n = 7 \\
            w = 0.005 [s]

        """
        src =  local_coupling*(y[1] - y[2])[:, 0]
        y_ = y.reshape(y.shape[:-1]).T
        c_ = c.reshape(c.shape[:-1]).T
        deriv = _numba_dfun_jr(y_, c_, src,
                               self.nu_max, self.r, self.v0, self.a, self.a_1, self.a_2, self.a_3, self.a_4,
                               self.A, self.b, self.B, self.J, self.mu
                               )
        return deriv.T[..., numpy.newaxis]


@guvectorize([(float64[:],) * 17], '(n),(m)' + ',()'*14 + '->(n)', nopython=True)
def _numba_dfun_jr(y, c,
                   src,
                   nu_max, r, v0, a, a_1, a_2, a_3, a_4, A, b, B, J, mu,
                   dx):
    sigm_y1_y2 = 2.0 * nu_max[0] / (1.0 + math.exp(r[0] * (v0[0] - (y[1] - y[2]))))
    sigm_y0_1 = 2.0 * nu_max[0] / (1.0 + math.exp(r[0] * (v0[0] - (a_1[0] * J[0] * y[0]))))
    sigm_y0_3 = 2.0 * nu_max[0] / (1.0 + math.exp(r[0] * (v0[0] - (a_3[0] * J[0] * y[0]))))
    dx[0] = y[3]
    dx[1] = y[4]
    dx[2] = y[5]
    dx[3] = A[0] * a[0] * sigm_y1_y2 - 2.0 * a[0] * y[3] - a[0] ** 2 * y[0]
    dx[4] = A[0] * a[0] * (mu[0] + a_2[0] * J[0] * sigm_y0_1 + c[0] + src[0]) - 2.0 * a[0] * y[4] - a[0] ** 2 * y[1]
    dx[5] = B[0] * b[0] * (a_4[0] * J[0] * sigm_y0_3) - 2.0 * b[0] * y[5] - b[0] ** 2 * y[2]

