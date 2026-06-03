# define simulation tool
from tvb.simulator.lab import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import time
import multiprocessing as mp
import matplotlib

def stim_generator(model,conn,stimulus_pulse,
                      dt = 0.01, G_couping = 0.045, seed=42, nsigma=0.045, R_biase=1, V_biase=2, Duration=1000, TR = 72):
    T_mon = 1.00
    nregions = conn.weights.shape[0]
    #Initial conditions
    V_0=np.random.uniform(low=0., high=0., size=((1,1,nregions,1)))
    R_0=np.random.uniform(low=0, high=0, size=((1,1,nregions,1)))
    init_cond=np.concatenate([V_0, R_0], axis=1)

    # Step 2: Define the brain connectivity (use default structural connectivity)
    conn = conn
    # Step 3: Define the coupling between brain regions
    coup = coupling.Linear(a=np.array([G_couping]), b=np.array([0.0]))
    # Step 4: Set up the integrator for solving the differential equations
    nsigma = nsigma
    Heun_noise = integrators.HeunStochastic(
                dt=dt,
                noise=noise.Additive(
                    nsig=np.array(
                        [nsigma*R_biase,nsigma*V_biase]
                    ), 
                    noise_seed=seed)
            )
    # Step 5: Define the monitors to record the output (raw neural activity)
    mon_raw = monitors.TemporalAverage(period=T_mon)  # Monitor the V and M variables
    # mon_raw = monitors.Raw()
    # mon_bold = monitors.Bold(period=TR)
    mon = (mon_raw,) 
    if monitors is None:
        print("Error: monitors object is None.")
    # Step 6: Create the simulator with the custom model
    sim = simulator.Simulator(
        model=model,         # Use the custom neural mass model
        connectivity=conn,  # Brain connectivity
        coupling=coup,          # Coupling between regions
        integrator=Heun_noise,      # Time integration
        monitors=mon,           # Monitors to observe the data
        stimulus = stimulus_pulse
    )
    sim.initial_conditions = init_cond
    sim.configure()  # Finalize configuration
    (raw_time, raw_data), = sim.run(simulation_length=Duration)
    from m1_data_process import tavg_to_bold
    bold_time, bold_data = tavg_to_bold(raw_time[:], raw_data[:,:,:,:], tavg_period=T_mon*10, svar=0, decimate=TR)
    return raw_time, raw_data, bold_time, bold_data


def simulation_runner(initial_dt, model, conn, stimulus_pulse, 
                      G, nsigma, Duration, TR, seed, output_dir,min_dt = 0.001):
    dt = initial_dt
    seed = seed
    while dt >= min_dt:        
        tic = time.time()
        tavg_t, tavg_d, bold_t, bold_d  = stim_generator(model,conn,stimulus_pulse,
                      dt = dt, G_couping = G, seed=seed, nsigma=nsigma, R_biase=1, V_biase=2, Duration=Duration, TR = TR)
        
        print(' ⤷ Simulation for dt= %4.4fms required %0.2f minutes.' % (dt, (time.time()-tic)/60))
        #check isnan, and if yes, break the loop and continue with next G
        if np.isnan(tavg_d).any():
            print('isnan detected, reducing dt')
            # divide dt by 2
            dt /= 2
            continue
        else:

            # create a name for a file to save the figure
            folder = output_dir
            if not os.path.exists(folder):
                os.mkdir(folder)
            if not os.path.exists(f'{folder}/data'):
                os.mkdir(f'{folder}/data')
            # save the figure
            # save the data
            np.savez(f'{folder}/data/seed_{seed}.npz', tavg_t=tavg_t, tavg_d=tavg_d, bold_t=bold_t, bold_d=bold_d)

            break
    else:
        print('Minimum dt reached, moving to next random seed')


def run_single_simulation(args):
    model, conn, stimulus_pulse, G_val, sigma_val, Duration, TR, seed, output_dir = args
    
    simulation_runner(initial_dt = 0.01, 
                      model= model, 
                      conn = conn, 
                      stimulus_pulse = stimulus_pulse, 
                      G = G_val, 
                      nsigma = sigma_val, 
                      Duration = Duration, 
                      TR = TR, 
                      seed = seed, 
                      output_dir = output_dir,
                      min_dt = 0.001)


