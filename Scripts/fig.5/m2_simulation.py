# define simulation tool
from tvb.simulator.lab import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import time
from neuromass_models import MontbrioPazoRoxin_m
from neuromass_models import MontbrioPazoRoxin
from m1_data_process import getFCD
from tvb.simulator.backend.nb_mpr import NbMPRBackend


def resting_generator(model,region_names, conn_weight, conn_length, regions_coordinate, 
                      dt = 0.01, G_couping = 0.53, seed=42, nsigma=0.03, R_biase=1, V_biase=2, Duration=1000, TR = 72):
    T_mon = 1.0
    condspeed= 3.9
    conn_weight = conn_weight/np.max(conn_weight)
    conn = connectivity.Connectivity(
        weights = conn_weight.T,
        tract_lengths = conn_length.T,
        region_labels=np.array(region_names),
        centres = regions_coordinate
    )
    conn.centres_spherical()
    nregions = len(conn.weights)   
    conn.speed = np.array([condspeed])
    conn.configure()

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
    )
    sim.initial_conditions = init_cond
    sim.configure()  # Finalize configuration
    runner = NbMPRBackend()
    (raw_time, raw_data), = runner.run_sim(sim, simulation_length=Duration)
    from m1_data_process import tavg_to_bold
    bold_time, bold_data = tavg_to_bold(raw_time[:], raw_data[:,:,:,:], tavg_period=T_mon*10, svar=0, decimate=TR)
    return raw_time, raw_data, bold_time, bold_data

def simulation_runner(initial_dt, model, G, nsigma, Duration, TR, seed, output_dir,
                      weights, lengths, regions_coordinate,region_names, 
                      min_dt = 0.0001, bold_init = 6, Tinit = 0, window_length = 20, overlap = 19):
    dt = initial_dt
    seed = seed
    while dt >= min_dt:        
        tic = time.time()
        tavg_t, tavg_d, bold_t, bold_d  = resting_generator(model = model, 
                            region_names=region_names, 
                            conn_weight= weights, 
                            conn_length= lengths,
                            regions_coordinate=regions_coordinate,
                            G_couping = G,
                            nsigma = nsigma,
                            Duration= Duration,
                            seed=seed,
                            dt = dt,
                            TR=TR)
        
        print(' ⤷ Simulation for dt= %4.4fms required %0.2f minutes.' % (dt, (time.time()-tic)/60))

        #check isnan, and if yes, break the loop and continue with next G
        if np.isnan(tavg_d).any():
            print('isnan detected, reducing dt')
            # divide dt by 2
            dt /= 2
            continue
        else:
            # Store the results if no NaN is detected
             # get BOLD
            fMRI = bold_d[bold_init:,0,:,0]

            # compute FC and FCD
            FC = np.corrcoef(fMRI,rowvar=False)
            corrFSC = np.corrcoef(weights, FC)
            corrFSC = corrFSC[0,1]
            FCD, t_bold = getFCD(fMRI, window_length, overlap)[:2]
            IsupdiagFCD = np.triu_indices(np.shape(FCD)[0], overlap) 
            varFCD = np.var(FCD[IsupdiagFCD])

            # plot a figure with 4 subplots, with the top row plots having at least twice larger height
            fig, axs = plt.subplots(2, 2, figsize=(5, 5), gridspec_kw={'height_ratios': [1, 1]})
            fig.suptitle('G= ' + '%.4f'%G + ', dt= ' + '%4.4f'%dt) # , fontsize=20

            # Plot neural activity time series
            axs[0, 0].imshow(tavg_d[Tinit:, 0, :, 0].T, aspect='auto', cmap='viridis', vmin=0, vmax=1)
            axs[0, 0].set_title('neuronal activity', fontsize=9)
            axs[0, 0].set_xlabel('Time [s]', fontsize=8)
            axs[0, 0].axes.get_yaxis().set_visible(False)

            # Visualize the BOLD
            axs[0, 1].imshow(bold_d[bold_init:,0,:,0].T, aspect='auto', cmap='viridis', vmin=0, vmax=1)
            axs[0, 1].set_title('BOLD', fontsize=9)
            axs[0, 1].set_xlabel('Time [s]', fontsize=8)
            axs[0, 1].axes.get_yaxis().set_visible(False)

            # Visualize the FC
            cs1 = axs[1, 0].imshow(FC, cmap='jet', aspect='equal', interpolation='none')
            axs[1, 0].set_title('BOLD FC, corr(FC,SC)=' + '%.3f'%corrFSC, fontsize=9)
            axcb1 = fig.colorbar(cs1, ax=axs[1, 0], shrink=0.6)

            # Visualize the FCD
            cs2 = axs[1, 1].imshow(FCD, cmap='jet', aspect='equal', interpolation='none', extent=[t_bold[0], t_bold[-1], t_bold[-1], t_bold[0]])
            axs[1, 1].set_title('FCD, var(FCD)=' + '%.5f'%np.var(FCD[IsupdiagFCD]), fontsize=9)
            axs[1, 1].set_xlabel('samples', fontsize=9)
            axcb2 = fig.colorbar(cs2, ax=axs[1, 1], shrink=0.6)

            plt.tight_layout()

            # create a name for a file to save the figure
            folder = output_dir
            if not os.path.exists(folder):
                os.mkdir(folder)
            if not os.path.exists(f'{folder}/figure'):
                os.mkdir(f'{folder}/figure')
            if not os.path.exists(f'{folder}/data'):
                os.mkdir(f'{folder}/data')
            filename = 'test' + 'G'+'%.4f'%G + 'sigma' '%.3f'%nsigma 
            # save the figure
            fig.savefig(f'{folder}/figure/{filename}.png')
            # save the data
            np.savez(f'{folder}/data/{filename}.npz', tavg_t=tavg_t, tavg_d=tavg_d, bold_t=bold_t, bold_d=bold_d, 
                     FC=FC, FCD=FCD, corrFSC=corrFSC, varFCD=varFCD)

            break
    else:
        print('Failed to resolve NaN with minimum dt, no solution found for G=', G) # If the loop exits without breaking, it means min_dt was reached without success
        corrFSC = np.nan
        varFCD = np.nan
    return corrFSC, varFCD


def run_single_simulation(args):

    i, j, G_val, sigma_val,region_names, model, weights,lengths, regions_coordinate, Duration, TR, seed, window_length, overlap, output_dir = args
    
    corrFSC, varFCD = simulation_runner(
        region_names=region_names,
        initial_dt=0.01,
        model=model,
        G=G_val,
        nsigma=sigma_val,
        Duration=Duration,
        TR=TR,
        weights=weights,
        lengths=lengths,
        seed=seed,
        output_dir=output_dir,
        regions_coordinate=regions_coordinate,
        min_dt=0.004,
        bold_init=14,
        Tinit=1000,
        window_length=window_length,
        overlap=overlap
    )

    print(f"Completed G{i}={G_val:.1f}, sigma{j}={sigma_val:.3f}")
    return i, j, corrFSC, varFCD


