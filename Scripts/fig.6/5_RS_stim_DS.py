# define simulation tool
from tvb.simulator.lab import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import time
from m3_simulation_stim import *
import multiprocessing as mp
import matplotlib

# read connectome data ------------------------------------------------------------------------------------------------
file_path = "./"
directional_conn = pd.read_csv(f"{file_path}/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv", header=None)
directional_conn = directional_conn.values
directional_conn_len = pd.read_csv(f"{file_path}/DirectionalConnectome_D99_HCPex_LR_length.csv", header=0, index_col=0)
directional_conn_len = directional_conn_len.values
nodirectional_conn = pd.read_csv(f"{file_path}/connectivity_HCP-YA100_streamline_Log10weighted.csv", header=0, index_col=0)
nodirectional_conn = nodirectional_conn.values
nodirectional_conn_len = pd.read_csv(f"{file_path}/connectivity_HCP-YA100_streamline_length_ranked.csv", header=0, index_col=0)
nodirectional_conn_len = nodirectional_conn_len.values
regions_coordinates = pd.read_csv(f"{file_path}/HCPmmp_HCPex_node_coordinate_ranked.csv")
region_names = np.array(regions_coordinates["HCP_MMP_name"].tolist())
regions_coordinate = regions_coordinates[['x','y','z']].values

# symmetrilized directional connectome
mask = np.zeros(directional_conn.shape)
mask[directional_conn > 0] = 1
mask = mask + mask.T
mask[mask > 0] = 1
symm_directional_conn = nodirectional_conn * mask
symm_directional_conn_len = nodirectional_conn_len * mask

# keep top 20% weights in nodirectional_conn
mask = np.zeros(directional_conn.shape)
mask[directional_conn > 0] = 1
threshold = np.percentile(nodirectional_conn[nodirectional_conn > 0], (100 - (np.sum(mask)/(directional_conn.shape[1]*directional_conn.shape[0]))*100) )
sparse_nodirectional_conn = np.where(nodirectional_conn >= threshold, nodirectional_conn, 0)
np.fill_diagonal(sparse_nodirectional_conn, 0)
mask = np.zeros(sparse_nodirectional_conn.shape)
mask[sparse_nodirectional_conn > 0] = 1
saparse_nodirectional_conn_len = nodirectional_conn_len * mask
np.fill_diagonal(saparse_nodirectional_conn_len, 0)
# ---------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    # define key elemets for simulation ----------------------------------------------------------------
    # define model parameters
    tau =np.array([1.])
    I =np.array([0.0])
    Delta= np.array([0.7])
    J =np.array([14.5])
    eta=np.array([-4.6])
    model = models.MontbrioPazoRoxin(J = J ,
                                tau = tau,
                                Delta = Delta,
                                eta = eta,
                                I = I,
                                cr = np.array([1.0]),
                                cv = np.array([0.0]),
                            )
    model.state_variable_range["V"] = np.array([-2., -2.])
    model.state_variable_range["r"] = np.array([0., 0.])

    # define connectivity directional
    weights = symm_directional_conn.copy() / np.max(symm_directional_conn)
    lengths = symm_directional_conn_len.copy()

    conn_weight = weights
    conn_length = lengths 
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

    # define stimulus
    dura_stim = 200.0 # duration of the stimulus (ms)
    pulse_train = equations.PulseTrain()
    pulse_train.parameters["tau"] = dura_stim # pulse width or pulse duration (ms)
    pulse_train.parameters["T"] = 100000.0 # pulse repetition period
    pulse_train.parameters["onset"] = 2000.0 # onset time  (ms)
    pulse_train.parameters["amp"] = 3.0 # amplitude of the pulse

    # based on the distance to select stimulation sites, select all nodes within 1mm as stimulation sites
    # stim_distribution = 0.0
    stim_distribution = 1.0
    # region_list = [0]
    region_list = [0, 23]
    for region in region_list: # test the first 10 regions
        region_seclect_p = region  # select the first region as an example 
        region_name = region_names[region_seclect_p]
        distances = np.linalg.norm(regions_coordinate - regions_coordinate[region_seclect_p], axis=1)
        stim_region = np.where(distances <= 1)[0] # select all regions within (1mm) as stimulation sites
        stim_region = np.append(stim_region, stim_region + 213) # # selec both hemispheres
        print(f"Stimulating region: {region_name}, index: {stim_region}, total {len(stim_region)} regions")
        stim_weights = np.zeros((nregions, 1))
        stim_weights[stim_region,0] = stim_distribution

        stimulus_pulse = patterns.StimuliRegion(temporal = pulse_train,
                                        connectivity = conn, 
                                        weight = stim_weights)
        # -----------------------------------------------------------------------------------------------
        seed_list = np.arange(1, 1001) # different seeds for different runs

        G_val = 0.036
        sigma_val = 0.050
        Duration = 8000  # duration of the simulation (ms)
        TR = 72  # repetition time of fMRI (ms)
        dir = "/data3/nh_data/20250204_DirectionalConnectome/dataset/simulation_output/20260323_stim/directed_symm"
        output_dir = f"{dir}/region_{region + 1}_{region_name}"
        if stim_distribution == 0:
            output_dir = f"{dir}/region_background"
        
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        params = [(model, conn, stimulus_pulse, G_val, sigma_val, Duration, TR, seed, output_dir)
                    for i, seed in enumerate(seed_list)]
        
        with mp.Pool(processes=10) as pool:  # 
            results = pool.map(run_single_simulation, params)