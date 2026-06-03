from tvb.simulator.lab import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import time
from neuromass_models import MontbrioPazoRoxin_m
from m1_data_process import getFCD
from m2_simulation import *
import multiprocessing as mp

# read connectivity 
file_path = "./"
directional_conn = pd.read_csv(f"{file_path}/DirectionalConnectome_D99_HCPex_LR_Log10weighted_20250715.csv", header=None)
directional_conn = directional_conn.values
directional_conn_len = pd.read_csv(f"{file_path}/DirectionalConnectome_D99_HCPex_LR_length_20250715.csv", header=0, index_col=0)
directional_conn_len = directional_conn_len.values
nodirectional_conn = pd.read_csv(f"{file_path}/connectivity_HCP-YA100_streamline_Log10weighted_20250715.csv", header=0, index_col=0)
nodirectional_conn = nodirectional_conn.values
nodirectional_conn_len = pd.read_csv(f"{file_path}/connectivity_HCP-YA100_streamline_length_ranked_20250715.csv", header=0, index_col=0)
nodirectional_conn_len = nodirectional_conn_len.values
regions_coordinates = pd.read_csv(f"{file_path}/HCPmmp_HCPex_node_coordinate_ranked.csv")
region_names = np.array(regions_coordinates["HCP_MMP_name"].tolist())
regions_coordinate = regions_coordinates[['x','y','z']].values

# symmetrilized directional connectome
def compute_symmetric_connectome(directional_conn, nodirectional_conn, nodirectional_conn_len):
    mask = np.zeros(directional_conn.shape)
    mask[directional_conn > 0] = 1
    mask = mask + mask.T
    mask[mask > 0] = 1
    symm_directional_conn = nodirectional_conn * mask
    symm_directional_conn_len = nodirectional_conn_len * mask
    return symm_directional_conn, symm_directional_conn_len

# keep same in nodirectional_conn
def compute_sparse_connectome(directional_conn, nodirectional_conn, nodirectional_conn_len):
    mask = np.zeros(directional_conn.shape)
    mask[directional_conn > 0] = 1
    threshold = -np.sort(-nodirectional_conn[nodirectional_conn > 0])[int(np.sum(mask)) - 1]
    sparse_nodirectional_conn = np.where(nodirectional_conn >= threshold, nodirectional_conn, 0)
    np.fill_diagonal(sparse_nodirectional_conn, 0)
    mask = np.zeros(sparse_nodirectional_conn.shape)
    mask[sparse_nodirectional_conn > 0] = 1
    saparse_nodirectional_conn_len = nodirectional_conn_len * mask
    np.fill_diagonal(saparse_nodirectional_conn_len, 0)
    return sparse_nodirectional_conn, saparse_nodirectional_conn_len

# read subject list
subject_path = "/data3/nh_data/20250204_DirectionalConnectome/dataset/data/human/subject_7T.txt"
with open(subject_path, 'r', encoding='utf-8') as file:
    subject_list = file.readlines()
subject_list = [line.strip() for line in subject_list]




if __name__ == '__main__':
    for sample in subject_list[35:60]: # test 3 subjects
        # read personalized connectivity
        sunject_file_path = f"/data3/nh_data/20250204_DirectionalConnectome/dataset/data/result/H{sample}"
        if os.path.exists(f"{sunject_file_path}/directional_masked_connectome_weight.csv") == False:
            print(f"{sample} directional connectome not exist, skip!")
            continue
        personalized_directional_conn = pd.read_csv(f"{sunject_file_path}/directional_masked_connectome_weight.csv", header=0, index_col=0)
        personalized_directional_conn = personalized_directional_conn.values
        personalized_directional_conn_len = pd.read_csv(f"{sunject_file_path}/directional_masked_connectome_length.csv", header=0, index_col=0)
        personalized_directional_conn_len = personalized_directional_conn_len.values

        personalized_nodirectional_conn = pd.read_csv(f"{sunject_file_path}/connectome_weight_ranked.csv", header=0, index_col=0)
        personalized_nodirectional_conn = personalized_nodirectional_conn.values
        personalized_nodirectional_conn_len = pd.read_csv(f"{sunject_file_path}/connectome_length_ranked.csv", header=0, index_col=0)
        personalized_nodirectional_conn_len = personalized_nodirectional_conn_len.values

        # compute symmetrized and sparse connectome
        personalized_symm_directional_conn, personalized_symm_directional_conn_len = compute_symmetric_connectome(personalized_directional_conn, personalized_nodirectional_conn, personalized_nodirectional_conn_len)
        personalized_sparse_nodirectional_conn, personalized_saparse_nodirectional_conn_len = compute_sparse_connectome(personalized_directional_conn, personalized_nodirectional_conn, personalized_nodirectional_conn_len)

        # define model
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

        # define parameters for simulation
        G = np.arange(0.030, 0.0451, 0.0005)
        sigma = np.arange(0.030, 0.0451, 0.005)

        weights = personalized_sparse_nodirectional_conn/ np.max(personalized_sparse_nodirectional_conn)
        lengths = personalized_saparse_nodirectional_conn_len.copy()
        Duration = 61 * 1000.0  # ms * 0.1
        TR = 72  # ms * 0.1
        seed = 42
        window_length = 42  # frames
        overlap = 40  # 
        output_dir = f"/data3/nh_data/20250204_DirectionalConnectome/dataset/simulation_output/20250916_HCA-YA100_individual/threshold/{sample}"

        corrFSCs = np.zeros((len(G), len(sigma)))
        varFCDs = np.zeros((len(G), len(sigma)))
        params = [(i, j, G_val, sigma_val,region_names, model, weights,lengths, regions_coordinate, Duration, TR, seed, window_length, overlap, output_dir)
                    for j, sigma_val in enumerate(sigma)
                    for i, G_val in enumerate(G)]

        
        with mp.Pool(processes=12) as pool:
            results = pool.map(run_single_simulation, params)

        for i, j, corr, var in results:
            corrFSCs[i, j] = corr
            varFCDs[i, j] = var
        
        np.savez(f'{output_dir}/data/paramters_sweep_result.npz', corrFSCs=corrFSCs, varFCDs=varFCDs)
