# Directional Human Connectome Inferred From Macaque Tracer Connectivity

This repository contains data and analysis code accompanying **Directed Human Structural Connectome Reveals Hierarchical Organization and Shapes Large-Scale Brain Dynamics**. The project estimates a directional human brain connectome by combining macaque tracer-derived connectivity with cross-species connectivity-blueprint mapping and human diffusion MRI tractography.

The released material is intended to support manuscript-level reproducibility and secondary reuse of the inferred connectome. The main data products are 426-node human structural connectivity matrices based on the HCP MMP 1.0 cortical atlas plus HCPex subcortical regions.

## Overview

The workflow maps directional connectivity information from the macaque brain to the human brain using tractography-based connectivity blueprints. In brief:

1. Connectivity blueprints are computed for human and macaque brain regions using major white-matter tracts that are homologous across species.
2. Human-to-macaque regional similarity is estimated from the blueprint profiles.
3. Macaque tracer-based directional connectivity, derived from the CocoMac/Modha-Singh compilation, is projected to the human parcellation.
4. The inferred directionality is combined with human diffusion MRI tractography weights and lengths to obtain directional human connectome matrices.

Please see the associated manuscript for methodological details, parameter choices, and validation analyses.

## Repository Structure

```text
.
|-- Data
|   |-- Nodes_informations.csv
|   |-- blueprints
|   |   |-- blueprint_human_HCPex_cortex_HCP-YA100.csv
|   |   |-- blueprint_human_HCPex_subcortex_HCP-YA100.csv
|   |   |-- blueprint_macaque_D99_cortex_BNA.csv
|   |   |-- blueprint_macaque_D99_subcortex_BNA.csv
|   |   `-- tract_order.txt
|   `-- connections
|       |-- DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv
|       |-- DirectionalConnectome_D99_HCPex_LR_length.csv
|       |-- DirectionalConnectome_D99_HCPex_left_dens20.csv
|       |-- DirectionalConnectome_D99_HCPex_right_dens20.csv
|       |-- DWI_connectivity_HCP-YA100_streamline.csv
|       |-- DWI_connectivity_HCP-YA100_streamline_length.csv
|       |-- Lh_connectionP.csv
|       |-- Rh_connectionP.csv
|       |-- ConnectionMatrix_Mac_D99.csv
|       |-- ConnectionList_ModhaDS_2010.csv
|       |-- ccep_probability.txt
|       `-- ccep_MNI-HCP-MMP1.txt
`-- Scripts
    |-- fig.1
    |-- fig.2,3
    |-- fig.4
    |-- fig.5
    `-- fig.6
```

## Main Data Products

### Node Information

`Data/Nodes_informations.csv` provides metadata for the 426 human nodes. It includes:

- node rank and atlas index;
- region names from the combined HCP MMP 1.0 + HCPex parcellation;
- network labels;
- MNI152-space coordinates;
- directed graph measures, including in-degree, out-degree, bidirectional degree, and trophic-level measures;
- cortical gradients, task-value annotations, and white-matter scores used in downstream analyses.

Matrix rows and columns should be interpreted according to the ranked node order in this file unless a matrix provides explicit row or column labels.

### Directional Human Connectome

`Data/connections/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv`

- Size: `426 x 426`.
- Direction: rows are source regions and columns are target regions.
- Values: directional connection weights, assigned from human diffusion MRI streamline counts as `log10(streamline count + 1)`.
- Intended use: primary weighted directional connectome, including left and right hemispheres.

`Data/connections/DirectionalConnectome_D99_HCPex_LR_length.csv`

- Size: `426 x 426` plus labels.
- Values: average streamline length between connected regions.
- Unit: millimeters.
- Intended use: tract-length/delay matrix for network modeling.

`Data/connections/DirectionalConnectome_D99_HCPex_left_dens20.csv` and `Data/connections/DirectionalConnectome_D99_HCPex_right_dens20.csv`

- Hemisphere-specific binary or thresholded directional matrices.
- The released matrices retain approximately 20% connection density.

### Human Diffusion MRI Connectome

`Data/connections/DWI_connectivity_HCP-YA100_streamline.csv`

- Human diffusion MRI streamline-count connectome.
- Based on the HCP-YA100 group data used in the manuscript workflow.

`Data/connections/DWI_connectivity_HCP-YA100_streamline_length.csv`

- Average streamline-length matrix corresponding to the diffusion MRI connectome.

### Macaque Tracer Connectivity

`Data/connections/ConnectionList_ModhaDS_2010.csv`

- Directed macaque connection list derived from the Modha and Singh compilation.
- Columns include source/target region identifiers and source-to-target or target-to-source directional flags.

`Data/connections/ConnectionMatrix_Mac_D99.csv`

- Macaque directional connectivity matrix after alignment to the D99-related region definitions used in this project.

### Connectivity Blueprints

The blueprint files in `Data/blueprints/` store region-by-tract connectivity profiles.

- Rows: brain regions.
- Columns: major white-matter tracts.
- Values: proportional tract connectivity profiles used to estimate cross-species regional similarity.

`tract_order.txt` defines the tract column order used by the blueprint scripts.

### CCEP Data

`Data/connections/ccep_probability.txt` and `Data/connections/ccep_MNI-HCP-MMP1.txt` contain CCEP-derived probability information and corresponding HCP MMP 1.0 region names used for validation analyses in `Scripts/fig.2,3/Compare_CCEP.R`.

## Scripts

The `Scripts/` directory contains analysis code organized by manuscript figure.

| Directory | Purpose |
| --- | --- |
| `Scripts/fig.1/` | Blueprint generation notes, MRtrix3 preprocessing scripts, cross-species direction inference, and supplementary weighting notebooks. |
| `Scripts/fig.2,3/` | CCEP comparison, connectome overlap, and graph-theory sensitivity analyses. |
| `Scripts/fig.4/` | Trophic-level analysis tools and notebook. |
| `Scripts/fig.5/` | Resting-state network-model simulations and analysis using TVB-based neural mass models. |
| `Scripts/fig.6/` | Stimulation simulations and analysis of stimulation-related signal features. |

The scripts are provided as research-analysis code. Some scripts contain local absolute paths from the original compute environment and should be edited before reuse.

## Basic Usage

### Python

```python
import pandas as pd

nodes = pd.read_csv("Data/Nodes_informations.csv")
weights = pd.read_csv(
    "Data/connections/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv",
    header=None,
)
lengths = pd.read_csv(
    "Data/connections/DirectionalConnectome_D99_HCPex_LR_length.csv",
    index_col=0,
)

print(nodes.shape)    # 426 nodes
print(weights.shape)  # 426 x 426
print(lengths.shape)  # 426 x 426
```

### R

```r
nodes <- read.csv("Data/Nodes_informations.csv")
weights <- read.csv(
  "Data/connections/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv",
  header = FALSE
)
lengths <- read.csv(
  "Data/connections/DirectionalConnectome_D99_HCPex_LR_length.csv",
  row.names = 1
)

dim(nodes)
dim(weights)
dim(lengths)
```

## Software Requirements

The released scripts use a mixture of R, Python, shell tools, and neuroimaging software. Based on the scripts in this repository, the following dependencies are required for full reproduction:

### Python

- `numpy`
- `pandas`
- `matplotlib`
- `networkx`
- `numba`
- `tvb` / The Virtual Brain simulator

### R

- `ggplot2`
- `pheatmap`
- `viridis`
- `igraph`
- `pracma`
- `expm`

### Neuroimaging / Command-Line Tools

- MRtrix3
- FSL / XTRACT
- FreeSurfer
- ANTs
- Workbench Command, if using HCP surface/volume mapping steps

Exact software versions should be added if required by the manuscript or journal policy.

## Reproducibility Notes

- The primary reusable connectome is stored in `Data/connections/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv`.
- The corresponding tract-length matrix is stored in `Data/connections/DirectionalConnectome_D99_HCPex_LR_length.csv`.
- Some scripts include hard-coded paths from the original analysis environment. Replace these paths with local paths before running.
- Figure 5 and Figure 6 simulation scripts include precomputed data files in their respective directories, but full reruns may require external subject-level diffusion or fMRI data that are not included here. All neuroimaging data used in this paper can be obtained at the following address: Huamn: https://www.humanconnectome.org/ ; Macaque: https://doi.org/10.57760/sciencedb.15197.
- The original HCP MMP 1.0 cortical labels were adjusted to match the left/right orientation used for the subcortical HCPex labels. See the manuscript for the exact atlas-handling procedure.

## Known Limitations

This connectome is an inferred directional structural connectome rather than a direct measurement of human axonal directionality. Important caveats include:

- **False negatives:** macaque tracer compilations may not contain all existing macaque connections, and tractography may miss some plausible human connections.
- **Cross-hemisphere projections:** directionality was inferred primarily from left-hemisphere mapping; cross-hemisphere connections were constrained using connection-weight thresholds.
- **Evolutionary differences:** human cortical expansions, especially in frontal and temporal regions, may contain connectivity patterns that are not fully captured by macaque-to-human mapping.
- **Threshold choice:** the density threshold used for the released binary directional matrices is a methodological choice and should be considered when comparing results across thresholds.

## References

Please cite the associated manuscript when using this repository:

> **[AUTHOR LIST TO BE ADDED]**. **[MANUSCRIPT TITLE TO BE ADDED]**. **[JOURNAL / PREPRINT SERVER TO BE ADDED]**, **[YEAR TO BE ADDED]**. **[DOI OR URL TO BE ADDED]**.

Key methodological references include:

- Glasser, M. F. et al. A multi-modal parcellation of human cerebral cortex. *Nature* 536, 171-178 (2016).
- Huang, C.-C., Rolls, E. T., Feng, J. & Lin, C.-P. An extended Human Connectome Project multimodal parcellation atlas of the human cortex and subcortical areas. *Brain Structure and Function* 227, 763-778 (2022).
- Mars, R. B. et al. Whole brain comparative anatomy using connectivity blueprints. *eLife* 7, e35237 (2018).
- Modha, D. S. & Singh, R. Network architecture of the long-distance pathways in the macaque brain. *Proceedings of the National Academy of Sciences* 107, 13485-13490 (2010).
- Warrington, S. et al. XTRACT - standardised protocols for automated tractography in the human and macaque brain. *NeuroImage* 217, 116923 (2020).
- Lu, Y. et al. Macaque Brainnetome Atlas: A multifaceted brain map with parcellation, connection, and histology. *Science Bulletin* 69, 2241-2259 (2024).


## License
Code is MIT licensed; data are released under CC BY 4.0

## Contact
For questions about the data or analysis workflow, please contact:
- Nan Huang: nan.HUANG@univ-amu.fr / nanhuang01@gmail.com
