# InfoThermo — Information thermodynamics of network responses

Research code accompanying **Information-thermodynamic principles unveil the response of complex networks to perturbations**, by Andrea Filippo Beretta and Manlio De Domenico.

The code integrates network dynamics, follows the response of small perturbations through a density-matrix representation, and computes entropy, entropy rate, informational potential and distance from an instantaneous Gibbs reference. It allows to reproduce the panels used in Figs. 1a, 1b and 2 of the manuscript.

## Contents

- [Repository structure](#Repository-structure)
- [Code and datasets](#code-and-datasets)
- [Requirements and workflow](#requirements-and-workflow)
- [Reproduction guide](docs/REPRODUCIBILITY.md)
- [Models and output formats](docs/MODELS_AND_DATA.md)
- [Citation and license](#citation-and-license)

## Repository structure

```text
.
├── README.md
├── docs/
│   ├── REPRODUCIBILITY.md       # Instructions for reproducing the figures
│   ├── MODELS_AND_DATA.md       # Dynamical models and output formats
│   └── LICENSE.md               # MIT License
│   
│
├── synth_biochemical/           # Fig. 1a: synthetic networks
│   ├── synthgen_ensemble.py     # Generate synthetic network ensembles
│   ├── timediff.ipynb          # Compute dynamics and entropy-rate peak times
│   └── plot_distr.ipynb        # Plot peak-time difference distributions
│
├── Neuronal/                    # Fig. 1b: C. elegans
│   ├── neuronal.npy            # Empirical adjacency matrix
│   ├── dataneuro3.py           # Simulate R = 0.10, 0.20, ..., 0.50
│   ├── dataneuro3_1.py         # Simulate R = 0.05, 0.15, ..., 0.45
│   ├── plots3.ipynb            # Combine both runs and generate the figure
│   ├── entropyplot.ipynb       # Visualize entropy trajectories
│   └── vizneuro.ipynb          # Additional single-run diagnostics
│
├── Noisy/                       # Fig. 2a: noisy voter dynamics
│   ├── noisy.npy              # Empirical adjacency matrix
│   ├── ERandCM.ipynb           # Generate ER and CM network ensembles
│   ├── dataNoisy.ipynb         # Compute trajectories and thermodynamics
│   └── plots2.ipynb            # Compare empirical, ER and CM results
│
├── Epidemics/                   # Fig. 2b: SIS dynamics
│   ├── epidemic.npy           # Empirical adjacency matrix
│   ├── ERandCM.ipynb           # Generate ER and CM network ensembles
│   ├── dataEpi.ipynb           # Compute trajectories and thermodynamics
│   └── plots2.ipynb            # Compare empirical, ER and CM results
│
└── Population/                 # Fig. 2c: population dynamics
    ├── pop.npy                # Empirical adjacency matrix
    ├── ERandCM.ipynb           # Generate ER and CM network ensembles
    ├── dataPop.ipynb           # Compute trajectories and thermodynamics
    └── plots2.ipynb            # Compare empirical, ER and CM results
```


## Code and datasets

| Folder | Analysis | Manuscript panel | Included input |
| --- | --- | --- | --- |
| `synth_biochemical/` | Biochemical dynamics versus diffusion on ER, BA and SBM networks | Fig. 1a | 
| `Neuronal/` | Neuronal dynamics versus diffusion on C. elegans | Fig. 1b | `neuronal.npy` |
| `Noisy/` | Noisy voter dynamics on the dolphin network and ER/CM ensembles | Fig. 2a | `noisy.npy` |
| `Epidemics/` | SIS dynamics on a contact network and ER/CM ensembles | Fig. 2b | `epidemic.npy` |
| `Population/` | Population dynamics on a trophic network and ER/CM ensembles | Fig. 2c | `pop.npy` |

This is a collection of Python scripts and Jupyter notebooks, not an installable Python package. The four empirical input matrices are included; simulation outputs are generated locally. 

## Requirements and workflow

Required: Python 3, NumPy, SciPy, NetworkX, Matplotlib, tqdm and Jupyter. **graph-tool is required** by the three `ERandCM.ipynb` notebooks to generate configuration-model ensembles.


Consult the [official graph-tool installation instructions](https://graph-tool.skewed.de/installation.html) for other installation methods. Use the same environment for Python and Jupyter.


| Workflow | Execution order within its folder |
| --- | --- |
| Synthetic | `synthgen_ensemble.py` → `timediff.ipynb` → `plot_distr.ipynb` |
| Neuronal | Both `dataneuro3.py` and `dataneuro3_1.py` → `plots3.ipynb` |
| Noisy voter | `ERandCM.ipynb` → `dataNoisy.ipynb` → `plots2.ipynb` |
| SIS | `ERandCM.ipynb` → `dataEpi.ipynb` → `plots2.ipynb` |
| Population | `ERandCM.ipynb` → `dataPop.ipynb` → `plots2.ipynb` |

Generated trajectories, observables and PDFs are not included and are regenerated during execution. Keep outputs inside their respective folders. The plot notebooks export panels and legends; final multi-panel page assembly is not automated.

## Citation and license

Andrea Filippo Beretta and Manlio De Domenico. *Information-thermodynamic principles unveil the response of complex networks to perturbations*. This project is licensed under the MIT License - see the [LICENSE](docs/LICENSE.md) file for details. 

## Contact

Andrea Filippo Beretta: andreafilippo.beretta@studenti.unipd.it. 
