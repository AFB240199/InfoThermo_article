# Models and output formats

## What the code computes

The system follows an ODE, dx/dt = v(x; A), on a fixed adjacency matrix. The Jacobian J along the trajectory governs infinitesimal perturbations. A statistical propagator is evolved using the flow equation and normalized to a density matrix rho. The symmetric part of J is D = (J + Jᵀ)/2.

The thermodynamic routines compute:

| Quantity | Definition | Common saved names |
| --- | --- | --- |
| Entropy | S = -Tr(rho log rho) | `S`, `entropy` |
| Entropy rate | -2[Tr(rho D log rho) - Tr(rho D)Tr(rho log rho)] | `Sdot`, `entropy_rate` |
| Informational potential | U = -2 Tr(rho D) | `E`, `energy` |
| Relative entropy to gamma | D_KL = t U - S + log Z_gamma | `div`, `divergence` |
| Free-energy difference | Delta F = D_KL / t | Computed in plotting notebooks |

Here gamma is proportional to exp(2tD). The matrix named `U` in the integration routines is the statistical propagator; it is distinct from the scalar potential denoted U in the manuscript.

## Dynamical models

Let (Ax)_i denote the i-th component of A times x. Powers and tanh act componentwise.

| Folder | Implemented dynamics |
| --- | --- |
| `synth_biochemical/` | dx/dt = F - Bx - R x ⊙ (Ax) |
| `Neuronal/` | dx/dt = -Bx + F tanh(x) + R A tanh(x); F = 0 in the supplied runs |
| `Epidemics/` | dx/dt = -Bx + R(1-x) ⊙ (Ax) |
| `Noisy/` | dx_i/dt = F - Bx_i + R(Ax)_i/k_i, with zero coupling for k_i = 0 |
| `Population/` | dx/dt = -B x^b + R A x^a |

For noisy voter dynamics k_i is the row sum of A. 

## Included empirical matrices

| File | Shape | Identification in the manuscript |
| --- | --- | --- |
| `Neuronal/neuronal.npy` | 279 × 279 | C. elegans nervous system |
| `Epidemics/epidemic.npy` | 95 × 95 | Contact network |
| `Noisy/noisy.npy` | 62 × 62 | Dolphin social network |
| `Population/pop.npy` | 183 × 183 | Trophic network |

The arrays are numeric NumPy `.npy` files.

## Generated data

### Synthetic workflow

`network_ensembles.npz` contains `A_er`, `A_ws`, `A_ba`, `A_sbm`, `A_rg` and seeds, together with generation parameters. The analysis uses ER, BA and SBM. Their trajectory archives contain `t`, `x`, `A`, `F`, `B` and `R`. Each `delta_t_bio_<model>.npz` contains `delta_t`, `t_peak_J` and `t_peak_L`. The plot uses absolute `delta_t` and Gaussian kernel density estimation.

### Neuronal workflow

Each `entropy_neuronal_vs_diffusive_ensemble*.npz` contains:

| Key | Meaning / shape |
| --- | --- |
| `R_grid` | Five coupling values |
| `taus`, `taus_samp` | Integration and sampled times |
| `S_neur`, `S_diff` | Entropies with shape (5, 30, number of sampled times) |
| `dt`, `stride`, `dt_eff` | Integration step, sampling stride and sampled interval |
| `z`, `p_neg`, `seeds_A` | Sign-assignment ensemble size, negative fraction and seeds |
| `F`, `B`, `eps` | Model and spectral-floor parameters |
| `L_abs` | Unsigned combinatorial Laplacian used for diffusion |

The diffusion result is calculated once per R and broadcast over the 30 realizations. The two scripts use the same sign assignments and initial condition at different coupling values. Failed calculations can remain NaN because exceptions are caught in the scripts.

### Empirical and null-model workflows

`ER_CM_ensemble.npz` contains `A_er` and `A_cm` with intended shape (100, N, N), `seeds`, `p_er` and `swap_factor`. Each workflow has its own file with this name.

Empirical trajectory files contain `t`, `x` with shape (N, T), `A` and model parameters. Ensemble trajectory files contain `x` with shape (z, N, T), times, seeds and model parameters.

Empirical thermodynamic files contain `t`, `S`, `Sdot`, `E` and `div`. Successful-ensemble archives contain `good_idx`, `seeds`, `taus`, `taus_samp`, `stride`, arrays `entropy`, `entropy_rate`, `energy`, `divergence` with shape (z_good, sampled times), and numerical/model parameters. `z_good` can be smaller than 100.

Failure archives store failed indices, seeds, steps, times and reasons. Partial archives are checkpoints of successful calculations, not automatic restart files. Final plots use ensemble mean and sample standard deviation; the neuronal manuscript plot instead uses median and quantiles.
