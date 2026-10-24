# Reproduction guide

## Scope and working directory

This guide describes how to reproduce the panels that compose the two figures of **Information-thermodynamic principles unveil the response of complex networks to perturbations**.
Start Jupyter in the relevant workflow folder and execute notebooks in the order below, with cells run from top to bottom. Use a fresh kernel for each notebook. Keep each folder's generated files separate: the three empirical workflows use the same filename `ER_CM_ensemble.npz` for different networks.


## Fig. 1a: synthetic networks

Working directory: `synth_biochemical/`.

1. Run `python synthgen_ensemble.py`. This creates `network_ensembles.npz`.
2. Execute `timediff.ipynb`. It computes biochemical trajectories and entropy-rate peak-time differences relative to diffusion.
3. Execute `plot_distr.ipynb` to plot the distributions of the absolute peak-time difference.

The generator uses N = 100, target average degree 8, and 100 realizations per model, with seeds 123–222. It generates ER, WS, BA, SBM and random geometric networks; the supplied analysis and manuscript panel use ER, BA and SBM.

Biochemical parameters are F = 0.5, B = 0.1 and R = 1, with a final time of 5 and 500 time points. The initial condition uses seed 42. The diffusion comparison uses the symmetric normalized Laplacian.

Intermediate outputs are `traj_bio_er/`, `traj_bio_ba/`, `traj_bio_sbm/` and `delta_t_bio_{er,ba,sbm}.npz`. The principal plot is `fig_kde_delta_t_bio_compare2.pdf`, with separate `legend_ER.pdf`, `legend_BA.pdf` and `legend_SBM.pdf` files.

## Fig. 1b: C. elegans

Working directory: `Neuronal/`. The supplied `neuronal.npy` is a 279 × 279 adjacency matrix. Two independent scripts split the calculation into five coupling values each, so a smaller run can be tested without computing all ten points.

| Script | Coupling values R | Output archive |
| --- | --- | --- |
| `dataneuro3.py` | 0.10, 0.20, 0.30, 0.40, 0.50 | `entropy_neuronal_vs_diffusive_ensemble.npz` |
| `dataneuro3_1.py` | 0.05, 0.15, 0.25, 0.35, 0.45 | `entropy_neuronal_vs_diffusive_ensemble1.npz` |

Both archives contain neuronal and diffusive entropy data. The scripts use identical dynamics, numerical settings and random seeds; only the coupling grids and output filenames differ in executable code.

### Full figure: both runs are required by the default visualizer

From `Neuronal/`, execute:

```bash
python dataneuro3.py
python dataneuro3_1.py
```

Each script also displays diagnostic plots. Close its plot windows if needed to let the process finish before starting the next command. The execution order of the two scripts is otherwise irrelevant.

Once both output archives exist, execute `plots3.ipynb`. Its default input list is:

```python
files = [
    "entropy_neuronal_vs_diffusive_ensemble.npz",
    "entropy_neuronal_vs_diffusive_ensemble1.npz",
]
```

For each R, the notebook calculates the total duration of positive sampled entropy increments for each valid realization, then the median and 16th–84th percentiles. It concatenates results from the two runs and sorts them by R, producing ten points from 0.05 to 0.50. It does not combine the runs into a 60-realization ensemble at each coupling. The default notebook fails with a missing-file error if only one archive is present.

The main output is `figs_entropy/timepos_entropy_positive_median_quantiles.pdf`. The same directory receives `timepos_entropy_positive_median_quantiles_WITH_legend_and_title2.pdf` and `legend_timepos_entropy_positive_median_quantiles2.pdf`.

### Testing and visualizing only one run

No simulation-code changes are needed: run only the script for the desired grid. In `plots3.ipynb`, replace the two-item `files` list with exactly one of the following.

After running `dataneuro3.py`:

```python
files = ["entropy_neuronal_vs_diffusive_ensemble.npz"]
```

After running `dataneuro3_1.py`:

```python
files = ["entropy_neuronal_vs_diffusive_ensemble1.npz"]
```

Then restart the notebook kernel and run all cells. The existing concatenation and sorting logic also works with a one-item list, and the figure contains five R values. No other plotting changes are required. Restore the two-item list to produce the full ten-point figure. Both modes use the same output filenames, so save any five-point plots separately if you want to retain them.

Apply the same `files` change if using `entropyplot.ipynb`. The latter plots entropy trajectories and checks that the input time grids agree when combining runs.

`vizneuro.ipynb` is a separate diagnostic visualizer: it already reads only `entropy_neuronal_vs_diffusive_ensemble1.npz` through its `NPZ_FILE` variable. It therefore needs only `dataneuro3_1.py`; to inspect the other run, set:

```python
NPZ_FILE = "entropy_neuronal_vs_diffusive_ensemble.npz"
```

Its mean/std plots differ from the median/quantile statistic used in Fig. 1b. The simulation scripts also produce their own diagnostic plots from a single run.

### Shared settings

Both scripts use F = 0, B = 1, 30 sign assignments with seeds 12345–12374, and a negative-edge fraction of 0.26. The final time is 15, with 1500 time points and entropy sampling stride 15. The initial condition uses seed 1. Diffusion is computed with J = -R L(|A|), using the combinatorial Laplacian of the unsigned adjacency matrix.

The manuscript plotting notebook counts sampled entropy increments above 1e-10 and multiplies their count by the sampling interval.
## Fig. 2: empirical networks and null models

| Folder | Input adjacency | Size | Dynamics notebook | Plot tag |
| --- | --- | --- | --- | --- |
| `Noisy/` | `noisy.npy` | 62 × 62 | `dataNoisy.ipynb` | `noisy` |
| `Epidemics/` | `epidemic.npy` | 95 × 95 | `dataEpi.ipynb` | `epi` |
| `Population/` | `pop.npy` | 183 × 183 | `dataPop.ipynb` | `pop` |

In each folder:

1. Execute `ERandCM.ipynb` to create `ER_CM_ensemble.npz`, containing `seeds`, `A_er`, `A_cm`, `p_er` and `swap_factor`. It generates 100 ER and 100 rewired CM realizations, using seeds 100–199. **This step imports `graph_tool.all`.**
2. Confirm that `ER_CM_ensemble.npz` is in the same folder as the empirical `.npy` input. The dynamics notebook reads the empirical matrix directly from that `.npy` file.
3. Execute the dynamics notebook to generate trajectories and thermodynamic observables.
4. Execute `plots2.ipynb` to export `fig_energy_and_free_energy_diff_<tag>_compare.pdf` and `legend_<tag>.pdf`. Additional entropy plots are also exported.

### Input availability

All four empirical `.npy` matrices are included. `ER_CM_ensemble.npz`, trajectories, thermodynamic archives and PDFs are generated by the workflows and do not need to be shipped as inputs.

### Parameters and outputs

| Dynamics | Parameters | Final time | Time points |
| --- | --- | --- | --- |
| Noisy voter | F = 0.2, B = 2, R = 1.8 | 10 | 2000 |
| SIS | B = 5, R = 1 | 0.5 | 500 |
| Population | a = 1, b = 2, B = 1, R = 1 | 4 | 500 |

All three use initial-condition seed 1 and ensemble thermodynamic sampling stride 5. The trajectory files use `dtraj_<tag>_d.npz` and `dtraj_<tag>_{er,cm}_ensemble.npz`. Empirical thermodynamics are saved as `dthermo_<tag>_d.npz`. Ensemble results are saved as:

- Noisy: `thermo_noisy_{er,cm}_good_final.npz`.
- SIS: `thermo_sis_epi_{er,cm}_good_final.npz`.
- Population: `thermo_pop_{er,cm}_good_final.npz`.

The ensemble notebooks also write partial progress and failure records. Check successful-realization counts before interpreting a nominal 100-realization ensemble. The neuronal script catches exceptions and can retain NaNs, so successful script completion alone does not guarantee complete results.

## Observable names and figure assembly

In thermodynamic output files, `energy` or `E` denotes the manuscript's informational potential U(t) = -2 Tr(rho D). `divergence` or `div` denotes the relative entropy to the instantaneous Gibbs reference. The plotting notebooks compute the free-energy difference as this divergence divided by time. The internal matrix variable `U` used during integration is the statistical propagator, not the scalar potential with the same letter in the paper.

The plotting notebooks export individual panels and legends. No script assembling the exact final manuscript layout was identified. For Fig. 2, arrange the outputs in manuscript order: Noisy, Epidemics, Population. Arial is requested with Helvetica/DejaVu Sans fallbacks; available fonts may change appearance.

Dataset identification follows from the manuscript.

## Troubleshooting

- `ModuleNotFoundError: graph_tool`: install graph-tool in the environment used by the Jupyter kernel, then restart the kernel. Installing it for a different Python interpreter is insufficient.
- Missing `ER_CM_ensemble.npz`: run `ERandCM.ipynb` first in that same folder.
- Missing trajectory or thermodynamic archive: execute the dynamics notebook from the beginning; plot notebooks do not generate simulation results.
- Missing neuronal entropy archive: run both neuronal scripts or use the single-run input list described above.
- Fewer successful realizations than expected: inspect `*_failures.npz`, the printed good-run counts and NaNs. A completed notebook is not sufficient evidence that all realizations succeeded.

