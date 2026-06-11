# analysis_bci

Offline MATLAB analysis tools for the BCI-VR pipeline. Two sub-packages:

```
analysis_bci/
├── matlab_simulation/   chunk-by-chunk replay of the full ROS pipeline
└── analysis_gdf/        quick EEG visualisation helpers
```

---

## matlab_simulation

Full offline re-implementation of every ROS processing node (FBCSP, artifact detector, sLDA, integrator), validated against ROS at MAE < 1e-6. Five entry points:

| Script | Purpose |
|---|---|
| `main_simulate` | Single GDF → full pipeline → per-trial P(c1) plot |
| `main_session_overview` | One or more GDFs → metrics + ERD/ERS heatmaps + topoplots. Saves SVG and `eval_single_<paradigm>_<basename>.mat`. |
| `main_compare_sessions` (deprecated name: `main_compare_paradigms`) | Loads `eval_single_*.mat` files recursively → per-paradigm aggregate table + bar charts |
| `main_hybrid_advantage` | Matched hybrid-vs-MI-vs-CVSA comparison + saved-trial breakdown |
| `main_browse_gdf` | Interactive scrollable viewer: classifier probabilities + integrator signal |

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate
main_session_overview
main_hybrid_advantage
main_browse_gdf
```

The companion YAML (rosparam dump saved by `bag_bci` alongside each GDF) is loaded automatically. For calibration recordings the YAML lives in `../parameters/` relative to the GDF; `load_params_yaml` searches there as a fallback.

`main_simulate`, `main_session_overview`, and `main_hybrid_advantage` save their figures as PNG into an `analysis_results/` subfolder created next to the selected GDF(s) (i.e. `<gdf_dir>/analysis_results/`).

Requirements: MATLAB R2019+ with Signal Processing Toolbox, [yamlmatlab](https://github.com/jiri-cigler/yamlmatlab), [BIOSIG](https://biosig.sourceforge.io) (`sload`).

See **`matlab_simulation/README.md`** for the full pipeline description, stage details, and variable glossary.

---

## analysis_gdf

Lightweight MATLAB scripts for quick inspection of raw GDF recordings, independent of the simulation pipeline:

| Script | Purpose |
|---|---|
| `check_psd.m` | Power spectral density per channel and paradigm band; useful for sanity-checking the EEG quality before training |
| `eeglab_gdf.m` | Loads a GDF into EEGLAB for visual inspection (requires EEGLAB on the MATLAB path) |
| `topo_erders.m` | ERD/ERS topoplot grids over time (EEGLAB `newtimef`), per paradigm and per file — see below |

### topo_erders.m

Loads one or more GDFs (any mix of MI/CVSA/Hybrid), runs EEGLAB preprocessing (resample to 128 Hz, 1-40 Hz bandpass, epoch around class onset, cue-only baseline removal `[-1500,0]` ms, CAR), and computes per-channel/per-class ERD/ERS (`newtimef`) over the CSP bands actually used by that recording's paradigm (loaded from the companion rosparam YAML).

For each paradigm found, output goes to `results_eeglab/<paradigm>/`:
- `all_<paradigm>` — average across all files of that paradigm, all channels
- `sel_<paradigm>` — same, but only CSP-selected channels are non-zero (MI bands → MI channels, CVSA bands → CVSA channels; hybrid keeps both separate, never mixed)
- `per_file/all_<filename>` and `per_file/sel_<filename>` — same pair for each individual file

Display convention per band origin:
- **MI bands**: ERD only — values above 0 are clipped to 0, one row per cue (e.g. "Cue 769", "Cue 770")
- **CVSA bands**: lateralization index — single row showing `P(cue1) - P(cue2)` instead of separate class rows
- **Hybrid**: both representations are produced (one for MI-origin bands, one for CVSA-origin bands)

**CAR**: implemented manually (not `pop_reref`) to match `rosneuro_filters_car/Car.hpp` exactly — the common average is computed over all channels *except* `Fp1`/`Fp2` (the EOG channels excluded in `car.yaml`), but is then subtracted from **every** channel, including `Fp1`/`Fp2`. EEGLAB's `pop_reref(..., 'exclude', ...)` does not match this (excluded channels are left completely untouched), hence the manual implementation.

---

## validate_features (in slda_bci)

`src/slda_bci/create_slda/validate_features.m` is a stage-by-stage comparison of the Python calibration notebook (`create_slda.ipynb`) vs `apply_processing.m`. It loads the `_training_features.mat` saved by Cell 10 of the notebook and compares pre-CSP, post-CSP, and log-feature power at aligned windows.

Expected residual MAE: Stage 1 ~5e-02, Stage 3 ~1.4e-02 — from the inherent 0–12 sample chunk-alignment offset between Python's onset-aligned windows and MATLAB's chunk-aligned ring buffer. This is not a bug; it reflects that `round(onset / chunk_size)` in MATLAB can differ by up to `chunk_size/2` samples from the exact onset position used in Python.
