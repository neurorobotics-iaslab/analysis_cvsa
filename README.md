# analysis_bci

Offline MATLAB analysis tools for the BCI-VR pipeline. Two sub-packages:

```
analysis_bci/
├── matlab_simulation/   chunk-by-chunk replay of the full ROS pipeline
└── analysis_gdf/        quick EEG visualisation helpers
```

---

## matlab_simulation

Full offline re-implementation of every ROS processing node (FBCSP, artifact detector, sLDA, integrator), validated against ROS at MAE < 1e-6. Three entry points:

| Script | Purpose |
|---|---|
| `main_simulate` | Single GDF → full pipeline → per-trial P(c1) plot |
| `main_batch_evaluate` | Multi-file metrics (accuracy, time-to-hit, confidence, artifact rate, …) + optional per-class trial figures. Saves `eval_<paradigm>_<basename>.mat`. |
| `main_compare_paradigms` | Loads three eval `.mat` files (MI / CVSA / Hybrid) → side-by-side table + bar charts + precision/recall scatter |

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate
main_batch_evaluate
main_compare_paradigms
```

The companion YAML (rosparam dump saved by `bag_bci` alongside each GDF) is loaded automatically. For calibration recordings the YAML lives in `../parameters/` relative to the GDF; `load_params_yaml` searches there as a fallback.

Requirements: MATLAB R2019+ with Signal Processing Toolbox, [yamlmatlab](https://github.com/jiri-cigler/yamlmatlab), [BIOSIG](https://biosig.sourceforge.io) (`sload`).

See **`matlab_simulation/README.md`** for the full pipeline description, stage details, and variable glossary.

---

## analysis_gdf

Lightweight MATLAB scripts for quick inspection of raw GDF recordings, independent of the simulation pipeline:

| Script | Purpose |
|---|---|
| `check_psd.m` | Power spectral density per channel and paradigm band; useful for sanity-checking the EEG quality before training |
| `eeglab_gdf.m` | Loads a GDF into EEGLAB for visual inspection (requires EEGLAB on the MATLAB path) |

---

## validate_features (in slda_bci)

`src/slda_bci/create_slda/validate_features.m` is a stage-by-stage comparison of the Python calibration notebook (`create_slda.ipynb`) vs `apply_processing.m`. It loads the `_training_features.mat` saved by Cell 10 of the notebook and compares pre-CSP, post-CSP, and log-feature power at aligned windows.

Expected residual MAE: Stage 1 ~5e-02, Stage 3 ~1.4e-02 — from the inherent 0–12 sample chunk-alignment offset between Python's onset-aligned windows and MATLAB's chunk-aligned ring buffer. This is not a bug; it reflects that `round(onset / chunk_size)` in MATLAB can differ by up to `chunk_size/2` samples from the exact onset position used in Python.
