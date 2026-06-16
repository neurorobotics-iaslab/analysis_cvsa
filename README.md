# analysis_bci

Offline MATLAB analysis tools for the BCI-VR pipeline. Two sub-packages:

```
analysis_bci/
├── matlab_simulation/   chunk-by-chunk replay of the full ROS pipeline
└── analysis_gdf/        quick EEG visualisation helpers
```

---

## matlab_simulation

Full offline re-implementation of every ROS processing node (FBCSP, artifact detector, sLDA, integrator), validated against ROS at MAE < 1e-6. Six entry points:

| Script | Input GDFs | Purpose |
|---|---|---|
| `main_simulate` | any paradigm | Single GDF → full pipeline → per-trial P(c1) plot |
| `main_session_overview` | any paradigm | One or more GDFs → trial accuracy, time metrics, sample accuracy, **CSP channel/band importance** (Fig 4), **sLDA-weighted feature importance** (Fig 5), **ERD/ERS ↔ CSP correlation** (Fig 6), and **hemisphere ROI time-courses** (Fig 7). Saves seven SVGs. |
| `main_compare_sessions` (deprecated name: `main_compare_paradigms`) | any paradigm | Loads `eval_single_*.mat` files recursively → per-paradigm aggregate table + bar charts |
| `main_hybrid_advantage_probs` | **hybrid only** | Classifier-probability analysis: per-stream mean P(target), frame accuracy, CVSA-fusion advantage over the leaky-integrator buffer — 5 figures |
| `main_hybrid_advantage_integ` | **hybrid only** | Counterfactual integrator analysis: runs the buffer 3× (hybrid / MI-only / CVSA-only) on the same sLDA streams — per-trial + aggregate + temporal advantage + deep analysis — 5 figures |
| `main_browse_gdf` | any paradigm | Interactive scrollable viewer: classifier probabilities + integrator signal |

> **Note**: `main_hybrid_advantage_probs` and `main_hybrid_advantage_integ` both require **hybrid GDF files** (recorded during a hybrid paradigm session). They will not produce meaningful results with MI-only or CVSA-only recordings.

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate
main_session_overview
main_hybrid_advantage_probs    % hybrid GDFs only
main_hybrid_advantage_integ    % hybrid GDFs only
main_browse_gdf
```

The companion YAML (rosparam dump saved by `bag_bci` alongside each GDF) is loaded automatically. For calibration recordings the YAML lives in `../parameters/` relative to the GDF; `load_params_yaml` searches there as a fallback.

`main_simulate`, `main_session_overview`, `main_hybrid_advantage_probs`, and `main_hybrid_advantage_integ` all save their figures as SVG (full-screen size) into a dedicated subfolder of `<gdf_dir>/analysis_results/`: `trial_simulation/`, `overview/`, `advantage_hybrid/`, and `hybrid_vs_unimodal/` respectively. `main_session_overview` produces seven SVGs in `overview/`: `overview_trial_accuracy.svg`, `overview_time_metrics.svg`, `overview_sample_accuracy.svg`, `overview_csp_importance.svg` (CSP filter weights — per-class channel topoplots with full 10-20 grid, selectivity, stacked band bars), `overview_slda_importance.svg` (sLDA-weighted importance — channel topoplot weighted by |sLDA coef|, band bars, feature selection heatmap comp × band), `overview_erd_csp_corr.svg` (ERD/ERS class discrimination vs CSP filter weight scatter, Pearson r per band), and `overview_roi_timecourse.svg` (Left/Right/Midline hemisphere ROI ERD/ERS time-courses, band-averaged, per class + lateralization). Figs 4–7 are skipped if no companion YAML is available. Each script has a `SHOW_FIGURES` flag at the top controlling whether the figures are also displayed on screen. `topo_erders.m` (in `analysis_gdf/`) follows the same pattern, saving into `<gdf_dir>/analysis_results/eeglab/`.

Requirements: MATLAB R2019+ with Signal Processing Toolbox, [yamlmatlab](https://github.com/jiri-cigler/yamlmatlab), [BIOSIG](https://biosig.sourceforge.io) (`sload`).

See **`matlab_simulation/README.md`** for the full pipeline description, stage details, and variable glossary.

---

## analysis_gdf

Lightweight MATLAB scripts for quick inspection of raw GDF recordings, independent of the simulation pipeline:

| Script | Purpose |
|---|---|
| `check_psd.m` | Power spectral density per channel and paradigm band; useful for sanity-checking the EEG quality before training |
| `eeglab_gdf.m` | Loads a GDF into EEGLAB for visual inspection (requires EEGLAB on the MATLAB path) |
| `topo_erders.m` | ERD/ERS topoplot grids + channel x time heatmaps + CSP-weight correlation, per paradigm and per file — see below |

### topo_erders.m

Loads one or more GDFs (any mix of MI/CVSA/Hybrid), runs EEGLAB preprocessing (resample to 128 Hz, 1-40 Hz bandpass, epoch around class onset, cue-only baseline removal `[-1500,0]` ms, CAR), and computes per-channel, per-trial ERD/ERS (Pfurtscheller-style band power, % change vs. pre-cue baseline) over the CSP bands actually used by that recording's paradigm (loaded from the companion rosparam YAML). Topoplots and the channel x time heatmaps (below) are both derived from this same per-trial band power, so they are always consistent with each other.

All figures are saved as SVG (full-screen size); the `SHOW_FIGURES` flag at the top of the script controls whether they are also displayed on screen.

For each paradigm found, output goes to `<gdf_dir>/analysis_results/eeglab/<paradigm>/`:
- `all_<paradigm>` — average across all files of that paradigm, all channels
- `sel_<paradigm>` — same, but only CSP-selected channels are non-zero (MI bands → MI channels, CVSA bands → CVSA channels; hybrid keeps both separate, never mixed)
- `per_file/all_<filename>` and `per_file/sel_<filename>` — same pair for each individual file

Display convention per band origin:
- **MI bands**: ERD only — values above 0 are clipped to 0, one row per cue (e.g. "Cue 769", "Cue 770")
- **CVSA bands**: lateralization index — single row showing `P(cue1) - P(cue2)` instead of separate class rows
- **Hybrid**: both representations are produced (one for MI-origin bands, one for CVSA-origin bands)

**CAR**: implemented manually (not `pop_reref`) to match `rosneuro_filters_car/Car.hpp` exactly — the common average is computed over all channels *except* `Fp1`/`Fp2` (the EOG channels excluded in `car.yaml`), but is then subtracted from **every** channel, including `Fp1`/`Fp2`. EEGLAB's `pop_reref(..., 'exclude', ...)` does not match this (excluded channels are left completely untouched), hence the manual implementation.

**Channel x time ERD/ERS heatmaps** (`heat_*` files, same folders as the topoplots): for each band, one figure with 2 rows (one per cue/class — never the lateralization-index collapsing used by the topoplots), each a channel-by-time heatmap from cue onset (t=0) to the end of the longest possible CF window, diverging blue(ERD)/white(0)/red(ERS) colormap with symmetric limits = `max(abs(value))`. `all_*` shows every channel; `sel_*` shows only CSP-selected channels (MI/CVSA: that paradigm's own CSP channels; Hybrid: union of MI+CVSA CSP channels, for every band).

Per-trial values are computed as Pfurtscheller-style band power (`bandpass` → instantaneous power → 200 ms moving average) as a percentage of that trial's own pre-cue baseline. Because evaluation trials have variable CF duration, each trial's samples after its own outcome event (897/898/899 latency) are set to NaN before averaging, so each time point in the heatmap is the mean over only the trials whose CF is still running at that time; NaN (no-data) regions render as the figure background.

**Topoplots use the same NaN-masking**: each `topo_*` time-interval column — "Cue+CF (0-5s)", "CF only (1.5-5s)", "Cue (0-1.5s)", then 1s-wide CF bins "CF (1.5-2.5s)", "CF (2.5-3.5s)", "CF (3.5-4.5s)", "CF (4.5-5s)" (last bin shorter if `epoch_limits(2) - CUE_DURATION_S` isn't a whole number) — is the mean of this same per-trial band power over that interval, omitting NaNs, so e.g. the "CF (4.5-5s)" column only averages over trials whose CF is still running at 4.5-5s. If NO trial of a class reaches a given interval, that channel/interval is shown as 0 (no ERD/ERS). Colour limits are computed from the same data, restricted to the CF window only `[CUE_DURATION_S, epoch_limits(2)]` (CF onset → max CF end) — i.e. every topoplot column is colour-mapped to the ERD/ERS range observed during continuous feedback, regardless of which interval it shows.

**ERD/ERS vs. CSP weight** (`corr_<paradigm>` files): one figure per paradigm, one panel per band that has CSP info — scatter of each channel's class discrimination `|ERD/ERS(class1) - ERD/ERS(class2)|` during continuous feedback (from `CUE_DURATION_S` to the max CF end) against that channel's CSP weight (sum of `|csp_matrices|` over components for that band; 0 for channels not selected by the CSP). The `|c1-c2|` difference (not the per-class average) is used for both MI and CVSA bands: it captures contralateral ERD/ERS patterns (e.g. C4 ERD for class1 / C3 ERD for class2 both appear as large `|c1-c2|`, whereas averaging the classes would cancel them out), and matches the lateralization index already used for CVSA topoplots. MI-origin bands use the MI CSP, CVSA-origin bands use the CVSA CSP (hybrid produces one panel per band, using the matching CSP). Each panel's title shows the Pearson correlation coefficient between this discrimination signal and CSP weight.

**ERD/ERS lateralization time-course by hemisphere ROI** (`roi_<paradigm>_<band>` and `per_file/roi_<filename>_<band>` files): one figure per band that has CSP info — channels selected by that band's CSP (MI-origin bands → MI CSP, CVSA-origin bands → CVSA CSP) are grouped into hemisphere ROIs (Left/Right/Midline, by the trailing digit/`z` of the 10-20 label: odd = left, even = right, `z` = midline; channels matching neither are dropped). One row per non-empty ROI, from cue onset (t=0) to the max CF end, dashed line = CF onset: left column overlays the per-class ERD/ERS time-courses (`c1(t)`, `c2(t)`); right column shows the lateralization time-course `c1(t) - c2(t)`. Bands without CSP info, or whose origin has no CSP-selected channels, are skipped.

---

## validate_features (in slda_bci)

`src/slda_bci/create_slda/validate_features.m` is a stage-by-stage comparison of the Python calibration notebook (`create_slda.ipynb`) vs `apply_processing.m`. It loads the `_training_features.mat` saved by Cell 10 of the notebook and compares pre-CSP, post-CSP, and log-feature power at aligned windows.

Expected residual MAE: Stage 1 ~5e-02, Stage 3 ~1.4e-02 — from the inherent 0–12 sample chunk-alignment offset between Python's onset-aligned windows and MATLAB's chunk-aligned ring buffer. This is not a bug; it reflects that `round(onset / chunk_size)` in MATLAB can differ by up to `chunk_size/2` samples from the exact onset position used in Python.
