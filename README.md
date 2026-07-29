# analysis_bci

Offline MATLAB analysis tools for the BCI-VR pipeline. Two sub-packages:

```
analysis_bci/
├── matlab_simulation/   chunk-by-chunk replay of the full ROS pipeline
└── analysis_gdf/        quick EEG visualisation helpers
```

---

## matlab_simulation

Full offline re-implementation of every ROS processing node (FBCSP, artifact detector, sLDA, integrator), validated against ROS at MAE < 1e-6. Ten entry points:

| Script | Input | Purpose |
|---|---|---|
| `main_simulate` | any paradigm GDF | Single GDF → full pipeline → per-trial P(c1) plot + `simulate_summary.mat` |
| `main_session_overview` | any paradigm GDF(s) | Trial accuracy, time metrics, sample accuracy, CSP/sLDA/ERD analysis, within-trial accuracy by fixed time bin, artifact rate vs outcome, per-class confusion matrix — 10 SVGs + `session_summary.mat`; console also reports ITR (bits/min) and a chance-level significance test per file/paradigm |
| `main_compare_sessions` (deprecated: `main_compare_paradigms`) | any paradigm GDF(s) | Loads `eval_single_*.mat` files recursively → per-paradigm aggregate table + bar charts |
| `main_trial_dynamics` | any paradigm GDF(s) | Within-session learning/fatigue + recording-order confound check: accuracy trend, 1st-vs-2nd-half, TTH trend (vs normalised trial position), plus accuracy vs chronological recording order across all files/paradigms — real GDF events only, no YAML needed — 4 figures + `trial_dynamics_summary.mat` |
| `main_threshold_sweep` | any paradigm GDF(s) | Integrator sensitivity: HIT-rate/TIMEOUT-rate heatmaps over a `(threshold_class1, threshold_class2)` grid, re-evaluating each trial's already-simulated buffer — deliberately not a classic ROC (two independent thresholds + TIMEOUT) — 1 image per paradigm + pooled `ALL` + per-file + `threshold_sweep_summary.mat` |
| `main_roc_analysis` | any paradigm GDF(s) | Classifier-level ROC/AUC, pre-integrator: pools the raw per-frame sLDA `P(class 1)` (fused `P` for Hybrid, treated as a hypothetical single classifier) and sweeps a scalar threshold — the genuine ROC/AUC counterpart to `main_threshold_sweep` — 1 image per paradigm + an overlay comparison (not pooled) + per-file + `roc_summary.mat` (keeps raw scores/labels for `main_group_analysis`) |
| `main_hybrid_advantage_probs` | **hybrid only** | Classifier-probability analysis: per-stream mean P(target), frame accuracy, CVSA-fusion advantage, cluster-based permutation test on the fusion effect — 6 figures + `hybrid_advantage_probs_summary.mat` |
| `main_hybrid_advantage_integ` | **hybrid only** | Counterfactual integrator analysis: runs the buffer 3× (hybrid / MI-only / CVSA-only) with the SAME fixed buffer/thresholds — per-trial + aggregate + temporal advantage + per-stream confusion matrix — 6 figures + `counterfactual_summary.mat`; console also reports ITR (bits/min) and a chance-level significance test per stream |
| `main_validate_counterfactual` | `session_summary.mat` + `counterfactual_summary.mat` | Compares real between-session accuracy vs offline simulation — validates the counterfactual methodology — 1 figure |
| `main_group_analysis` | root folder (recursive scan) | Aggregates `session_summary.mat`, `counterfactual_summary.mat`, `hybrid_advantage_probs_summary.mat`, `topo_erders_summary.mat`, and `roc_summary.mat` across **all subjects** found under the root → group-level accuracy (both TIMEOUT=fail and decided-trials-only conventions, plus false-positive rate), counterfactual-advantage, per-paradigm/per-stream best-worst breakdowns, CVSA-fusion-mechanism (+ cross-subject cluster-based permutation test and a Stouffer meta-analytic combination of each subject's own within-session significance test, sidestepping the group sign-flip test's small-cohort power floor), ERD/ERS-grounding, classifier-ROC, and cross-subject performance-correlate figures with subject-level statistics + paired Cohen's d alongside every group-level p-value — 10 figures + `group_summary.mat` |
| `main_browse_gdf` | any paradigm GDF | Interactive scrollable viewer: classifier probabilities + integrator signal |

> **Note**: `main_hybrid_advantage_probs` and `main_hybrid_advantage_integ` both require **hybrid GDF files**. `main_validate_counterfactual` requires the `.mat` outputs of `main_session_overview` and `main_hybrid_advantage_integ` for one subject. `main_group_analysis` requires the same outputs (plus `main_hybrid_advantage_probs`'s, `topo_erders`'s, and `main_roc_analysis`'s) but pooled across **multiple subjects** under a common root folder (`recordings/<subject>/...`) — run the per-file scripts first for each subject. `main_trial_dynamics` has no such dependency — it only needs the GDF files themselves.

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate                  % saves simulate_summary.mat
main_session_overview          % saves session_summary.mat
main_trial_dynamics            % any paradigm — no YAML needed, saves trial_dynamics_summary.mat
main_threshold_sweep           % any paradigm — saves threshold_sweep_summary.mat
main_roc_analysis              % any paradigm — saves roc_summary.mat
main_hybrid_advantage_probs    % hybrid GDFs only — saves hybrid_advantage_probs_summary.mat
main_hybrid_advantage_integ    % hybrid GDFs only — saves counterfactual_summary.mat
main_validate_counterfactual   % loads the two .mat files above, one subject
main_group_analysis            % no GUI: hardcoded default root, all subjects (or pass a subjects_filter cellstr)
main_browse_gdf

% Batch runners, hardcoded RECORDINGS_ROOT + SUBJECTS list at the top, no GUI prompts:
batch_run_subjects             % edit SUBJECTS = {'a1','a2',...}, then run — calls run_subject_analysis per (subject, day)
batch_group_analysis           % same SUBJECTS list — calls main_group_analysis + group_topo_erders restricted to it
```

`main_group_analysis(root_dir, subjects_filter, show_figures)` has no GUI folder picker: `root_dir` defaults to a hardcoded path at the top of the script (edit it directly, or pass one) and `subjects_filter` is an optional cellstr (e.g. `{'a1','a2','a3'}`) to restrict aggregation to those subjects; omit it to use everyone found under the root. `group_topo_erders` (see below) follows the same convention.

The companion YAML (rosparam dump saved by `bag_bci` alongside each GDF) is loaded automatically. For calibration recordings the YAML lives in `../parameters/` relative to the GDF; `load_params_yaml` searches there as a fallback. `main_trial_dynamics` is the only GDF-based script that does **not** need the YAML — it works purely from GDF event codes.

Every GDF-based script saves its figures as SVG (full-screen size) into a subfolder of `<gdf_dir>/analysis_results/` **named after the script itself**: `simulate/`, `session_overview/`, `trial_dynamics/`, `threshold_sweep/`, `roc_analysis/`, `hybrid_advantage_probs/`, `hybrid_advantage_integ/`. `main_session_overview` also saves `session_summary.mat` in `session_overview/` and produces ten SVGs: `overview_trial_accuracy.svg`, `overview_time_metrics.svg`, `overview_sample_accuracy.svg`, `overview_csp_importance.svg`, `overview_slda_importance.svg`, `overview_erd_csp_corr.svg`, `overview_roi_timecourse.svg`, `overview_timebin_accuracy.svg` (within-trial classifier accuracy by fixed time bin, with per-bin trial count and a low-N warning marker — distinct from `main_trial_dynamics`'s trial-order analysis), `overview_artifact_rate.svg` (artifact rate during CF, HIT vs MISS/TIMEOUT, per paradigm), and `overview_confusion_matrix.svg` (real outcome by target class, pooled per paradigm — flags per-class bias an aggregate accuracy would hide). Figs 4–7 are skipped if no companion YAML is available. The console additionally reports, per file and per paradigm TOTAL, ITR (bits/trial and bits/min, Wolpaw formula) and an exact one-sided binomial chance-level significance test (decided trials only, p0=1/n_cls). `main_trial_dynamics` saves `trial_dynamics_summary.mat` in `trial_dynamics/`. `main_threshold_sweep` saves `threshold_sweep_summary.mat` in `threshold_sweep/`. `main_roc_analysis` saves `roc_summary.mat` in `roc_analysis/`, including the raw pooled per-frame scores/labels behind each paradigm's curve (not just fpr/tpr/auc) so `main_group_analysis` can re-pool a subject's sessions. `main_hybrid_advantage_probs` saves `hybrid_advantage_probs_summary.mat` in `hybrid_advantage_probs/`. `main_hybrid_advantage_integ` produces a sixth figure, `confusion_matrix.svg` (per-stream target-class × outcome breakdown — since all three streams share the same buffer/thresholds, any per-class bias difference reflects the classifier signal, not the integration settings), reports ITR + chance-level per stream on the console, and saves `counterfactual_summary.mat` in `hybrid_advantage_integ/` (now including `itr_bits_trial`, `itr_bits_min`, `chance_p`, each `[1x3]`). `main_validate_counterfactual` saves `counterfactual_validation.svg` under a sibling `validate_counterfactual/` folder. `main_group_analysis` scans a root folder recursively (subject = first path component under the root) and saves `group_summary.mat` + ten SVGs under `<root>/group_analysis/`. Each GDF-based script has a `SHOW_FIGURES` flag at the top. `topo_erders.m` (in `analysis_gdf/`) saves into `<gdf_dir>/analysis_results/topo_erders/`, including `topo_erders_summary.mat` (per paradigm × band: mean ERD/ERS per class, discrimination, CSP-weight Pearson r) for cross-subject use by `main_group_analysis`.

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
| `group_topo_erders.m` | Cross-subject (group-average) ERD/ERS topoplots, reusing `topo_erders.m`'s per-subject output — see below |

### topo_erders.m

Loads one or more GDFs (any mix of MI/CVSA/Hybrid), runs EEGLAB preprocessing (resample to 128 Hz, 1-40 Hz bandpass, epoch around class onset, cue-only baseline removal `[-1500,0]` ms, CAR), and computes per-channel, per-trial ERD/ERS (Pfurtscheller-style band power, % change vs. pre-cue baseline) over the CSP bands actually used by that recording's paradigm (loaded from the companion rosparam YAML). Topoplots and the channel x time heatmaps (below) are both derived from this same per-trial band power, so they are always consistent with each other.

All figures are saved as SVG (full-screen size); the `SHOW_FIGURES` flag at the top of the script controls whether they are also displayed on screen.

For each paradigm found, output goes to `<gdf_dir>/analysis_results/topo_erders/<paradigm>/`:
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

**`topo_erders_summary.mat`**: saved directly under `analysis_results/topo_erders/` (one level above the per-paradigm subfolders). One row per (paradigm, band): `mean_erd_c1`, `mean_erd_c2` (averaged over that band's CSP-selected channels during the CF window), `discrimination` (`|mean_erd_c1 - mean_erd_c2|`), `csp_r` (the same Pearson r shown in the `corr_<paradigm>` panel for that band), `n_csp_channels`. Rows with no CSP info for that (origin, band) are skipped. Used by `main_group_analysis` to pool ERD/ERS discrimination and CSP-weight grounding across subjects.

**`topo_erders_channels.mat`**: also saved in the same folder — the full per-channel, per-time grand-average band power behind the topoplots (this subject's own average across their files of that paradigm/band), plus channel labels and CSP masks. Richer than the scalar-only `topo_erders_summary.mat`; this is what `group_topo_erders.m` needs to redraw a topoplot averaged across subjects. **Re-run `topo_erders.m` for any subject analysed before this file existed** — old runs only have `topo_erders_summary.mat`.

### group_topo_erders.m

Cross-subject (group-average) ERD/ERS topoplots. Run **after** `topo_erders.m` for every subject you want included. Recursively scans a root folder for every `topo_erders_channels.mat` (subject = first path component under the root, same convention as `main_group_analysis`) and reconstructs the same 7-column topoplot grid (Cue+CF, CF only, Cue, then four 1s CF bins) at the cohort level, drawn with `topo_map.m` (`matlab_simulation/utils/`) — **no EEGLAB dependency**, even though the script lives next to `topo_erders.m` in `analysis_gdf/` for symmetry. Channel labels must match across all subjects (fixed hardware montage) — errors clearly on a mismatch.

```matlab
group_topo_erders()                                  % default recordings root, ALL subjects
group_topo_erders(root_dir)                          % given root, ALL subjects
group_topo_erders(root_dir, {'a1','a2','a3'})        % only these subjects
group_topo_erders(root_dir, {'a1','a2','a3'}, true)  % also show figures on screen
group_topo_erders([], {'a1','a2'})                   % default root, subset of subjects
```

No GUI folder picker — `root_dir` defaults to a hardcoded path (`DEFAULT_ROOT` at the top of the script; edit it directly if your `recordings/` root differs) rather than prompting. The subject list is an explicit opt-in filter: pass a cellstr of subject IDs to average only those subjects; omit it (or pass `{}`) to use every subject found under the root.

If a subject contributes more than one file (e.g. several recording days), that subject's own entries are averaged first; the cross-subject average is a **simple (unweighted) mean over subjects** — every subject counts equally, matching how `main_group_analysis` treats most per-subject metrics. Both `all_<paradigm>_<origin>_<band>.svg` (every channel) and `sel_<paradigm>_<origin>_<band>.svg` (only channels selected by **any** subject's CSP for that origin — union across the cohort) are produced, same display convention per band-origin as `topo_erders.m` (MI = ERD only, CVSA = lateralization index; hybrid produces both, never mixed). Output under `<root>/group_topo_erders/<paradigm>/`. Console reports, per (paradigm, band-origin, frequency) group, how many of the requested/found subjects contributed and the CSP-mask footprint.

**Colour scale e significatività per canale**: i limiti colore di `all_*`/`sel_*` sono il 95° percentile di `|valore|` nella sola finestra CF (non il max/min grezzo) — un canale/istante rumoroso satura al colore più forte invece di schiacciare tutta la scala. Ogni cella (canale × intervallo × cue) passa anche un t-test a un campione cross-soggetto (H0: media=0); i canali con p<0.05 ricevono un anello nero sovrapposto al punto colorato — distingue un pattern medio genuino e consistente tra soggetti da un canale semplicemente rumoroso.

**Quick reference — cosa dovresti vedere in ogni figura di `group_topo_erders/<paradigma>/`:**

| Figura | Cosa rappresenta | Risultato sperato | Domanda a cui risponde |
|---|---|---|---|
| `all_<paradigma>_<origine>_<banda>.svg` | Topoplot ERD/ERS medio di gruppo (tutti i canali 10-20), 7 colonne temporali (Cue+CF, CF only, Cue, poi 4 bin CF da 1s); anello nero = canale con t-test cross-soggetto significativo (p<0.05) | Gli anelli neri si concentrano in zone neurofisiologicamente plausibili (es. C3/C4 per MI, parieto-occipitali per CVSA), non sparsi a caso sullo scalpo | La risposta ERD/ERS media del gruppo è genuina e consistente tra soggetti, o è rumore che si media via? |
| `sel_<paradigma>_<origine>_<banda>.svg` | Come `all_*` ma solo sui canali selezionati dal CSP di **almeno un** soggetto (unione cross-soggetto) | Il pattern ERD/ERS più marcato (e gli anelli di significatività) cade proprio sui canali che il CSP ha scelto | Il CSP ha selezionato canali che coincidono con dove il segnale neurofisiologico reale è più forte? |
| `heat_<paradigma>_<origine>_<banda>.svg` | Heatmap canale × tempo dell'ERD/ERS medio di gruppo (colormap divergente blu/bianco/rosso), una riga per classe, da cue onset a fine CF | Banda blu (ERD) che si intensifica dopo l'onset del CF nei canali attesi, evoluzione temporale liscia — non un mosaico di rumore senza struttura | Come evolve nel tempo il pattern ERD/ERS medio del gruppo, non solo la sua media su una finestra? |
| `corr_<paradigma>.svg` | Un pannello per banda con dati CSP: scatter canale-per-canale tra discriminazione ERD/ERS `\|c1−c2\|` durante CF e peso CSP di quel canale, r di Pearson in titolo | r positiva e "pulita" (pochi outlier che la trascinano) | A livello di gruppo, il CSP pesa di più i canali che discriminano di più dal punto di vista neurofisiologico? |
| `roi_<paradigma>_<banda>.svg` | Andamento temporale ERD/ERS medio di gruppo per ROI emisferica (Left/Right/Midline, canali CSP-selezionati raggruppati per lato), una riga per ROI: tracce per classe + lateralizzazione `c1(t)−c2(t)` | Lateralizzazione col segno atteso (es. MI mano destra → ERD più forte in C3, quindi `c1(t)-c2(t)` di segno coerente) stabile nel tempo durante il CF | La lateralizzazione emisferica attesa dalla neurofisiologia del compito è visibile nella media di gruppo? |

**Nota**: a differenza delle figure di `main_group_analysis`, qui non c'è un test statistico di gruppo esplicito riassunto in una singola stella/p — il segnale di "affidabilità" da guardare è l'anello nero per-canale su `all_*`/`sel_*` (t-test cross-soggetto) e la pulizia (assenza di outlier dominanti) dello scatter in `corr_*`.

---

## validate_features (in slda_bci)

`src/slda_bci/create_slda/validate_features.m` is a stage-by-stage comparison of the Python calibration notebook (`create_slda.ipynb`) vs `apply_processing.m`. It loads the `_training_features.mat` saved by Cell 10 of the notebook and compares pre-CSP, post-CSP, and log-feature power at aligned windows.

Expected residual MAE: Stage 1 ~5e-02, Stage 3 ~1.4e-02 — from the inherent 0–12 sample chunk-alignment offset between Python's onset-aligned windows and MATLAB's chunk-aligned ring buffer. This is not a bug; it reflects that `round(onset / chunk_size)` in MATLAB can differ by up to `chunk_size/2` samples from the exact onset position used in Python.
