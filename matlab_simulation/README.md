# matlab_simulation — Offline MATLAB replay of the ROS BCI pipeline

This package replays a **GDF recording** chunk-by-chunk through a pure-MATLAB
re-implementation of every ROS node in the live pipeline
(`artifacts_bci`, `processing_bci` FBCSP, `slda_bci`, `rosneuro_integrator`),
and produces a per-trial plot of the integrated signal during continuous
feedback. The final test you'll typically look at is: **does the integrated
P(class 1) reach the launch-file threshold inside the CF window?**

Everything is function-based (no MATLAB classes), single-file-per-stage,
verbose at every step (`log_step` prefixes every line with `[sim]`).

The core stage functions (`apply_processing`, `detect_artifacts`,
`apply_slda`, the integrator step, the per-class normalisation) are
**ported directly from the validated reference tests** in the workspace:

| Stage | Reference test (validated against ROS) | MAE vs ROS |
|---|---|---|
| Artifact detector  | [`src/artifacts_bci/test/test_artifacts.m`](../../artifacts_bci/test/test_artifacts.m) | ~ 1e-7 |
| FBCSP (CAR + BP + ringbuf + CSP + power) | [`src/slda_bci/test/test_slda.m`](../../slda_bci/test/test_slda.m) (GDF mode) | < 1e-6 |
| sLDA (log + sigmoid) | [`src/slda_bci/test/test_slda.m`](../../slda_bci/test/test_slda.m) | < 1e-6 |
| Integrator (leaky WTA + fusion + normalize) | [`src/test_pipeline/src/full_pipeline.m`](../../test_pipeline/src/full_pipeline.m) | < 1e-6 |
| Per-class normalisation | [`src/feedback_bci_vr/src/Training.cpp`](../../feedback_bci_vr/src/Training.cpp) `normalize_input` | exact port |

`apply_processing.m` is also used by [`src/slda_bci/create_slda/validate_features.m`](../../slda_bci/create_slda/validate_features.m), which runs a stage-by-stage comparison of the calibration notebook (Python) vs the MATLAB implementation. Because `apply_processing` has been validated against ROS at ~1e-7 MAE, any discrepancy found by `validate_features` is attributed to the Python side.

So this simulator is meant to behave like ROS to numerical noise, not just
in spirit.

---

## 1. How to run

### Recommended: batch runner

`run_subject_analysis` is the standard entry point for running all analyses on a subject.
Recordings are expected in the layout:

```
recordings/<subject>/<day>/evaluation/*.gdf   ← steps 1-9 (including topo_erders)
recordings/<subject>/<day>/calibration/*.gdf  ← topo_erders only
```

Subject folders must start with **`a`**; other names are ignored.

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
run_subject_analysis       % opens uigetdir — pick any of:
```

| Folder selected | Behaviour |
|---|---|
| `recordings/` | Finds all `a*` subjects, processes **all** their days in batch |
| `recordings/a001/` | Processes that subject; if multiple days → `listdlg` to pick |
| `recordings/a001/20250101/` | Finds `evaluation/` inside, uses it |
| `recordings/a001/20250101/evaluation/` | Used directly |

GDFs are read **only** from the immediate folder — sub-folders such as `test/`, `game/` are ignored automatically.

Steps run in order for each `evaluation/` folder:
1. `main_simulate` — all GDFs
2. `main_session_overview` — all GDFs
3. `main_trial_dynamics` — all GDFs
4. `main_threshold_sweep` — all GDFs
5. `main_roc_analysis` — all GDFs: classifier-level ROC/AUC (pre-integrator); saves raw pooled scores/labels per paradigm for `main_group_analysis`
6. `main_hybrid_advantage_probs` — hybrid GDFs only (skipped if none)
7. `main_hybrid_advantage_integ` — hybrid GDFs only (skipped if none)
8. `main_validate_counterfactual` — auto-skipped if step 2 or 7 `.mat` outputs are missing
9. `topo_erders` (from `analysis_gdf/`) — ERD/ERS topoplots (requires EEGLAB)

Then, if a sibling `calibration/` folder exists and contains GDFs:

- `topo_erders` — same script, run on calibration GDFs

After processing all subjects, run the group aggregator:

```matlab
main_group_analysis   % uigetdir → select recordings/
```

### Individual scripts (interactive or programmatic)

All `main_*.m` scripts are **MATLAB functions**, so they can be called either
interactively (no arguments → GUI picker) or programmatically (pass paths directly):

```matlab
% Interactive — shows GUI file picker as before
main_session_overview

% Programmatic — pass directory and file list
main_session_overview('/path/to/evaluation', {'file1.gdf', 'file2.gdf'})
main_hybrid_advantage_integ('/path/to/evaluation', {'subj.hybrid.gdf'})

% validate_counterfactual takes two .mat paths directly
main_validate_counterfactual('/path/.../session_summary.mat', '/path/.../counterfactual_summary.mat')
```

`show_figures` (default `false`) can be passed as an optional third argument to display figures on screen in addition to saving them.

---

Eleven individual analysis scripts + three batch runners (`run_subject_analysis`, `batch_run_subjects`, `batch_group_analysis`). `main_hybrid_advantage_probs` and `main_hybrid_advantage_integ` **require hybrid GDF files**. `main_trial_dynamics` works on any paradigm and needs no companion YAML. `main_threshold_sweep` and `main_roc_analysis` work on any paradigm (or a mix) and need the companion YAML, since they reuse the full simulated pipeline.

### `main_simulate` — single-file, visual inspection

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate
```

Single GDF → companion YAML → full pipeline → one figure with all CF trials.
Use this to inspect a specific recording in detail. The figure is saved as SVG
(full-screen size) under `<gdf_dir>/analysis_results/simulate/`; the
`SHOW_FIGURES` flag at the top of the script controls whether it is also shown
on screen. Also saves `simulate_summary.mat` in that folder (per-file: basename, paradigm, n_trials, n_hit, n_miss, n_to, n_pass_sim) for cross-subject use by `main_group_analysis`.

### `main_session_overview` — multi-file cross-paradigm summary

```matlab
main_session_overview
```

Multi-select GDFs of any paradigm (MI, CVSA, Hybrid mixed). Paradigm is inferred
from the filename (`hybrid` > `cvsa` > `mi`). Files are sorted MI → CVSA → Hybrid
with group separators in the figures.

**Ten figures:**

| Figure | Content | Source |
|---|---|---|
| Fig 1 | Trial accuracy — two subplots: `897/(897+898+899)` and `897/(897+898)` (no timeout). Per-file bar + per-paradigm dashed mean line. | Real GDF events |
| Fig 2 | Time metrics — TTH (HIT trials, mean ± std + dots) and time-to-MISS (MISS trials). | Real GDF events |
| Fig 3 | Sample accuracy (3 subplots): MI classifier (full CF), CVSA classifier (first 3 s), Hybrid fused (full CF). Each subplot shows all / HIT / MISS breakdown plus per-class split. | Offline simulation |
| Fig 4 | **CSP filter weights** (one row per CSP type: MI / CVSA). Full-scalp topoplots via `topo_map(bg_zero=true)` — all standard 10-20 channels shown (selected = black-bold, rest = gray). Col 1: class-1 channel importance (% of total filter weight). Col 2: class-2 channel importance. Col 3: per-channel selectivity `(w1−w2)/(w1+w2)` — red = class-1-only, blue = class-2-only. Col 4: **stacked band bar** — bar height = total |CSP filter| energy in that band (% of all bands); stack split = which class dominates it (blue = class-1 components, orange = class-2 components). Skipped if no companion YAML. | CSP model (from YAML) |
| Fig 5 | **sLDA-weighted feature importance** (one row per type: MI / CVSA). Restricted to sLDA-selected features, weighted by `|sLDA coefficient|`. Col 1: channel topoplot — `Σ_k |coef_k| · |CSP_filter_k(ch)|` per channel (full 10-20 grid shown). Col 2: band bar — same sum per frequency band. Col 3: feature selection heatmap (CSP component × band) — color = `|coef|` for selected features, gray = excluded by sLDA. Skipped if no companion YAML. | CSP + sLDA models (from YAML) |
| Fig 6 | **ERD/ERS ↔ CSP weight correlation** (one panel per CSP type × frequency band). Scatter of per-channel ERD/ERS class discrimination `|ERD(c1) − ERD(c2)|` during CF vs normalised CSP filter weight `Σ|W|/max`, Pearson r in title, linear fit line. Tests the claim: channels heavily weighted by the CSP spatial filter also show stronger neurophysiological class separation. Skipped if no companion YAML. | Raw GDF signal + CSP model |
| Fig 7 | **Hemisphere ROI ERD/ERS time-courses** (one row per CSP type × non-empty ROI). CSP-selected channels grouped by 10-20 label suffix into Left (odd digit), Right (even digit), Midline (z). ERD/ERS band-averaged; baseline = cue period [−1.5s, 0s]. Left col: per-class traces ± across-channel SEM. Right col: lateralization c1(t)−c2(t). Gray solid = CF onset, blue dashed = CVSA influence window (3s), gray dotted = cue onset (−1.5s). Skipped if no companion YAML. | Raw GDF signal + CSP model |
| Fig 8 | **Within-trial classifier accuracy by fixed time bin** (one panel per paradigm group). Each trial's CF window split into fixed `BIN_WIDTH_S`-wide (default 0.5s) ABSOLUTE time bins from CF onset, up to `MAX_TIME_S` (default 5s) — NOT fractional position, so a 2s trial and a 5s trial both contribute meaningfully to the same "0-0.5s" bin; accuracy = fraction of valid frames where `argmax(P_sLDA)==target`, mean ± SEM per bin, one line per stream (MI/CVSA/Fused). The trial count contributing to each bin is annotated above every point — evaluation trials end at variable times, so later bins naturally have fewer trials still running; a hollow marker flags bins below 15% of the group's trial count as a weak/low-N average (same convention as `main_hybrid_advantage_integ`'s `temporal_significance.svg`). Answers "does the classifier get more often correct as the trial progresses?" — distinct from `main_trial_dynamics.m`, which looks at trial-ORDER effects across the session rather than within-trial timing. | Offline simulation |
| Fig 9 | **Artifact rate vs outcome** (one panel per paradigm group). Mean per-trial artifact rate during CF (fraction of frames the artifact gate froze the integrator), bar per outcome (HIT vs MISS/TIMEOUT), point-biserial correlation between artifact rate and HIT/not-HIT with a permutation p-value. Distinguishes "the integrator stalled because it was frozen" from "the classifier pointed the wrong way" as a MISS failure mode. | Offline simulation |
| Fig 10 | **Confusion matrix** (one panel per paradigm group). Rows = target class (real GDF onset code), columns = real outcome (HIT/MISS/TIMEOUT), pooled across all files of that paradigm; cell = count (row %). An aggregate accuracy (Fig 1) can hide a strong per-class asymmetry — e.g. class 1 90% HIT vs class 2 50% HIT — that matters in practice (a VR direction that almost never triggers). Skipped if no file has simulated class/outcome data. | Real GDF events + offline simulation (target class) |

Console output per file: a per-trial outcome table (HIT/MISS/TIMEOUT + time for that
event), hit rate, sample accuracy (all / HIT / MISS / per-class), within-trial accuracy
by fixed time bin (`SA by time bin (0-5.0s, step 0.5s): MI=[62%, 71%, 78%, 85%, ...]`), integrator
parameters, top-3 CSP channels per class (`[CSP-MI] top ch c769: C3>FC5>CP5 | c770: C4>FC6>CP6`),
sLDA feature count (`[sLDA-MI] 12/24 features selected by sLDA`), and
ERD trial counts (`[ERD-MI] c1=N c2=N trials`).
Set `VERBOSE_DIAGNOSTIC = true` (default `false`) for an additional per-trial diagnostic
table — very verbose, off by default.

Also per file, and pooled per paradigm as a `TOTAL` row: **ITR** (bits/trial and bits/min,
Wolpaw formula — `P=hit_tot` counting TIMEOUT as a failed trial, `N=n_cls` classes,
`T=`mean time-to-outcome over all trials) and a **chance-level significance test** (exact
one-sided binomial, decided trials only: `n=n_hit+n_miss`, `k=n_hit`, `p0=1/n_cls`, computed
via `binom_test_upper` — no Statistics Toolbox required) answering "is this accuracy
significantly above what guessing among `n_cls` classes would achieve?"

The ten figures are saved as SVG (full-screen size) into `<gdf_dir>/analysis_results/session_overview/`:
`overview_trial_accuracy.svg`, `overview_time_metrics.svg`, `overview_sample_accuracy.svg`,
`overview_csp_importance.svg`, `overview_slda_importance.svg`, `overview_erd_csp_corr.svg`,
`overview_roi_timecourse.svg`, `overview_timebin_accuracy.svg`, `overview_artifact_rate.svg`,
and `overview_confusion_matrix.svg`.
The `SHOW_FIGURES` flag at the top of the script controls whether they are also shown on screen.

**CSP weight interpretation** (Fig 4): component order `'alternate'` (MNE convention) means odd
filter rows (1, 3, 5, …) maximise class-1 variance; even rows (2, 4, 6, …) maximise class-2
variance. Channel importance = `Σ_bands Σ_{comp∈class} |W(comp, ch)|`, normalised per class.
Stacked band bar: height = total filter energy per band; color split reveals whether a frequency
range mainly helps discriminate class 1 or class 2.

**sLDA weight interpretation** (Fig 5): `w_ch(channel) = Σ_k |coef_k| · |CSP_row_k(channel)|`
sums the contribution of each sLDA-selected feature (weighted by discriminative power) to every
scalp channel. This is the direct spatial signature of the classifier's linear decision boundary —
heavier channels are literally the ones the classifier looks at most. The feature heatmap shows
exactly which (frequency band, CSP component) pairs were retained by feature selection and how
strongly they enter the decision.

After all files are processed, a **session summary** prints, per paradigm: each
file's HIT/MISS/TO counts, accuracy, and mean TTH/Tmiss/Tto, followed by a `TOTAL`
row pooling all trials of that paradigm across files (experiment-wide accuracy and
mean times).

Also saves `session_summary.mat` in `session_overview/` — required by `main_validate_counterfactual` and `main_group_analysis`. Contains per-file metrics (paradigm, n_hit, n_miss, n_to, hit_rate, tth_vals, t_miss_vals, to_vals) for all files processed in that run.

---

### `main_trial_dynamics` — within-session learning/fatigue trends

```matlab
main_trial_dynamics   % GUI file picker — any paradigm, mixed MI/CVSA/Hybrid
```

**Question this answers**: within a single evaluation session, do later trials go better or worse than earlier trials? E.g. in a 20-trial MI session, is trial 18 more likely to be a HIT than trial 2 — and is it faster? This is a trial-**order** effect (learning/fatigue across the session), not the within-trial ERD/ERS time-course (that is `topo_erders.m`).

Uses **only real GDF events** (897/898/899 relative to 781) — no CSP/sLDA/artifact simulation, no companion YAML required, so it runs even on files missing their calibration model. For each file, trial index is normalised to fractional session progress (0 = first trial, 1 = last trial), so sessions of different length pool correctly within the same paradigm.

**Why this matters beyond itself**: it is a validity/confound check for every other cross-paradigm comparison in this package. If MI/CVSA/Hybrid are always recorded in the same order within a day and performance drifts across that order, an apparent "paradigm X is worse" finding elsewhere (`main_session_overview`, `main_validate_counterfactual`, `main_group_analysis`) could really be a time-of-day/fatigue effect confounded with recording order rather than a true paradigm effect — Fig 4 below checks exactly this.

**Four figures** (saved under `<gdf_dir>/analysis_results/trial_dynamics/`):

| Figure | Content |
|---|---|
| `trial_dynamics_trend.svg` | **Main figure.** One combined plot, one line per paradigm: accuracy binned into 5 quintiles of session progress, mean ± binomial SE per bin. Significance comes from a point-biserial correlation between raw trial position and HIT/not-HIT (every trial, not the binned means) with a permutation p-value — shown in the legend per paradigm. Positive r = improves over the session (learning), negative r = degrades (fatigue). |
| `trial_dynamics_half_split.svg` | One panel per paradigm found: paired per-file dots (1st-half accuracy → 2nd-half accuracy) + paradigm mean line; sign-flip permutation test on the paired delta (2nd half − 1st half). A simple file-level cross-check of the Fig 1 trend. |
| `trial_dynamics_tth_trend.svg` | One panel per paradigm: time-to-HIT vs normalised trial position (all HIT trials pooled across that paradigm's files), linear fit, Pearson r with a permutation p-value (shuffling trial position). Negative r = getting faster over the session = learning. Same question as Fig 1, but on speed instead of accuracy. |
| `trial_dynamics_recording_order.svg` | Accuracy vs **chronological recording order across all selected files, any paradigm** (timestamp parsed from the `bag_bci` `<subject>.<YYYYMMDD>.<HHMMSS>.<...>` filename convention), one marker shape per paradigm, linear fit, Pearson r with permutation p-value. Skipped if a timestamp can't be parsed from every filename or fewer than 3 files are selected. |

**Console output**: per-file 1st-half/2nd-half accuracy table; per-paradigm accuracy-trend statistics (r, p-value); per-paradigm half-split statistics (mean delta, p-value); per-paradigm TTH-trend statistics (r, p-value); chronological order table + recording-order trend (r, p-value) + **mean recording-order index per paradigm** — if these means are far apart across paradigms, flag possible order confounding before trusting any cross-paradigm comparison done elsewhere.

Also saves `trial_dynamics_summary.mat` — a flat per-trial table (file_label, paradigm, trial_idx, n_trials, frac, outcome_code, dt) plus `file_info` (file-level chronological table: file_label, paradigm, timestamp, acc) — under `trial_dynamics/`.

**Interpretation caveat**: with few trials per file the half-split is noisy; pool more files of the same paradigm (and ideally more subjects via a multi-subject extension) before drawing conclusions from a single session.

---

### `main_hybrid_advantage_probs` — classifier-probability analysis of hybrid recordings

```matlab
main_hybrid_advantage_probs    % GUI file picker — hybrid GDF files only
```

> **Requires hybrid GDF files.** Designed to answer: **is the fused classifier driving the integrator better than
MI or CVSA alone — faster to the correct direction, fewer frames pointing the
wrong way?**

Accepts only **hybrid** GDF files. Runs the pipeline once and computes, for every
valid CF frame of each trial, the following per-stream metrics over the CF window:

| Metric | Definition |
|---|---|
| **mean P(target)** | mean over valid CF frames of P(target class) |
| **frame accuracy** (acc) | fraction of valid CF frames where P(target) > 0.5 |
| **acc_buf** | fraction of valid CF frames where buffer(target) > p_rest (= 0.5) |
| **TTH** | first frame (in seconds) where P(target) > 0.5 — lower = reacts faster |
| **false-direction rate** | fraction of valid CF frames where P(non-target) > 0.5 — lower = fewer misleading frames |

Streams: **MI**, **CVSA** (full CF), **CVSA-inf** (first `cvsa_influence` seconds only),
**fused** (Bayesian LOP), **buffer** (leaky WTA integrator output).

**Six figures:**

| Figure | Content |
|---|---|
| Fig 1 | **Mean P(target) scatter** — one point per trial per stream (5 panels: MI / CVSA / CVSA-inf / fused / buffer). Color = outcome (green HIT / orange MISS / red TO); marker = cued class. Dashed lines = per-outcome group means. |
| Fig 2 | **Frame accuracy scatter** — same layout as Fig 1, y-axis = fraction frames where P > 0.5. |
| Fig 3 | **Does CVSA fusion help? (binary view)** — 2×2 layout: (1,1) **frame-level rescue/hurt bar chart** — counts, over all trials within the CVSA-influence window, of CF frames where `P_MI(target)` and `P_fused(target)` agree/disagree on being `>0.5`: *rescued* (MI wrong → fused correct), *hurt* (MI correct → fused wrong), *both correct*, *both wrong*. Title shows the net effect (`rescued − hurt`) and a verdict ("CVSA helped"/"CVSA hurt"/"neutral"); (1,2) per-trial **CVSA-fusion advantage** = mean(P_fused(target) − P_MI(target)) over the CVSA-influence window — positive = fusion pushed toward the target class more than MI alone, negative = fusion pulled away; (2,1)/(2,2) mean ± SEM of P_MI(target) and P_fused(target) over CF time, one panel per cued class (each truncated to that class's longest trial, so the SEM band renders for both classes), with a vertical line marking where the CVSA prior fades to zero (`cvsa_influence`). |
| Fig 4 | **CF trial duration & CVSA effect on HIT trials, by hit timing** — two panels: (1) per-trial CF duration scatter colored by outcome, with per-outcome dashed mean lines and a dotted line at `cvsa_influence` (split out from Fig 3 into its own figure); (2) **HIT trials only** (MISS/TIMEOUT excluded) — mean `P_MI(target) → P_fused(target)` slopegraph (± SEM, delta annotated), pooling CF frames within the first `cvsa_influence` seconds, split into trials that reached HIT *within* `cvsa_influence` (`duration < cvsa_influence`) vs *after* it (`duration ≥ cvsa_influence`) — tests whether fusion helped fast hits cross threshold sooner (positive delta) vs corrected/hindered MI early in slower hits. |
| Fig 5 | **Does CVSA help, broken down by MI/CVSA agreement** — every CF frame within the CVSA-influence window is classified into one of 4 categories by whether `P_MI(target)>0.5` and `P_CVSA(target)>0.5`: *agree-correct* (both right), *MI correct, CVSA wrong*, *MI wrong, CVSA correct* (the "rescue" case), *agree-wrong* (both wrong). Three panels: (1) scatter of `P_MI(target)` vs `P_fused(target)` per frame, colored by category, with the `y=x` diagonal and crosshairs at 0.5; (2) slopegraph of mean `P_MI(target) → P_fused(target)` per category (± SEM), with the mean delta annotated above each pair — this is the "small improvement" view; (3) for the two disagreement categories only, `P_fused(target) − P_MI(target)` vs the cosine-annealed CVSA weight `α(t)`, showing how the rescue/cost effect scales with `α`. |
| Fig 6 | **Cluster-based permutation test** (Maris & Oostenveld style, `utils/cluster_permutation_test.m`, no Statistics Toolbox) on `P_fused(target)-P_MI(target)` pooled across trials of both classes within the CVSA-influence window — tests whether the CVSA-fusion effect is a genuine, temporally-localised period of influence rather than noise, correcting for the multiple-comparisons problem of testing every frame independently. Per-frame statistic = mean/SEM across trials; candidate clusters = contiguous frames where `\|t\|≥2.0` (a fixed conventional cluster-forming threshold, **not** calculated from this session's actual degrees of freedom — it only decides candidate cluster boundaries, not the final significance); cluster p-value = permutation on the max `\|cluster mass\|` under whole-trial sign-flipping (preserves temporal autocorrelation within a trial). Two panels: pooled mean±SEM P_MI/P_fused with significant cluster(s) shaded, and the t-statistic curve itself with the cluster-forming threshold marked. Each cluster also reports `mean_delta`/`cohen_d`/`ci_lo`/`ci_hi` (trial-level mean of `P_fused-P_MI`, Cohen's d, and a percentile bootstrap 95% CI — all collapsed over just that cluster's own time range, not the whole CVSA-influence window) printed alongside its p-value in the console cluster list. |

**Console output** (per trial, then summaries):
- Per-trial table: `#, cue, result, mMI, aMI, mCV, aCV, mCVi, aCVi, mFus, aFus, mBuf, aBuf`
- Per-class summary: mean of all metrics split by cued class
- Per-outcome summary: mean of all metrics split by HIT/MISS/TIMEOUT
- CVSA-fusion advantage summary: count of trials where fusion helped (>0) vs hurt (<0), the mean delta, and a **within-session sign-flip permutation test + Cohen's d** on the trial-level advantage (trial is the unit, 2000 sign-flips, no Statistics Toolbox) — this per-subject test is well-powered (uses every trial in the session) and is what `main_group_analysis`'s Stouffer meta-analysis combines across subjects, instead of relying on that group test's own coarse 2^n_subj sign-flip floor
- CVSA-fusion frame-level effect (binary): rescued/hurt/both-correct/both-wrong frame counts within the CVSA-influence window, plus the net effect and verdict
- Cluster-based permutation test: list of candidate clusters (time range, mass, p-value), flagged significant at p<0.05
- HIT trials by hit timing: trial/frame counts, mean `P_MI -> P_fused`, and delta for HIT-within-`cvsa_influence` vs HIT-after-`cvsa_influence` groups
- CVSA vs MI agreement breakdown: per-category frame counts, mean `P_MI -> P_fused`, and delta; plus rescue effect (MI wrong, CVSA correct) and cost effect (MI correct, CVSA wrong)

Saves `hybrid_advantage_probs_summary.mat` (per-file: mean P/acc per stream, fusion advantage, rescue/cost frame counts and deltas, `fus_adv_pval`/`fus_adv_cohend` (within-session sign-flip test + Cohen's d on the trial-level fusion advantage), plus `cluster_time_axis`/`cluster_diff_mean`/`cluster_n_trials`/`cluster_cvsa_inf` — this session's own mean `P_fused-P_MI` curve, used by `main_group_analysis` for a cross-subject cluster test).

**Interpreting mean P:**

| Value | Interpretation |
|---|---|
| > 0.5 + HIT | Classifier strong, integrator fast |
| > 0.5 + MISS | Classifier correct direction but integrator too slow → tune `k_gain`/`buffer_size` |
| ≈ 0.5 + MISS | Near-chance; integration goes nowhere |
| < 0.5 + MISS | Classifier points wrong direction → genuine classifier failure |

The five figures are saved as SVG (full-screen size) into `<gdf_dir>/analysis_results/hybrid_advantage_probs/`;
the `SHOW_FIGURES` flag at the top of the script controls whether they are also shown on screen.
Also saves `hybrid_advantage_probs_summary.mat` in that folder (per-file: mean P/accuracy per stream,
CVSA-fusion advantage, frame-level rescue/hurt counts, rescue/cost deltas) for cross-subject use by
`main_group_analysis`.

---

### `main_hybrid_advantage_integ` — counterfactual: hybrid vs MI-only vs CVSA-only

```matlab
main_hybrid_advantage_integ    % GUI file picker — hybrid GDF files only
```

> **Requires hybrid GDF files.** Designed to answer: **given the exact same classifier outputs recorded during a
hybrid session, would trials have ended differently (HIT/MISS/TIMEOUT, and at
what time) if the integrator buffer had been driven by MI alone or CVSA alone,
instead of the Bayesian-fused signal?**

Accepts only **hybrid** GDF files. For each file, the pipeline is run once to
get `p_mi_aligned`, `p_cvsa_aligned`, `art_flags`, and `int_cfg` (identical to
`main_simulate`/`main_hybrid_advantage`), then `integrate_signal` is called
**three times** with the same `int_cfg` (same `buffer_size`, `k_gain`,
`thresholds`, `classes`, `init_val`) — only the driving stream differs:

| Run | Stream fed to the leaky-WTA buffer |
|---|---|
| Hybrid | Bayesian LOP fusion of MI+CVSA (as actually run online) |
| MI-only | raw MI sLDA output |
| CVSA-only | raw CVSA sLDA output |

This gives three directly-comparable `trials` arrays sharing the same
`n_pre`/`n_cf`/`artifact`/`target_class` per trial — only `.raw`/`.integrated`/
`.normalized*` differ. The **simulated outcome** of each counterfactual is
derived from its own integrator buffer (`integrated(:,i) >= thresholds(i) - 5e-3`,
first class to cross wins; neither crossing within the CF window = TIMEOUT) —
this is independent of the *real* GDF outcome, which only reflects the hybrid
run.

**Per-file figures** (saved under `per_file/`):

| Figure | Content |
|---|---|
| `signals_<basename>.svg` | One panel per 781 trial — overlaying the leaky-integrator P(c1) control signal for all three streams (Hybrid/MI-only/CVSA-only) plus the thresholds and `p_rest` line. Background colour = REAL hybrid outcome (HIT=light green/MISS=light red/TIMEOUT=light yellow); the title shows the *simulated* outcome and time-to-event for each stream. |
| `advantage_<basename>.svg` | 2×3 per-trial breakdown: (1,1) mean sLDA P(target) per trial (raw classifier signal); (1,2) **mean integrator buffer(target class) per trial** — the actual control-signal level for the cued class relative to the win threshold; (1,3) simulated time-to-outcome per trial; (2,1) simulated-outcome heatmap (trial × stream); (2,2) outcome counts per stream; (2,3) **per-trial control-signal advantage** (Hybrid − MI-only) and (Hybrid − CVSA-only) — positive = fusion pushed the buffer higher for the cued class that trial. |

**Aggregate figures** (saved at the top level):

| Figure | Layout | Content |
|---|---|---|
| `summary_hybrid_vs_unimodal.svg` | 3×3 | Outcome counts; accuracy (title states the Hybrid−MI/CVSA delta); mean TTH (title states speed advantage); mean time-to-MISS; **mean time-to-TIMEOUT** (restored); two "who hits?" concordance bars with rescued/cost/net annotated; mean buffer(target class) ± SEM with sign-flip permutation significance brackets and p-values in title; per-trial buffer-advantage histogram with count-above-zero annotation. |
| `temporal_significance.svg` | 3×1 | (1) Mean ± SEM buffer(target class) aligned to CF onset for all three streams; (2) **Hybrid−MI delta trajectory** — mean ± SEM band, zero line = no advantage, coloured shading = p<0.05, peak annotated with magnitude and time; (3) same for **Hybrid−CVSA**. Directly answers: IS there an advantage, HOW MUCH, and WHEN. Grey patches = fewer than 15% of trials still running. |
| `advantage_deep.svg` | 2×2 | (1,1) **Cumulative HIT fraction over time** with final accuracy-gap text box; thin dashed lines per class; (1,2) per-trial scatter MI-only vs Hybrid mean buffer(target) — count above diagonal annotated in title; (2,1) same scatter for CVSA-only vs Hybrid; (2,2) accuracy bars with **bootstrap 95% CI** (2000 resamples, asymmetric). |
| `confusion_matrix.svg` | 1×3 | **Per-stream confusion matrix** — one panel per stream (Hybrid/MI-only/CVSA-only): rows = target class, columns = simulated outcome (HIT/MISS/TIMEOUT), cell = count (row %). Since the SAME buffer/thresholds drive all three streams, any per-class bias difference across panels reflects the classifier signal itself (e.g. CVSA-only skewed toward one attention direction while MI-only is balanced), not the integration settings — a class-bias check for the "fixed MI-tuned parameters make the three paradigms comparable" assumption. |

**Console output**: per-trial table; aggregate counts/accuracy/mean times; "who hits?" breakdown; per-class accuracy; **sign-flip permutation test** on mean buffer(target) advantage (Hybrid vs MI-only, Hybrid vs CVSA-only); **Cohen's d** (paired effect size) on buffer advantage; **bootstrap 95% CI** on per-stream accuracy and on accuracy deltas; **ITR** (bits/trial and bits/min, Wolpaw formula — `P=acc(s)` counting TIMEOUT as failure, `N=n_cls`, `T=`mean time-to-outcome for that stream) and a **chance-level significance test** (exact one-sided binomial, decided trials only) per stream — answers whether the MI-tuned integration settings put each stream above chance, and how their communication rate compares once trial speed is folded in.

**Statistics (no Stats Toolbox required)**: the sign-flip permutation test (`sign_flip_test` local helper, 2000 permutations for global tests, 500 for per-frame temporal test), the bootstrap CI (percentile method, 2000 resamples), and the chance-level test (`binom_test_upper`, exact binomial via `gammaln`) are all implemented without MATLAB's Statistics and Machine Learning Toolbox. The temporal significance figure is explicitly labelled as pointwise/uncorrected — for small n (40 trials) it is exploratory; sustained runs of low p-values are more meaningful than isolated frames.

All figures are saved as SVG (full-screen size) under
`<gdf_dir>/analysis_results/hybrid_advantage_integ/`; the `SHOW_FIGURES` flag at the
top of the script controls whether they are also shown on screen.

Also saves `counterfactual_summary.mat` in `hybrid_advantage_integ/` — required by `main_validate_counterfactual` and `main_group_analysis`. Contains aggregate accuracy, TTH, bootstrap 95% CI, and outcome counts for the three streams (index 1=Hybrid, 2=MI-only, 3=CVSA-only), plus `itr_bits_trial`, `itr_bits_min`, and `chance_p` (each `[1x3]`, same stream order).

---

### `main_validate_counterfactual` — real session vs offline simulation

```matlab
main_validate_counterfactual   % GUI picks session_summary.mat then counterfactual_summary.mat
```

Loads the `.mat` outputs of `main_session_overview` (all three paradigm sessions) and `main_hybrid_advantage_integ` (hybrid counterfactual), then compares per paradigm:

| Comparison | What it tests |
|---|---|
| Real MI session acc vs simulated MI-only (from hybrid data) | Is the counterfactual a valid proxy for a dedicated MI session? |
| Real CVSA session acc vs simulated CVSA-only (from hybrid data) | Same for CVSA |
| Real Hybrid session acc vs simulated Hybrid | Sanity check — should agree closely |

**One figure** (`counterfactual_validation.svg`, saved under a sibling `validate_counterfactual/` folder next to `hybrid_advantage_integ/`):
- Grouped bar chart: real (solid) vs simulated (lighter + dashed) accuracy per paradigm, with bootstrap CI and Δ annotation
- Scatter: simulated vs real accuracy for the three paradigms with identity line and ±10% tolerance band
- TTH comparison: same grouped bar layout for mean time-to-HIT

**Interpretation**: points near the identity line → simulation valid → the counterfactual claim ("hybrid beats MI-only, given identical classifier streams") is representative of true between-session performance. Large divergence → dual-task cognitive cost degrades the MI signal in hybrid sessions.

**Console output**: per-paradigm table of real acc, simulated acc, delta, and n_trials; TTH comparison; interpretation guidance.

---

### `main_threshold_sweep` — HIT/TIMEOUT rate over a (threshold_class1, threshold_class2) grid

```matlab
main_threshold_sweep   % GUI file picker — any paradigm, mixed MI/CVSA/Hybrid
```

**Question this answers**: how sensitive is the experiment's outcome to the choice of the two per-class integrator thresholds (`int_cfg.thresholds`)? Re-evaluates every trial's **already-simulated** integrator buffer (`integrate_signal.m` is run exactly once per file — its dynamics don't depend on `thresholds` at all, only the pass/fail check does) against every `(threshold_class1, threshold_class2)` combination on a fixed grid, instead of a classic ROC/AUC: with two independent per-class thresholds and a TIMEOUT outcome, there's no single binary positive/negative label, so a 2D heatmap over the threshold plane is the natural generalisation.

**Outcome rule** per `(th1, th2)`, mirrors `Training.cpp::is_target_hit` / `main_hybrid_advantage_integ.m`'s `classify_trial_outcome`: within the CF window, the first class `i` whose `integrated(:,i) >= th(i) - 5e-3` is reached wins — HIT if that class is the cue's target, MISS if it's the other class, TIMEOUT if neither threshold is ever reached. Implemented with `cummax` per trial so the whole grid sweep needs no pipeline re-run, just index lookups against each trial's already-computed buffer trace.

**Grid**: fixed `0.55:0.05:1.00` on both axes regardless of the recording's real threshold, so every saved image has the same scale and is directly comparable across files/paradigms. The real recording-time threshold is marked with a red square on every image (for pooled groups spanning files with slightly different real thresholds, the marker is the mean).

**Paradigm handling**: each paradigm actually present in the selection gets its **own** image — MI/CVSA/Hybrid are never silently merged, since they can behave very differently (e.g. Hybrid's CVSA-fusion window). A pooled `ALL` image is also produced. Additionally, **every individual GDF file** gets its own image in a `per_file/` subfolder.

**Output**, under `<gdf_dir>/analysis_results/threshold_sweep/`:

| File | Content |
|---|---|
| `threshold_sweep_<paradigm>.svg` | One per paradigm present — HIT rate heatmap + TIMEOUT rate heatmap, side by side |
| `threshold_sweep_ALL.svg` | Same, pooled across all selected files/paradigms |
| `per_file/threshold_sweep_<basename>.svg` | Same, for one individual GDF file |
| `threshold_sweep_summary.mat` | All grids (`sweep_by_group`, `sweep_by_file`), the shared `th_vals`, and `real_thresholds_per_file` |

**Console output**: per paradigm/file, the best HIT-rate cell on the grid and the HIT/TIMEOUT rate at the (nearest grid point to the) real recording-time threshold.

---

### `main_roc_analysis` — classifier-level ROC/AUC (pre-integrator)

```matlab
main_roc_analysis   % GUI file picker — any paradigm, mixed MI/CVSA/Hybrid
```

**Question this answers**: how separable are the two classes according to the classifier alone, independently of the integrator/buffer/threshold machinery? Complements `main_threshold_sweep`, which sweeps the *integrator's* buffer thresholds and is deliberately not a classic ROC (two independent per-class thresholds + a TIMEOUT outcome don't reduce to one binary label). This script instead pools the raw per-frame sLDA probability `P(class 1)` across CF frames and sweeps every possible scalar threshold in `[0,1]` — a genuine ROC/AUC, computed with a hand-rolled non-parametric routine (`utils/compute_roc_curve.m`, no Statistics Toolbox required).

**Signal used per paradigm**, reusing `trials(t).raw` from `integrate_signal.m` (already paradigm-correct):
- MI paradigm → raw MI sLDA `P(c)`
- CVSA paradigm → raw CVSA sLDA `P(c)`
- **Hybrid paradigm → the cosine-annealed Bayesian-fused `P(c)` (`bayesian_fuse.m`), treated as the output of a single hypothetical classifier** — the only sense in which "ROC of the hybrid" is defined, since MI-only/CVSA-only already have their own paradigm's files.

**Frame pooling**: excludes the `N_PRE` reset frame and any artifact-flagged frame (corrupted EEG upstream of the classifier, not a genuine class-conditional sample — exactly the frames the real system already ignores via the artifact gate). Alongside the ALL-frames pool, a **windowed** pool is also built per paradigm group: only frames with time-since-CF-onset `< cvsa_influence` (frame `k` of a trial's CF, `k=1..n_cf`, is at `t=(k-1)/framerate` — `k=1` is the CF onset chunk, matching `bayesian_fuse.m`'s own `t=0` convention). This is the ONLY window where a fused/Hybrid classifier can actually differ from pure MI (`bayesian_fuse.m` forces `alpha=0` — `P_fused==P_MI` exactly — from `cvsa_influence` onward), so the ALL-frames AUC dilutes any early classifier-level advantage with the later, by-construction-identical portion of every trial; `main_group_analysis`'s Fig 13c is built from this windowed pool. Printed to console per group as an extra line (`windowed (first Xs of CF): n_frames=... AUC=...`) alongside the ALL-frames AUC.

**Paradigm handling**: each paradigm present gets its own ROC image, plus every individual GDF file gets its own image in `per_file/`. Unlike `threshold_sweep`'s pooled `ALL` image, pooling raw scores across MI-only/CVSA-only/Hybrid-fused would mix different classifiers on different probability scales into a meaningless single AUC — instead `roc_ALL_paradigms.svg` **overlays each paradigm's own ROC curve** (own AUC each) on one plot, the valid "general" cross-paradigm comparison. Every curve also marks the natural `p=0.5` operating point (matches `main_session_overview`'s frame-accuracy convention) with sensitivity/specificity annotated.

**Output**, under `<gdf_dir>/analysis_results/roc_analysis/`:

| File | Content |
|---|---|
| `roc_<paradigm>.svg` | One per paradigm present |
| `roc_ALL_paradigms.svg` | Overlay of each paradigm's own ROC curve (NOT a pooled curve) |
| `per_file/roc_<basename>.svg` | Same, for one individual GDF file |
| `roc_summary.mat` | `fpr`/`tpr`/`thr`/`auc` per group and per file, **plus the raw pooled `scores`/`labels` per paradigm group** — kept so `main_group_analysis` can validly re-pool a subject's sessions (same deployed classifier/calibration) into a per-subject ROC; **plus** `window_s`/`auc_win`/`scores_win`/`labels_win` per group — the same frames restricted to the first `cvsa_influence` seconds of each trial's CF, used for the windowed classifier comparison (Fig 13c) |

**Console output**: per paradigm/file, AUC plus sensitivity/specificity at the `p=0.5` operating point.

---

### `main_group_analysis` — multi-subject aggregation

```matlab
main_group_analysis()                                          % default recordings root (hardcoded), ALL subjects
main_group_analysis('/path/to/recordings')                     % given root, ALL subjects
main_group_analysis('/path/to/recordings', {'a1','a2','a3'})   % only these subjects
main_group_analysis('/path/to/recordings', {'a1','a2'}, true)  % also show figures on screen
main_group_analysis([], {'a1','a2'})                            % default root, subset of subjects
main_group_analysis('/path/to/recordings', {'a1','a2'}, false, 'well')  % tag every output file 'well'
```

No GUI folder picker: `root_dir` defaults to a hardcoded `DEFAULT_ROOT` at the top of the script (edit it directly, or pass a path) rather than prompting. `subjects_filter` is an explicit opt-in filter (cellstr of subject IDs); omit it (or pass `{}`) to use every subject found under the root.

Run this **after** `run_subject_analysis` has finished for all subjects.
It reads the `.mat` files already saved by the per-subject scripts — it does not re-run any pipeline.

**`group_tag`** (optional 4th argument, default `'all'`): stitched into every output filename (`'<num>_<tag>_<name>.svg'`, `'group_summary_<tag>.mat'`) and into the in-figure caption cross-references between sibling figures (e.g. "full test in `06_well_speed_itr.svg`"). Each call still analyzes exactly ONE group of subjects — running several groups (e.g. a manually-assigned "well"/"bad" performer split, on top of the default "all") is orchestrated by the *caller*, calling `main_group_analysis` once per group with a different `subjects_filter`/`group_tag`; see `batch_group_analysis.m` below, the intended entry point for this. The function now also **returns** `group_summary` so the caller can collect it (e.g. to feed a "well" and a "bad" run into `well_vs_bad_comparison.m`, described below).

**Batch scripts** (`batch_run_subjects.m` and `batch_group_analysis.m`, both in this folder): edit `RECORDINGS_ROOT` and `SUBJECTS = {'a1','a2',...}` at the top of each, then run directly — no GUI, no arguments to pass. `batch_run_subjects` calls `run_subject_analysis(<day>/evaluation)` once per (subject, day) found under each listed subject (always a direct evaluation folder, so it never triggers `run_subject_analysis`'s own multi-day `listdlg`). `batch_group_analysis` calls `main_group_analysis(RECORDINGS_ROOT, SUBJECTS, ..., 'all')` and `group_topo_erders`, both restricted to the same `SUBJECTS` list. It also defines `WELL_SUBJECTS`/`BAD_SUBJECTS` (cellstr, e.g. `{'a1','a3'}`; leave `{}` to skip) — manually-assigned "strong"/"weak" performer subsets — and, for each non-empty one, calls `main_group_analysis(RECORDINGS_ROOT, WELL_SUBJECTS, ..., 'well')` / `(..., BAD_SUBJECTS, ..., 'bad')`, skipping gracefully (no error) if a list doesn't match any subject found on disk. If both `well` and `bad` ran, it then calls `well_vs_bad_comparison` on the two returned summaries. Keep `SUBJECTS` in sync between `batch_run_subjects.m` and `batch_group_analysis.m` to analyse and aggregate the same cohort.

Recursively scans the selected root for every `session_overview/session_summary.mat`, `hybrid_advantage_integ/counterfactual_summary.mat`, `hybrid_advantage_probs/hybrid_advantage_probs_summary.mat`, `topo_erders/topo_erders_summary.mat`, and `roc_analysis/roc_summary.mat` found anywhere underneath it (expects the layout `<root>/<subject>/.../analysis_results/...`). Subject ID = first path component directly under the root. If a subject has multiple matching `.mat` files (e.g. several evaluation sessions), real-session counts are summed, counterfactual accuracy is pooled weighted by `n_trials`, fusion-mechanism/ERD metrics are averaged (frame counts summed), and ROC is built by **re-pooling that subject's raw per-frame `(score, label)` pairs across all their sessions** (valid, since the same deployed classifier/calibration produced every one of a subject's frames) before computing per-subject statistics.

Every real-session metric (accuracy, timeout rate, TTH/T-miss/T-timeout) is shown **two ways** in the same figure: one grouped-bar per subject (subject is the statistical unit for the sign-flip tests elsewhere in this script), plus a final **"GRAND"** group built by **pooling every subject's trials together** — the population-level estimate, which can differ from the mean of the per-subject bars when subjects contribute unequal trial counts.

**Two accuracy conventions**, computed and labelled side by side everywhere (per-subject arrays, console tables, group-level tests, figures, `group_summary.mat`) — never silently picking one:
- `*_acc` = `n_hit / n_total` (**TIMEOUT counted as a failure**) — the original convention
- `*_acc_dec` = `n_hit / (n_hit + n_miss)` (**decided trials only**, TIMEOUT excluded) — "how good is the classifier/integrator when it actually reaches a decision", independent of how often it fails to decide at all
- `*_fp_rate` = `n_miss / n_total` (**false positive rate**) — the WRONG class's threshold was reached with confidence, distinct from TIMEOUT (neither threshold reached)
- `*_to_rate` = `n_to / n_total` (already existed)

Both real-session and counterfactual (simulated) data get this same 4-way split, each with a trial-pooled `GRAND` value.

Best/worst subject ranking is done **per paradigm** (real: MI/CVSA/Hybrid) and **per stream** (simulated: Hybrid/MI-only/CVSA-only) — not a single global pick — ranked on decided-trial accuracy, with **ties kept together** (more than one subject can be marked best/worst). Printed as console ranking tables and marked with a gold star (best) / red triangle (worst) on the corresponding figures.

**Eighteen figures** (numbered 01-12 plus 04b/13a/13b/14a/14b/15, saved under `<root>/group_analysis/`, each filename prefixed `NN_` with real-data figures first, then counterfactual, then the purely statistical ones — same order as this list; the redundant `group_` word was dropped from filenames since the numeric prefix already sorts them):
- `01_real_session_overview.svg` — two panels: per-subject real-session **accuracy** (TIMEOUT = fail) and **timeout rate** per paradigm (MI/CVSA/Hybrid), grouped bars + trial-pooled `GRAND` group (visually separated by a vertical line and a thicker bar edge)
- `02_real_times.svg` — one panel per paradigm: per-subject **TTH / T-miss / T-timeout** (mean seconds, pooling that subject's trials of that outcome), grouped bars + trial-pooled `GRAND` group
- `03_real_accuracy_breakdown.svg` — **decided-trial accuracy** (TIMEOUT excluded) on top (full width), **false-positive rate** and **timeout rate** below (side by side), per paradigm, per subject + `GRAND`; gold star / red triangle mark the best/worst subject *for that paradigm* (ties shown together)
- `04_counterfactual_advantage.svg` — left panel: per-subject counterfactual (**SIMULATED**, not real) accuracy (Hybrid/MI-only/CVSA-only), TIMEOUT = fail; right panel: per-subject Hybrid−MI and Hybrid−CVSA accuracy deltas, plus group mean ± SEM and a **sign-flip permutation test across subjects** (not trials) — once N>1 subject is available this is the statistically meaningful unit, unlike the per-trial test already reported inside `main_hybrid_advantage_integ` — plus its **paired Cohen's d** (subject is the unit): how large the advantage is, not just whether it's non-zero, shown in the panel title alongside the p-value stars. SIMULATED means: the SAME integrator/buffer/thresholds are driven counterfactually by the Bayesian-fused, raw MI-only, or raw CVSA-only signal on the SAME Hybrid-session trials — isolating the fusion algorithm's own effect from any real-session differences in trial count, timing, or thresholds.
- `04b_time_in_correct_zone.svg` — complements `04`'s binary hit/miss/timeout view: for EVERY trial (hit/miss/timeout alike, whole actual CF duration, no outcome-based truncation), what fraction of that trial's own duration did the target-class control signal spend on the correct side of 0.5? RAW (pre-integrator sLDA/fused `P(target)`) and INT (post-integrator `buffer(target)`, the actual VR control signal) versions, each with per-subject bars + a paired Hybrid-MI-only/Hybrid-CVSA-only delta scatter + group sign-flip test + Cohen's d. Computed by `main_hybrid_advantage_integ.m`, saved as `frac_correct_raw`/`frac_correct_int` in `counterfactual_summary.mat`; skipped if not yet regenerated with these fields.
- `05_simulated_accuracy_breakdown.svg` — same layout as `03_real_accuracy_breakdown.svg`, but for the 3 counterfactual streams instead of the 3 real paradigms; best/worst marked *per stream*; same accuracy-vs-timeout-rate scatter panel as `03`.
- `06_speed_itr.svg` — 2×2: is Hybrid not just more accurate but also faster and a better bits/min channel? Left column = real sessions, right column = counterfactual; top row = per-subject Δ TTH scatter (Hybrid−MI, Hybrid−CVSA, negative = faster) + group mean±SEM with sign-flip stars; bottom row = the same for ITR (Wolpaw bits/min, `N_CLS=2` assumed). `T` (the ITR denominator) is pooled from TTH/T-miss/T-timeout weighted by each category's own trial count. None of the other figures test TTH/ITR deltas for significance, so this is the only place these two questions get an explicit answer.
- `07_early_late_keepup.svg` — **counterfactual only** (real MI-only/CVSA-only sessions don't share trial indices with Hybrid, so this cross-stream comparison on the SAME trial isn't possible there): splits Hybrid's own HIT trials into "early" (resolved within `cvsa_influence`) vs "late" (resolved after, alpha~0 already) and checks whether MI-only/CVSA-only ALSO hit on those exact same trials — "would the unimodal stream have kept up?" Two slopegraph panels + group-level sign-flip test on the paired delta. A low MI keep-up rate specifically in the "early" group (not "late", where Hybrid~MI already) would mean CVSA specifically enabled the fast resolutions. Skipped if no subject's `counterfactual_summary.mat` has the raw per-trial `oc_hyb`/`oc_mi`/`oc_cvsa` arrays yet.
- `08_fusion_mechanism.svg` — left panel: per-subject CVSA-fusion advantage `mean(P_fused-P_MI)` with group mean ± SEM, sign-flip test, and Cohen's d; right panel: per-subject rescue-delta vs cost-delta slopegraph with a group-level sign-flip test + Cohen's d on `(rescue-cost)` across subjects — tests whether the frame-level fusion mechanism described per-file by `main_hybrid_advantage_probs` generalises across the cohort
- `09_rescue_hurt_breakdown.svg` — "on average, in how many SAMPLES does the CVSA help vs hurt vs not matter?" Three stacked panels reusing `main_hybrid_advantage_probs`'s rescued/hurt/both-correct/both-wrong frame classification (MI-alone-correct vs fused-correct), one per frame window: within `cvsa_influence` (all trials), whole trial (no time restriction), and from `cvsa_influence` onward (title states how many of the total trials actually contribute a post-window frame). Per-subject % of classified frames + trial-pooled GRAND. Skipped if no subject's `hybrid_advantage_probs_summary.mat` has the `_all`/`_post` window variants yet.
- `10_fusion_cluster_test.svg` — **cross-subject cluster-based permutation test** (same `utils/cluster_permutation_test.m` as `main_hybrid_advantage_probs`'s Fig 6, one level up): each subject contributes their own mean `P_fused(target)-P_MI(target)` curve (interpolated onto a common time grid, subject = permutation unit) — tests whether the temporally-localised CVSA-fusion effect found per-session generalises across the cohort, not just noise in one session. Two panels: cross-subject mean ± SEM curve with significant cluster(s) shaded, and the t-statistic curve. Any significant cluster also reports `mean_delta`/`cohen_d`/`ci_lo`/`ci_hi` localised to its own time range (subject is the unit here) in the console cluster list. Skipped if no subject's `hybrid_advantage_probs_summary.mat` has been regenerated with cluster-test data yet (re-run `main_hybrid_advantage_probs` after upgrading).
- **Console-only companion, no new figure — group-level CVSA-fusion meta-analysis (Stouffer)**: combines each subject's own **within-session sign-flip test p-value** (`fus_adv_pval`, computed by `main_hybrid_advantage_probs`, trial is the unit) into one group-level one-sided p-value testing "does CVSA-fusion help" (H1: mean advantage > 0). Exists because the group sign-flip test in `08_fusion_mechanism.svg` above has a permutation null of only 2^n_subj sign patterns — with a handful of subjects it structurally cannot reach p<0.05 no matter how strong the effect (a floor on the test's own resolution, not evidence the effect is weak). The meta-analysis sidesteps this by using each subject's own well-powered (many-trial) within-session evidence instead of first collapsing every subject to a single mean. Method: each session's two-sided `fus_adv_pval` is converted to one-sided (in the "CVSA helps" direction, using that session's own `fus_adv_mean` sign — a session pointing the wrong way gets a weak/large one-sided p rather than being dropped), a subject's own multiple sessions (if any) are combined weighted by `sqrt(n_trials)`, then subjects are combined with **equal weight** (subject is the unit, matching every other group-level test here) via Stouffer's Z method (hand-rolled `norminv`/`normcdf` via `erfcinv`/`erfc`, both base MATLAB — no Statistics Toolbox). Valid at any cohort size, including large ones — it stays a genuine complement to the group sign-flip test (not just a small-n workaround) since it uses the full continuous per-subject evidence rather than one collapsed mean per subject. Printed per-subject (within-subject combined p, total trial count) plus the final combined p; saved in `group_summary.mat` as `fus_adv_meta_p_subj`/`fus_adv_meta_n_subj` (per subject) and `p_fusadv_meta_group` (combined). Requires `hybrid_advantage_probs_summary.mat` regenerated with `fus_adv_pval` (older summaries are skipped, not an error, same backward-compatibility convention as the cluster-test fields).
- `11_cvsa_help_integrator.svg` — the SAME per-subject + Stouffer-meta-analysis pattern as the console companion above, but on `main_hybrid_advantage_integ`'s `buf_adv_mi_win_pval` (Hybrid vs MI-only mean **integrator buffer** advantage, restricted to the first `cvsa_influence` seconds of CF) instead of the raw sLDA-probability fusion advantage — directly answers "does CVSA help" on the actual counterfactual control signal rather than the pre-integrator probability. Left panel: per-subject mean advantage bar (coloured by that subject's own within-session significance, `buf_adv_mi_win_pval<0.05`) plus a **GROUP** bar annotated with BOTH the naive cross-subject sign-flip test (`p_bufadvwin_naive_group`, the same low-powered `2^n_subj`-floor test used everywhere else in this figure set, subject is the unit) and the Stouffer meta-analysis (`p_bufadvwin_meta_group`, no small-n floor, combines each subject's own well-powered within-session p-value) — shown side by side and explicitly labelled `naive:`/`meta:` so the two numbers are never confused with each other (before this fix the bar showed only the meta p-value next to an error bar sized for the naive test, which read as an unexplained mismatch). Right panel: the per-subject p-values themselves (`buf_adv_mi_win_pval`, i.e. "was CVSA's help significant for THIS subject alone, within their own trials") plotted directly against the p=0.05 line. Skipped if no subject's `counterfactual_summary.mat` has `buf_adv_mi_win_pval` yet.
- `12_erd_across_subjects.svg` — left panel: per-subject ERD/ERS class discrimination (`|ERD(c1)-ERD(c2)|`, averaged over that subject's bands), grouped by task, with group mean dashed lines — this is the **mean ERD across subjects, per task**; right panel: per-subject CSP-weight correlation (Pearson r), grouped by task, with a group-level sign-flip test on r across subjects — tests whether the neurophysiological grounding shown per-file by `topo_erders.m` holds across the cohort
- `13a_roc_analysis.svg` — **ALL CF frames**. One panel per paradigm present (MI/CVSA/Hybrid): each subject's own re-pooled ROC curve (thin, drawn from a small fixed 7-colour palette × 3 line styles cycled across subjects rather than one unique hue per subject — readable at any cohort size instead of turning into "twenty thousand colours") plus the **cross-subject macro-average curve** (thick, ± SEM band) — each subject's curve is interpolated onto a common FPR grid before averaging, since pooling raw scores directly across subjects would mix different classifiers'/calibrations' probability scales (mirrors the "ALL paradigms" reasoning in `main_roc_analysis`, one level up). The full per-subject legend (bottom-right) is only drawn once, on the last populated panel, to avoid repeating it three times. **4th panel**: three pairwise AUC comparisons — Hybrid vs MI, Hybrid vs CVSA, MI vs CVSA — the SAME paired sign-flip test on the per-subject AUC delta + Cohen's d for all three, subject is the unit. No vs-chance test here (all three AUCs are visibly well above 0.5 already; showing that as an extra statistic on every bar alongside the pairwise brackets was a source of confusion, not clarity).
- `13b_calibration_dprime.svg` — **ALL CF frames**, complements `13a`: AUC only tests whether scores RANK the two classes correctly, not whether the probability VALUES themselves are trustworthy — relevant because both the integrator's threshold-crossing and the Bayesian fusion's `P_CVSA^alpha` term use the actual value, not just its rank. Same per-subject pooled (score,label) frames as `13a`'s ROC. Three panels: (1) Brier score (`mean((P(target)-label)^2)`, LOWER=better) per paradigm with the same 3 pairwise brackets as `13a`'s 4th panel; (2) d-prime (`z(hit rate)-z(false alarm rate)` at P=0.5, HIGHER=better) per paradigm, same layout; (3) reliability diagram — scores binned into 10 equal-width bins, mean predicted P(target) ("confidence") vs observed fraction of true target ("accuracy"), cross-subject mean per paradigm, diagonal = perfect calibration. Skipped together with `13a` if no ROC data is found.
- `14a_roc_window.svg` — **WINDOWED** exact mirror of `13a` (per-paradigm ROC curves + macro-average + pairwise AUC panel), restricted to frames within the first `cvsa_influence` seconds of each trial's CF only — the ONLY window where a fused/Hybrid classifier can actually differ from pure MI (`bayesian_fuse.m` forces `alpha=0`, i.e. `P_fused==P_MI` exactly, from `cvsa_influence` onward, so the ALL-frames version dilutes any early classifier-level advantage with that later, by-construction-identical portion of every trial). Every panel title says "WINDOW ONLY: first Xs of CF" so it's never confused with `13a`. Computed by `main_roc_analysis.m` (`scores_win`/`labels_win`/`window_s` in `roc_summary.mat`'s `roc_by_group`); skipped if not yet regenerated with these fields.
- `14b_calibration_dprime_window.svg` — **WINDOWED** exact mirror of `13b` (Brier, d-prime, reliability diagram), same window as `14a`.
- `15_performance_correlates.svg` — four scatter panels (subject = one point), each with Pearson r + permutation p: (1) counterfactual Hybrid-MI advantage vs real Hybrid-MI advantage — does the simulated advantage predict the real one? (individual-level cross-check of `main_validate_counterfactual`'s pooled comparison); (2) fusion-mechanism advantage vs real Hybrid accuracy; (3) ERD/ERS-CSP grounding (Hybrid) vs real Hybrid accuracy; (4) **CVSA-only decided-trial accuracy vs Hybrid-MI counterfactual advantage** — directly supports the "CVSA helps even when imperfect" thesis: if the advantage stays positive even for subjects with weak CVSA-only accuracy, the benefit isn't contingent on CVSA being good. Console also prints the subject with the weakest CVSA-only accuracy and their advantage explicitly.

**Quick reference — what to expect from each `group_analysis/` figure:**

| Figure | Cosa rappresenta | Risultato sperato | Domanda a cui risponde |
|---|---|---|---|
| `01_real_session_overview.svg` | Accuratezza reale (TIMEOUT=fallimento) e tasso di timeout per paradigma (MI/CVSA/Hybrid), per soggetto + GRAND | Hybrid con accuratezza più alta e timeout più basso di MI e CVSA da soli | Nelle sessioni vere, il sistema Hybrid funziona meglio degli altri due? |
| `02_real_times.svg` | TTH / T-miss / T-timeout medi per paradigma, per soggetto + GRAND (solo media±SEM, nessun test) | Hybrid con TTH più basso (più veloce) | Quanto tempo impiegano in media le tre modalità a produrre un esito? |
| `03_real_accuracy_breakdown.svg` | Accuratezza sulle sole prove decise (TIMEOUT escluso) + falsi positivi + timeout rate, best/worst marcati per paradigma | Hybrid tra i migliori su accuratezza-decise, basso su falsi positivi/timeout | Quando il sistema decide è affidabile, e quanto spesso non decide/decide male? |
| `04_counterfactual_advantage.svg` | Accuratezza SIMULATA (stesso integratore/soglie guidato da Hybrid-fuso/MI-only/CVSA-only sulle stesse prove Hybrid) + delta con test di gruppo | Delta Hybrid−MI e Hybrid−CVSA positivo e significativo | L'algoritmo di fusione in sé (isolato da differenze tra sessioni) produce più HIT? |
| `04b_time_in_correct_zone.svg` | % del tempo di OGNI prova (hit+miss+timeout, 5s interi) in cui il segnale (RAW pre-integratore / INT post-integratore) sta dal lato giusto (>0.5), per stream, con test di gruppo | Hybrid con % più alta di MI-only/CVSA-only, sia RAW che INT | Anche quando l'esito binario non cambia, il segnale Hybrid sta "dalla parte giusta" più a lungo? |
| `05_simulated_accuracy_breakdown.svg` | Stessa scomposizione di `03_real_accuracy_breakdown.svg` ma sui 3 stream counterfactual, con lo stesso pannello scatter accuratezza-vs-timeout | Hybrid migliore anche qui; CVSA-only verosimilmente debole da sola | Come sopra, ma isolato nel mondo simulato/controllato |
| `06_speed_itr.svg` | Delta TTH e delta ITR (bit/min), sia reale che counterfactual, con test di significatività | Hybrid più veloce (TTH negativo) e con ITR più alto, in entrambi i mondi | Hybrid è non solo più accurato ma anche più veloce / un canale di comunicazione migliore? |
| `07_early_late_keepup.svg` | Tra le prove Hybrid-HIT, se MI-only/CVSA-only (simulati sulle stesse prove) avrebbero vinto anch'essi, separando presto (entro `cvsa_influence`) da tardi. x-tick e titolo mostrano anche che % di TUTTE le prove Hybrid (non solo gli HIT) cadono in ciascun gruppo | MI-only "tiene il passo" poco nelle prove veloci, molto in quelle tardive — ma non garantito: l'integratore ha memoria su tutto il trial, quindi anche un HIT "tardivo" può essere in parte merito di una spinta CVSA iniziale | Il contributo della CVSA si concentra nelle prove che si risolvono velocemente? |
| `08_fusion_mechanism.svg` | Vantaggio medio `P_fused-P_MI` per soggetto + rescue vs cost | Vantaggio positivo e significativo, rescue > cost | A livello di probabilità grezza del classificatore, la fusione con la CVSA spinge verso la classe giusta? |
| `09_rescue_hurt_breakdown.svg` | % di sample rescued / hurt / entrambi-giusti / entrambi-sbagliati, su 3 finestre temporali | Rescued > hurt nella finestra attiva, che collassano entrambi a ~0% dopo `cvsa_influence` | Quanto spesso, a livello di singolo campione, la CVSA aiuta vs danneggia vs è ininfluente? |
| `10_fusion_cluster_test.svg` | Test a permutazione su cluster temporali del vantaggio `P_fused-P_MI`, cross-soggetto | Un cluster temporale significativo, concentrato nei primi secondi | L'effetto della fusione è un vero periodo di influenza localizzato nel tempo, o rumore? |
| `11_cvsa_help_integrator.svg` | Stesso tipo di test di `08_fusion_mechanism.svg` ma sul segnale DOPO l'integratore (buffer). Pannello sx: barra per soggetto + barra GROUP che mostra **due** p-value etichettati distintamente — `naive` (sign-flip cross-soggetto, lo stesso test a bassa potenza usato altrove) e `meta` (Stouffer, combina il p-value ben potenziato di ciascun soggetto, nessun limite di risoluzione a basso N). Pannello dx: p-value **per singolo soggetto** (this-subject-only, non di gruppo) contro la linea 0.05 | Vantaggio positivo; meta significativo anche con pochi soggetti; maggioranza dei soggetti sotto 0.05 nel pannello dx | Il contributo della CVSA sopravvive fino al segnale che guida davvero l'esito (non solo alla probabilità grezza)? |
| `12_erd_across_subjects.svg` | Discriminazione ERD/ERS per compito + correlazione coi pesi CSP, cross-soggetto | Correlazione positiva e significativa | Il classificatore sfrutta segnale EEG genuino, o rumore/artefatto? |
| `13a_roc_analysis.svg` | **TUTTI i frame CF.** Curva ROC media cross-soggetto + AUC, per paradigma (colori per soggetto da una tavolozza fissa a 7 colori × 3 stili linea, non un colore unico per soggetto; legenda in basso a dx), più un 4° pannello con 3 confronti a coppie: Hybrid-vs-MI, Hybrid-vs-CVSA, MI-vs-CVSA (nessun test vs caso: già ovvio dal grafico che tutti e tre stanno sopra 0.5) | Hybrid con AUC significativamente più alta di entrambi gli unimodali; MI vs CVSA tipicamente n.s. | Hybrid discrimina meglio dei due unimodali? E i due unimodali sono diversi tra loro? |
| `13b_calibration_dprime.svg` | **TUTTI i frame CF.** Brier score (calibrazione, meno=meglio) e d-prime (sensibilità, più=meglio) per paradigma con gli stessi 3 confronti a coppie di `13a_roc_analysis.svg`, più un diagramma di affidabilità (confidenza predetta vs frequenza osservata, 10 bin) | Hybrid con Brier più basso e d-prime più alto; punti vicini alla diagonale nel diagramma di affidabilità | Le probabilità del classificatore sono anche BEN CALIBRATE (non solo ben ordinate come dice l'AUC)? |
| `14a_roc_window.svg` | **SOLO la finestra `cvsa_influence`** (primi ~2.5s della CF, l'unica finestra dove P_fused può davvero differire da P_MI): mirror esatto di `13a_roc_analysis.svg` (curve ROC + pannello AUC a coppie) ma solo su questi frame | Vantaggio Hybrid più marcato qui che nella versione "tutti i frame", almeno vs CVSA | Il contributo della fusione sul classificatore stesso (non l'integratore) è significativo proprio dove la CVSA è ancora attiva? |
| `14b_calibration_dprime_window.svg` | **SOLO la finestra `cvsa_influence`.** Mirror esatto di `13b_calibration_dprime.svg` (Brier, d-prime, diagramma di affidabilità) ma solo su questi frame | Come sopra, isolato dalla diluizione della parte tardiva del trial (dove fused=MI per costruzione) | Le probabilità del classificatore sono ben calibrate proprio nella finestra dove la CVSA è ancora attiva? |
| `15_performance_correlates.svg` | 4 scatter cross-soggetto: vantaggio counterfactual vs reale, vantaggio fusione vs accuratezza reale, grounding ERD-CSP vs accuratezza reale, qualità CVSA-only vs vantaggio Hybrid-MI | Correlazioni positive; il 4° pannello NON deve andare a zero per i soggetti con CVSA-only debole | Le altre metriche predicono davvero la prestazione reale? Il vantaggio regge anche quando la CVSA da sola è scarsa? |

**Console output**: per-subject real-session tables in **both accuracy conventions** (TIMEOUT=fail, and decided-trials-only + false-positive rate), each with a trailing `GRAND` row, trial-pooled; timeout-rate and times tables; per-subject counterfactual accuracy tables (both conventions); per-paradigm and per-stream **subject ranking tables** (decided-trial accuracy, ties tagged `<-- BEST`/`<-- WORST`); per-subject fusion-mechanism table; per-subject ERD-discrimination/CSP-r table; per-subject ROC-AUC table (MI/CVSA/Hybrid); group-level paired sign-flip tests (mean delta, p-value, stars) **plus paired Cohen's d** (subject is the unit; `negligible`/`small`/`medium`/`large` label for `|d|>0.2`/`0.5`/`0.8`, Cohen 1988) for accuracy — **in both conventions** — (Hybrid vs MI-only/CVSA-only, real and counterfactual) plus a Friedman omnibus test (MI vs CVSA vs Hybrid, both conventions), the fusion mechanism (fusion advantage, rescue-cost — also with Cohen's d), ERD/ERS CSP-weight correlation, classifier AUC vs chance (0.5), and the four cross-subject performance correlates above — all across subjects.

Also saves `group_summary_<tag>.mat` under `<root>/group_analysis/` (`<tag>` = `all`/`well`/`bad`, see `group_tag` above): subjects list; per-subject `real_acc`/`real_acc_dec`/`real_fp_rate`/`real_to_rate`/`real_n_hit`/`real_n_miss`/`real_to_n`/`real_n`/`real_tth`/`real_t_miss`/`real_to_time` matrices `[n_subj x 3]`; the trial-pooled `grand_acc`/`grand_acc_dec`/`grand_fp_rate`/`grand_to_rate`/`grand_to_n`/`grand_n`/`grand_tth`/`grand_t_miss`/`grand_to_time` vectors `[1x3]`; per-subject counterfactual `cf_acc`/`cf_acc_dec`/`cf_fp_rate`/`cf_to_rate`/`cf_n_hit`/`cf_n_miss`/`cf_n_to`/`cf_n`; trial-pooled `grand_cf_acc`/`grand_cf_acc_dec`/`grand_cf_fp_rate`/`grand_cf_to_rate`; per-subject fusion-mechanism vectors; `cluster_time_grid`/`cluster_diff_subj` (per-subject mean `P_fused-P_MI` curve, interpolated onto a common time grid) plus the cross-subject cluster test's `cluster_clusters_group`/`cluster_tstat_group`/`cluster_sig_mask_group`; per-subject ERD-discrimination/CSP-r matrices; per-subject `roc_auc_subj` `[n_subj x 3]`; the cross-subject `roc_fpr_grid`/`roc_mean_tpr`/`roc_sem_tpr`; `idx_best_real`/`idx_worst_real`/`idx_best_cf`/`idx_worst_cf` (1x3 cells of tied subject-index vectors, by decided-trial accuracy); group-level p-values in both conventions (including `p_auc_group`, `p_friedman_real`/`p_friedman_real_dec`) **plus their paired Cohen's d companions** (`d_real_hyb_mi_cohen`/`d_real_hyb_cvs_cohen`/`d_real_mi_cvs_cohen` and `_dec` twins, `d_mi_group_cohen`/`d_cvs_group_cohen` and `_dec` twins, `d_fusadv_group_cohen`, `d_rescue_group_cohen`); `p_bufadvwin_naive_group`/`d_bufadvwin_naive_cohen` (the naive cross-subject sign-flip companion to `p_bufadvwin_meta_group`, shown side by side on `11_cvsa_help_integrator.svg`'s GROUP bar); and the four cross-subject correlate r/p pairs (`r_cf_vs_real`, `r_fus_vs_real`, `r_erd_vs_real`, `r_cvsa_quality`).

### `well_vs_bad_comparison` — unpaired "well" vs "bad" performer comparison

```matlab
well_vs_bad_comparison(gs_well, gs_bad, out_dir, show_figures)
```

Own file, not part of `main_group_analysis.m`. Takes the two `group_summary` structs returned by a `main_group_analysis(..., WELL_SUBJECTS, ..., 'well')` and a `main_group_analysis(..., BAD_SUBJECTS, ..., 'bad')` run (see `batch_group_analysis.m`, the intended caller) and directly compares the two manually-assigned groups — **unpaired**, since "well" and "bad" are disjoint subjects, not the same subjects under different conditions (unlike every other statistical test in this package, which is paired/sign-flip across the SAME subjects). Answers: does the fusion advantage / classifier quality genuinely differ between strong and weak performers — e.g. does Hybrid fusion *compensate* weak performers (bigger advantage in "bad") or only pay off for subjects who are already good (bigger advantage in "well")?

Six panels, each with a hand-rolled **two-sample permutation test** (pooled-label shuffle, 5000 permutations, no Statistics Toolbox) + **unpaired Cohen's d**: real Hybrid decided-trial accuracy, counterfactual Hybrid-MI advantage, counterfactual Hybrid-CVSA advantage, fusion advantage (mean `P_fused-P_MI`), ROC AUC (Hybrid), ROC AUC (MI-only). Per-subject points (jittered) + group mean ± SEM + significance bracket in each panel; console prints the same six comparisons as a table.

Saves `16_well_vs_bad_comparison.svg` under the given `out_dir` (typically `<root>/group_analysis/`, no group tag since it spans both groups).

### `main_browse_gdf` — interactive scrollable viewer

```matlab
main_browse_gdf
```

Runs the full pipeline on one GDF, then shows two synchronised panels in a scrollable figure:
1. sLDA classifier outputs — MI (green), CVSA (orange), fused (purple)
2. Leaky-integrator buffer output per class, with threshold lines; the `init_val` reset point is shown one chunk before each event-781 onset

All GDF markers are drawn as colour-coded vertical lines. Artifact frames are shaded red. Controls: slider + ← → keys to scroll, **Win (s)** field to change window width.

Requirements (all entry points):
- MATLAB R2019+ with the **Signal Processing Toolbox** (`butter`, `filter`)
- **yamlmatlab** on the path — https://github.com/jiri-cigler/yamlmatlab
- **BIOSIG** (`sload`) on the path — https://biosig.sourceforge.io

---

## 2. Pipeline overview

```
GDF (raw EEG + events)
   │
   ├─► apply_processing  ──── features per chunk ──► apply_slda ─► P(c1), P(c2)
   │      CAR + LP+HP per band + ring buffer + CSP + mean power
   │
   ├─► detect_artifacts  ──── artifact flag per chunk ──────────────┐
   │      CAR + EOG bandpass + peaks HP + ring buffer + thresholds  │
   │                                                                ▼
   └─► integrate_signal  ◄──── (P(c1),P(c2)) + artifact flag, per chunk
          per-event-781 leaky binary WTA integrator + (hybrid) Bayesian
          fusion + ROS-style per-class normalisation
                          │
                          ▼
                    plot_trials  (one panel per 781 event)
```

Both the **processing** stream and the **artifact** stream live on the
same chunk-index axis (`k = 1..n_chunks`), with `chunk_size = samplerate / framerate`.
The artifact ring buffer fills earlier than the processing one (250 vs
512 samples in the default config), but by the time the first event 781
fires in the GDF the warmup is long over, so the streams are trivially
aligned by chunk index.

---

## 3. File layout

```
matlab_simulation/
├── run_subject_analysis.m           # BATCH RUNNER: evaluation/ → steps 1-9; calibration/ → topo_erders
├── batch_run_subjects.m             # BATCH RUNNER: hardcoded RECORDINGS_ROOT + SUBJECTS, no GUI → run_subject_analysis per (subject, day)
├── batch_group_analysis.m           # BATCH RUNNER: hardcoded RECORDINGS_ROOT + SUBJECTS (+ optional WELL_SUBJECTS/BAD_SUBJECTS), no GUI → main_group_analysis (per group) + well_vs_bad_comparison + group_topo_erders
├── main_simulate.m                  # function: single-file pipeline replay → per-trial plot + simulate_summary.mat
├── main_session_overview.m          # function: multi-GDF (all paradigms) → 10 figures + session_summary.mat
├── main_trial_dynamics.m            # function: any paradigm, GDF events only → learning/fatigue figures + trial_dynamics_summary.mat
├── main_hybrid_advantage_probs.m    # function: HYBRID GDFs only → mean P + frame acc + CVSA-fusion + cluster test, 6 figures + hybrid_advantage_probs_summary.mat
├── main_hybrid_advantage_integ.m    # function: HYBRID GDFs only → counterfactual hybrid vs MI-only vs CVSA-only, 6 figures + counterfactual_summary.mat
├── main_validate_counterfactual.m   # function: .mat inputs → real session vs offline simulation, 1 figure
├── main_threshold_sweep.m           # function: any paradigm → (th1, th2) HIT/TIMEOUT heatmaps + threshold_sweep_summary.mat
├── main_roc_analysis.m              # function: any paradigm → classifier-level ROC/AUC (pre-integrator) + roc_summary.mat
├── main_group_analysis.m            # function: hardcoded default root (no GUI), optional subjects_filter + group_tag, recursive scan → group-level figures + group_summary_<tag>.mat (one group per call; returns group_summary)
├── well_vs_bad_comparison.m         # function: two group_summary structs ('well'/'bad') → unpaired comparison, 1 figure (16_well_vs_bad_comparison.svg)
├── main_browse_gdf.m                # interactive scrollable viewer (script, not function)
├── io/
│   ├── load_gdf.m             # signal [N x C], header (Label, SampleRate, EVENT.*), basename
│   ├── load_params_yaml.m     # full rosparam-dump struct sibling to the GDF
│   ├── load_csp.m             # extracts processing_fbcsp_<paradigm>.CspCfg.params
│   ├── load_slda.m            # follows slda_node_<paradigm>.path_slda_model and reads sLDACfg.params
│   └── load_hybrid_streams.m  # shared loader: GDF+YAML -> p_mi_aligned/p_cvsa_aligned/art_flags/int_cfg (used by main_hybrid_advantage_probs/_integ)
├── processing/
│   └── apply_processing.m     # FBCSP chunk loop; rescales EVENT.POS/.DUR to chunks
├── artifacts/
│   └── detect_artifacts.m     # artifact flag stream
├── classifier/
│   └── apply_slda.m           # log + sigmoid on every valid feature row
├── integrator/
│   ├── integrate_signal.m     # per-trial integration around each event 781
│   └── bayesian_fuse.m        # hybrid prior fusion (cosine-annealed LOP)
├── plotting/
│   ├── plot_trials.m          # one panel per trial, P(c1) view
│   └── plot_trials_streams.m  # one panel per trial, Hybrid/MI-only/CVSA-only control-signal overlay (main_hybrid_advantage_integ)
└── utils/
    ├── read_yaml.m            # thin yamlmatlab wrapper
    ├── log_step.m             # `[sim] ...` printf used everywhere
    ├── to_vec.m, to_mat.m, to_strcell.m   # coerce yamlmatlab output
    ├── parse_filters_band.m   # "8.0 10.0; 10.0 12.0; ..." -> [n_bands x 2]
    ├── resolve_channels.m     # case-insensitive name -> index lookup against header.Label
    ├── topo_scatter.m         # scalp map at 10-20 positions: scatteredInterpolant bg + head outline
    ├── compute_roc_curve.m    # non-parametric ROC/AUC (no Statistics Toolbox); shared by main_roc_analysis and main_group_analysis
    └── align_streams.m        # (unused at runtime — kept as a reference helper)
```

---

## 4. The companion YAML

The file that `bag_bci` saves next to every GDF is the **rosparam dump**
of the whole launch. It looks like this (truncated):

```yaml
acquisition: { samplerate: 512, framerate: 20, ... }
RingBufferCfg:          { params: { size: 512 } }   # processing ringbuf
RingBufferCfgArtifact:  { params: { size: 256 } }   # artifact ringbuf
CarCfg:                 { params: { EOG_ch_names: [Fp1, Fp2] } }
ArtifactCfg:            { params: { th_hEOG: 7500, th_vEOG: 7500, th_peaks: 15000,
                                    freq_low_EOG: 10, freq_high_EOG: 1,
                                    freq_high_peaks: 1, filterOrder_EOG: 4,
                                    filterOrder_peaks: 4, EOG_ch_names: [Fp1, Fp2] } }
processing_fbcsp_mi:
  do_car: true
  filter_order: 4
  CspCfg:
    params:
      bands: [[8,10],[10,12],[12,14],[8,14],[14,20]]
      selected_channels: [FC5, FC1, C3, ...]
      csp_matrices: [ <band 1>, <band 2>, ... ]   # each [n_components x n_selected]
slda_node_mi:
  path_slda_model: /home/paolo/bci_vr_ws/src/slda_bci/models/mi/slda_mi_test.yaml
integrator:
  paradigm: mi
  classes: [769, 770]
  buffer_size: 40
  increment: 1            # 1 = SOFT, 0 = HARD
  init_val: [0.5, 0.5]
  k_gain: 1.5
  thresholds: [1.0, 0.8]   # raw integrated thresholds, one per class
```

`load_params_yaml.m` reads the whole thing into a single struct. The CSP
matrix is **embedded directly** (`processing_fbcsp_<paradigm>.CspCfg`),
while the sLDA model is **referenced by path** (`slda_node_<paradigm>.path_slda_model`)
and pulled in by `load_slda.m`.

---

## 5. Per-stage details

### `apply_processing.m` — FBCSP

For every chunk (`chunk_size = round(samplerate / framerate)` samples):

1. **CAR** on the chunk if `do_car`: subtract the per-sample mean over all
   channels except the ones in `CarCfg.params.EOG_ch_names`.
2. **Per band**: stateful causal Butterworth low-pass at `band(hi)`, then
   high-pass at `band(lo)`, both order 4, using **MATLAB `butter` ba-form +
   `filter`** (identical to ROS `rtfilter`). Filter state (`zi`) is carried
   across chunks; ICs start at zero. The Python calibration notebook
   (`create_slda.ipynb` Cell 4) must also use `lfilter` ba-form (not
   `sosfilt` SOS) to match this implementation — see `slda_bci/README.md §5`.
3. **Ring buffer push**: NaN-initialised buffer of size `RingBufferCfg.params.size`
   (= `samplerate` = 1 s). New chunk samples replace the oldest. `isfull()`
   ⇔ no NaN anywhere — matches `rosneuro::RingBuffer`.
4. When full, for each band: select `csp.selected_channels` columns from
   the buffer, multiply by `csp_matrices{b}.'`, compute mean power per
   component as `sum(x.^2) / bufsize` (same expression as
   `test_slda.m`).
5. Flatten column-major into a row of length `n_components * n_bands`,
   write into `features(k, :)`.

The returned `header_out` has `EVENT.POS` and `EVENT.DUR` **rescaled to
the chunk timeline** (`round(POS/chunk_size)`, `round(DUR/chunk_size)`),
so downstream stages can read them as chunk indices directly.

Output rows for chunks before the ring buffer fills are `NaN` — sLDA
filters them out and the integrator never reaches them anyway (event
781 is always far past warmup).

### `detect_artifacts.m`

Runs in parallel on the raw signal with its own (smaller) ring buffer
(`RingBufferCfgArtifact.params.size`, default `samplerate/2` = 0.5 s):

1. **CAR** on the chunk using non-EOG channels.
2. EOG path: `chunk_car → LP(freq_low_EOG) → HP(freq_high_EOG)` (stateful).
3. Peaks path: `chunk_car → HP(freq_high_peaks)` (stateful).
4. Push both into circular ring buffers of size `bufsize_art`.
5. When full, compute:
   - `hEOG = max(abs(buf_eog(:,c1) − buf_eog(:,c2)))`
   - `vEOG = max(abs((buf_eog(:,c1) + buf_eog(:,c2))/2 − buf_eog(:,c3)))`
     if there is a third EOG reference; otherwise `max(abs((c1+c2)/2))`.
   - `peaks = max(abs(buf_peaks(:, non_eog)), [], 'all')`.
6. `has_artifact = (hEOG > th_hEOG) || (vEOG > th_vEOG) || (peaks > th_peaks)`.

### `apply_slda.m`

For every row of `features` that is not all-NaN:

```
log_feats = log(features_aligned_to_slda_bands)
score     = log_feats * slda.weights' + slda.intercept   # scalar per row
P(c2)     = 1 / (1 + exp(-score))
P(c1)     = 1 - P(c2)
```

The band-order match (CSP band order vs sLDA band order) is done by
numerical comparison of the `[lo hi]` pairs with 1e-3 tolerance.

### `integrate_signal.m`

For **every event 781** in `header_chunks.EVENT.TYP` (one trial each),
we run a fresh integrator state for the chunks of the CF window. The
CF window is `[POS, POS + DUR]` **inclusive** on both ends (so the
number of integrated chunks is `n_cf = DUR + 1`). Plus **one extra
frame at the very front** — the reset publish (`p_rest = 0.5`).
Nothing is shown outside the CF.

State (per trial, reset at the start):
- `p_prev` — scalar probability of class 1 carried across chunks.
              Initialised to `init_val(1) = 0.5`.
- `frame_count` — chunks elapsed since reset; drives the hybrid CVSA
              prior decay `alpha(t)`.

Per CF chunk `c = start_chunk + j − 1`:

1. **Pick the input** for the binary integrator:
   - `mi`:    `p_in = p_mi(c, 1)` (sLDA P(class 1))
   - `cvsa`:  `p_in = p_cvsa(c, 1)`
   - `hybrid`: Bayesian fusion of MI and CVSA with cosine-annealed prior
     `alpha = 0.5 * (1 + cos(pi * min(t, H) / H))`, where
     `t = frame_count / framerate` and `H = int_cfg.cvsa_influence` (default 3.0 s).
     See `bayesian_fuse.m` for cosine-annealed LOP details. `p_in = fused(1)`.
2. **Leaky binary integrator step**
   (this is `step_integrator` from the validated test, which in turn is
   equivalent to the C++ `Buffer` plugin for the binary case with init=0.5):
   ```
   p_max = max(p_in, 1 - p_in)
   vel   = min(|p_max - 0.5| * 2 * k_gain, 1)        # SOFT mode
   step  = vel / buffer_size
   p_prev = clip( p_prev + sign(p_in - 0.5) * step, 0, 1 )
   ```
3. If `art_flags(c)` is true the buffer **freezes** (no step), exactly
   like the C++ plugin when `has_artifact == true`.
4. Write `integrated(k, :) = [p_prev, 1 - p_prev]`.

Edge frame at the front of the trial:
- **`N_PRE = 1` reset frame** (`t = −1/framerate`): the value
  `resetIntegrator()` publishes to `/raw` and `/normalized` when the 781
  event fires — just `[p_rest, 1 − p_rest] = [0.5, 0.5]`. No input,
  no step. This is the "0.5 sample" you see at the left of every panel.
  It does **not** count toward PASS.

The first integration step happens at `t = 0` (chunk `POS`); the last
integration step happens at `t = (n_cf − 1) / framerate` (chunk `POS + DUR`).

#### ROS-style normalisation

`training_node` ([`feedback_bci_vr/src/Training.cpp`](../../feedback_bci_vr/src/Training.cpp))
normalises **each class independently** with its own threshold after receiving `integrated/raw`:

```
for each class i with threshold thr_i > p_rest:
    slope_i = (1 − p_rest) / (thr_i − p_rest)
    norm_i  = clip( p_rest + (raw_i − p_rest) * slope_i, 0, 1 )
```

So `normalized(:, target_class) >= 1.0` ⇔ `integrated(:, target_class) >= thr(target_class)`.
PASS is checked directly as `integrated(:, target_class) >= thresholds(target_class)`,
matching Training.cpp evaluation mode (`is_target_hit`: `raw[i] >= thresholds[i]`).
The `normalized >= 1.0` form is mathematically equivalent but the raw comparison is canonical.

`integrate_signal` exposes two normalised outputs:

| Field | What it is | Why |
|---|---|---|
| `tr.normalized` | `[n x 2]` — class-1 and class-2 normalised independently, **exactly** as ROS does it | PASS check |
| `tr.normalized_pc1` | `[n x 1]` — single P(c1) curve that uses **both** thresholds piecewise: above `p_rest` uses `thr(c1)`, below `p_rest` uses `1 − thr(c2)` | Visualisation only |

The reason for `normalized_pc1`: with `thresholds = [1.0, 0.8]`,
`normalized(:, 1)` is **identity** in the upper half (`slope = 1`), so
plotting it on top of `integrated(:, 1)` shows the two curves
overlapping. The piecewise `normalized_pc1` view applies the same
`normalize_input` formula symmetrically around `p_rest`, so the
class-2 threshold (`0.8` in P(c2), i.e. `0.2` in P(c1)) reshapes the
lower half of the curve — and you can actually see the stretch.

### `plot_trials.m`

One subplot per 781 trial. Everything is in the **P(class 1) reference
frame** because that is also the only number the integrator consumes.

| Element | Meaning |
|---|---|
| Orange scatter (`tr.raw(:, 1)`) | Raw sLDA P(c1) per frame (for hybrid: post-fusion P(c1)) |
| Solid blue line (`tr.integrated(:, 1)`) | Output of the leaky integrator, P(c1) view |
| Dashed green line (`tr.normalized_pc1`) | Normalised P(c1) using both thresholds piecewise |
| Black dashed horizontal line at `thr(c1)` | Class-1 win boundary in P(c1) (e.g. 1.0) |
| Black dashed horizontal line at `1 − thr(c2)` | Class-2 win boundary in P(c1) (e.g. 0.2) |
| Grey dotted horizontal line at `p_rest = 0.5` | Prior / reset value |
| Two grey dotted vertical lines | CF window boundaries: `t = 0` (first CF chunk = `POS`) and `t = (n_cf − 1)/framerate` (last CF chunk = `POS + DUR`). The reset frame sits at `t < 0` on the left |
| Light-grey patches | Chunks where `art_flags(c) = true` (integrator was frozen) |

Title: `trial k  target=cN (event_code)  PASS|miss`. PASS is computed
**only over the CF window** — the reset frame at `t < 0` doesn't count.

---

## 6. Variables glossary

The names below are reused consistently across all stages.

### Signal / timing

| Name | Meaning |
|---|---|
| `samplerate` (`fs`) | EEG sampling rate from `acquisition.samplerate` (e.g. 512 Hz) |
| `framerate` | Processing rate from `acquisition.framerate` (e.g. 20 Hz) |
| `chunk_size` | `round(samplerate / framerate)` — samples per processing chunk |
| `bufsize` (`bufsize_proc`) | Processing ring-buffer size (`RingBufferCfg.params.size`, 1 s) |
| `bufsize_art` | Artifact ring-buffer size (`RingBufferCfgArtifact.params.size`, 0.5 s) |
| `n_chunks` | `floor(N_samples / chunk_size)` |
| `first_valid_chunk` | `ceil(bufsize / chunk_size)` — first chunk where the ring buffer is full |

### Channels

| Name | Meaning |
|---|---|
| `header.Label` | Channel labels from the GDF (`hdr.Label`) |
| `eog_names` | `CarCfg.params.EOG_ch_names` — excluded from CAR |
| `csp.selected_channels` | Subset of channels CSP operates on (from CSP yaml) |
| `csp_ch` | Indices of `csp.selected_channels` resolved into `header.Label` |

### CSP / FBCSP

| Name | Meaning |
|---|---|
| `csp.bands` | `[n_bands x 2]` of `[lo hi]` pairs |
| `csp.csp_matrices{b}` | `[n_components x n_selected]` per band |
| `csp.n_bands` / `n_components` / `n_selected` | shape metadata |
| `features` | `[n_chunks x (n_components * n_bands)]` raw mean-power (column-major flatten) |

### sLDA

| Name | Meaning |
|---|---|
| `slda.weights` | `[1 x n_features]` |
| `slda.intercept` | scalar |
| `slda.bands` | `[n_bands x 2]` (used to reorder feature columns if CSP and sLDA disagree) |
| `slda.classes` | numeric class labels (e.g. `[769 770]`) |
| `p_mi`, `p_cvsa` | `[n_chunks x 2]` sLDA outputs `[P(c1), P(c2)]` |

### Integrator

| Name | Meaning |
|---|---|
| `int_cfg.classes` | trial-onset event codes, paradigm-specific (e.g. `[769,770]` for MI) |
| `int_cfg.init_val` | `[0.5, 0.5]` — buffer state at reset |
| `int_cfg.thresholds` | `[thr_c1, thr_c2]` — raw integrated thresholds, one per class |
| `int_cfg.buffer_size` | denominator of the per-step delta (40 by default) |
| `int_cfg.k_gain` | multiplier in the SOFT step formula (1.5 by default) |
| `int_cfg.increment` | 0 = HARD step, 1 = SOFT step (binary impl always uses SOFT) |
| `p_rest` | `init_val(1) = 1/n_classes = 0.5` for binary |
| `p_prev` | scalar P(c1) carried across CF chunks within a trial |
| `frame_count` | chunks elapsed since reset (for hybrid `alpha(t)`) |

### Trial struct (one entry per 781 event)

| Field | Shape | Meaning |
|---|---|---|
| `start_chunk` | scalar | `POS` of the 781 event on the chunk timeline |
| `n_cf` | scalar | number of CF chunks = `DUR + 1` (CF spans `[POS, POS+DUR]` inclusive) |
| `n_pre` | scalar | `1` — reset frame at the front (the `p_rest` publish) |
| `onset_code` | scalar | most-recent onset event before this 781 |
| `target_class` | scalar (1 or 2) or NaN | 1-based index into `int_cfg.classes` |
| `raw` | `[n x 2]` | raw / fused sLDA per frame (`NaN` on the reset frame) |
| `integrated` | `[n x 2]` | integrator state per frame (`[0.5, 0.5]` on the reset frame) |
| `normalized` | `[n x 2]` | per-class ROS-style normalisation |
| `normalized_pc1` | `[n x 1]` | P(c1) view using both thresholds piecewise |
| `artifact` | `[n x 1]` logical | whether the integrator was frozen on this frame |
| `pass` | logical | `integrated(:, target_class) >= thresholds(target_class)` somewhere inside the CF window |

`n = n_pre + n_cf = 1 + (DUR + 1)`. Frame `1` is the reset publish at
`t = −1/framerate`; frames `2 … n` are the CF chunks at
`t = 0, 1/framerate, …, (n_cf−1)/framerate`.

---

## 7. What "matches ROS" means here

- Per-band causal Butterworth via `butter()` + stateful `filter()` carrying
  `zi` across chunks, zero ICs at startup — identical to `rtfilter`.
- NaN-initialised ring buffer + `~any(isnan(...))` fill check, matching
  `rosneuro::RingBuffer::isfull()`.
- One chunk per frame (`chunk_size = round(samplerate / framerate)`),
  same per-frame order as ROS (artifact → fbcsp → slda → integrator).
- GDF events delivered to the integrator on the chunk that contains them,
  so event 781 resets the leaky buffer and the hybrid CVSA-prior timer at
  exactly the right sample.
- Bayesian fusion: cosine-annealed LOP with `alpha = 0.5*(1+cos(pi*t/T))` where `T = int_cfg.cvsa_influence` (default 3.0 s). α(0)=1 (full CVSA), α(T/2)=0.5 (equal weight), α(T)=0 (pure MI). No plateau. Symmetric disagreement naturally yields near-uniform output; at α=0 output equals pure MI.
- Leaky binary WTA equivalent to the n-class `Buffer` plugin for the
  symmetric binary case with `init = [0.5, 0.5]`.
- Per-class linear-stretch normalisation, line-by-line port of
  `Integrator::normalize_input`.
