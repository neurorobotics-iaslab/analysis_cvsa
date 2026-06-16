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

Six entry points, all self-contained (pick GDFs via file dialogs).
`main_hybrid_advantage_probs` and `main_hybrid_advantage_integ` **require hybrid GDF files** — they analyse the Bayesian fusion of MI and CVSA and will not produce meaningful results on MI-only or CVSA-only recordings.

### `main_simulate` — single-file, visual inspection

```matlab
cd /home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation
main_simulate
```

Single GDF → companion YAML → full pipeline → one figure with all CF trials.
Use this to inspect a specific recording in detail. The figure is saved as SVG
(full-screen size) under `<gdf_dir>/analysis_results/trial_simulation/`; the
`SHOW_FIGURES` flag at the top of the script controls whether it is also shown
on screen.

### `main_session_overview` — multi-file cross-paradigm summary

```matlab
main_session_overview
```

Multi-select GDFs of any paradigm (MI, CVSA, Hybrid mixed). Paradigm is inferred
from the filename (`hybrid` > `cvsa` > `mi`). Files are sorted MI → CVSA → Hybrid
with group separators in the figures.

**Seven figures:**

| Figure | Content | Source |
|---|---|---|
| Fig 1 | Trial accuracy — two subplots: `897/(897+898+899)` and `897/(897+898)` (no timeout). Per-file bar + per-paradigm dashed mean line. | Real GDF events |
| Fig 2 | Time metrics — TTH (HIT trials, mean ± std + dots) and time-to-MISS (MISS trials). | Real GDF events |
| Fig 3 | Sample accuracy (3 subplots): MI classifier (full CF), CVSA classifier (first 3 s), Hybrid fused (full CF). Each subplot shows all / HIT / MISS breakdown plus per-class split. | Offline simulation |
| Fig 4 | **CSP filter weights** (one row per CSP type: MI / CVSA). Full-scalp topoplots via `topo_map(bg_zero=true)` — all standard 10-20 channels shown (selected = black-bold, rest = gray). Col 1: class-1 channel importance (% of total filter weight). Col 2: class-2 channel importance. Col 3: per-channel selectivity `(w1−w2)/(w1+w2)` — red = class-1-only, blue = class-2-only. Col 4: **stacked band bar** — bar height = total |CSP filter| energy in that band (% of all bands); stack split = which class dominates it (blue = class-1 components, orange = class-2 components). Skipped if no companion YAML. | CSP model (from YAML) |
| Fig 5 | **sLDA-weighted feature importance** (one row per type: MI / CVSA). Restricted to sLDA-selected features, weighted by `|sLDA coefficient|`. Col 1: channel topoplot — `Σ_k |coef_k| · |CSP_filter_k(ch)|` per channel (full 10-20 grid shown). Col 2: band bar — same sum per frequency band. Col 3: feature selection heatmap (CSP component × band) — color = `|coef|` for selected features, gray = excluded by sLDA. Skipped if no companion YAML. | CSP + sLDA models (from YAML) |
| Fig 6 | **ERD/ERS ↔ CSP weight correlation** (one panel per CSP type × frequency band). Scatter of per-channel ERD/ERS class discrimination `|ERD(c1) − ERD(c2)|` during CF vs normalised CSP filter weight `Σ|W|/max`, Pearson r in title, linear fit line. Tests the claim: channels heavily weighted by the CSP spatial filter also show stronger neurophysiological class separation. Skipped if no companion YAML. | Raw GDF signal + CSP model |
| Fig 7 | **Hemisphere ROI ERD/ERS time-courses** (one row per CSP type × non-empty ROI). CSP-selected channels grouped by 10-20 label suffix into Left (odd digit), Right (even digit), Midline (z). ERD/ERS band-averaged; baseline = cue period [−1.5s, 0s]. Left col: per-class traces ± across-channel SEM. Right col: lateralization c1(t)−c2(t). Gray solid = CF onset, blue dashed = CVSA influence window (3s), gray dotted = cue onset (−1.5s). Skipped if no companion YAML. | Raw GDF signal + CSP model |

Console output per file: a per-trial outcome table (HIT/MISS/TIMEOUT + time for that
event), hit rate, sample accuracy (all / HIT / MISS / per-class), integrator
parameters, top-3 CSP channels per class (`[CSP-MI] top ch c769: C3>FC5>CP5 | c770: C4>FC6>CP6`),
sLDA feature count (`[sLDA-MI] 12/24 features selected by sLDA`), and
ERD trial counts (`[ERD-MI] c1=N c2=N trials`).
Set `VERBOSE_DIAGNOSTIC = true` (default `false`) for an additional per-trial diagnostic
table — very verbose, off by default.

The seven figures are saved as SVG (full-screen size) into `<gdf_dir>/analysis_results/overview/`:
`overview_trial_accuracy.svg`, `overview_time_metrics.svg`, `overview_sample_accuracy.svg`,
`overview_csp_importance.svg`, `overview_slda_importance.svg`, `overview_erd_csp_corr.svg`,
and `overview_roi_timecourse.svg`. The `SHOW_FIGURES` flag at the top of the script controls
whether they are also shown on screen.

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

Does **not** save `.mat` files — it is a quick overview, not an archiving tool.

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

**Five figures:**

| Figure | Content |
|---|---|
| Fig 1 | **Mean P(target) scatter** — one point per trial per stream (5 panels: MI / CVSA / CVSA-inf / fused / buffer). Color = outcome (green HIT / orange MISS / red TO); marker = cued class. Dashed lines = per-outcome group means. |
| Fig 2 | **Frame accuracy scatter** — same layout as Fig 1, y-axis = fraction frames where P > 0.5. |
| Fig 3 | **Does CVSA fusion help? (binary view)** — 2×2 layout: (1,1) **frame-level rescue/hurt bar chart** — counts, over all trials within the CVSA-influence window, of CF frames where `P_MI(target)` and `P_fused(target)` agree/disagree on being `>0.5`: *rescued* (MI wrong → fused correct), *hurt* (MI correct → fused wrong), *both correct*, *both wrong*. Title shows the net effect (`rescued − hurt`) and a verdict ("CVSA helped"/"CVSA hurt"/"neutral"); (1,2) per-trial **CVSA-fusion advantage** = mean(P_fused(target) − P_MI(target)) over the CVSA-influence window — positive = fusion pushed toward the target class more than MI alone, negative = fusion pulled away; (2,1)/(2,2) mean ± SEM of P_MI(target) and P_fused(target) over CF time, one panel per cued class (each truncated to that class's longest trial, so the SEM band renders for both classes), with a vertical line marking where the CVSA prior fades to zero (`cvsa_influence`). |
| Fig 4 | **CF trial duration & CVSA effect on HIT trials, by hit timing** — two panels: (1) per-trial CF duration scatter colored by outcome, with per-outcome dashed mean lines and a dotted line at `cvsa_influence` (split out from Fig 3 into its own figure); (2) **HIT trials only** (MISS/TIMEOUT excluded) — mean `P_MI(target) → P_fused(target)` slopegraph (± SEM, delta annotated), pooling CF frames within the first `cvsa_influence` seconds, split into trials that reached HIT *within* `cvsa_influence` (`duration < cvsa_influence`) vs *after* it (`duration ≥ cvsa_influence`) — tests whether fusion helped fast hits cross threshold sooner (positive delta) vs corrected/hindered MI early in slower hits. |
| Fig 5 | **Does CVSA help, broken down by MI/CVSA agreement** — every CF frame within the CVSA-influence window is classified into one of 4 categories by whether `P_MI(target)>0.5` and `P_CVSA(target)>0.5`: *agree-correct* (both right), *MI correct, CVSA wrong*, *MI wrong, CVSA correct* (the "rescue" case), *agree-wrong* (both wrong). Three panels: (1) scatter of `P_MI(target)` vs `P_fused(target)` per frame, colored by category, with the `y=x` diagonal and crosshairs at 0.5; (2) slopegraph of mean `P_MI(target) → P_fused(target)` per category (± SEM), with the mean delta annotated above each pair — this is the "small improvement" view; (3) for the two disagreement categories only, `P_fused(target) − P_MI(target)` vs the cosine-annealed CVSA weight `α(t)`, showing how the rescue/cost effect scales with `α`. |

**Console output** (per trial, then summaries):
- Per-trial table: `#, cue, result, mMI, aMI, mCV, aCV, mCVi, aCVi, mFus, aFus, mBuf, aBuf`
- Per-class summary: mean of all metrics split by cued class
- Per-outcome summary: mean of all metrics split by HIT/MISS/TIMEOUT
- CVSA-fusion advantage summary: count of trials where fusion helped (>0) vs hurt (<0), and the mean delta
- CVSA-fusion frame-level effect (binary): rescued/hurt/both-correct/both-wrong frame counts within the CVSA-influence window, plus the net effect and verdict
- HIT trials by hit timing: trial/frame counts, mean `P_MI -> P_fused`, and delta for HIT-within-`cvsa_influence` vs HIT-after-`cvsa_influence` groups
- CVSA vs MI agreement breakdown: per-category frame counts, mean `P_MI -> P_fused`, and delta; plus rescue effect (MI wrong, CVSA correct) and cost effect (MI correct, CVSA wrong)

**Interpreting mean P:**

| Value | Interpretation |
|---|---|
| > 0.5 + HIT | Classifier strong, integrator fast |
| > 0.5 + MISS | Classifier correct direction but integrator too slow → tune `k_gain`/`buffer_size` |
| ≈ 0.5 + MISS | Near-chance; integration goes nowhere |
| < 0.5 + MISS | Classifier points wrong direction → genuine classifier failure |

The five figures are saved as SVG (full-screen size) into `<gdf_dir>/analysis_results/advantage_hybrid/`;
the `SHOW_FIGURES` flag at the top of the script controls whether they are also shown on screen.
No `.mat` files are saved.

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

**Console output**: per-trial table; aggregate counts/accuracy/mean times; "who hits?" breakdown; per-class accuracy; **sign-flip permutation test** on mean buffer(target) advantage (Hybrid vs MI-only, Hybrid vs CVSA-only); **Cohen's d** (paired effect size) on buffer advantage; **bootstrap 95% CI** on per-stream accuracy and on accuracy deltas.

**Statistics (no Stats Toolbox required)**: the sign-flip permutation test (`sign_flip_test` local helper, 2000 permutations for global tests, 500 for per-frame temporal test) and the bootstrap CI (percentile method, 2000 resamples) are both implemented without MATLAB's Statistics and Machine Learning Toolbox. The temporal significance figure is explicitly labelled as pointwise/uncorrected — for small n (40 trials) it is exploratory; sustained runs of low p-values are more meaningful than isolated frames.

All figures are saved as SVG (full-screen size) under
`<gdf_dir>/analysis_results/hybrid_vs_unimodal/`; the `SHOW_FIGURES` flag at the
top of the script controls whether they are also shown on screen. No `.mat` files
are saved.

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
├── main_simulate.m                  # single-file: GUI → pipeline → per-trial plot
├── main_session_overview.m          # multi-GDF (all paradigms): trial acc + TTH + sample acc + CSP + sLDA analysis, 5 figures
├── main_hybrid_advantage_probs.m    # HYBRID GDFs only: mean P + frame acc + CVSA-fusion advantage, 5 figures
├── main_hybrid_advantage_integ.m    # HYBRID GDFs only: counterfactual hybrid vs MI-only vs CVSA-only, per-trial + aggregate figures
├── main_browse_gdf.m                # interactive scrollable viewer: classifier probabilities + integrator signal
├── io/
│   ├── load_gdf.m             # signal [N x C], header (Label, SampleRate, EVENT.*), basename
│   ├── load_params_yaml.m     # full rosparam-dump struct sibling to the GDF
│   ├── load_csp.m             # extracts processing_fbcsp_<paradigm>.CspCfg.params
│   └── load_slda.m            # follows slda_node_<paradigm>.path_slda_model and reads sLDACfg.params
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
