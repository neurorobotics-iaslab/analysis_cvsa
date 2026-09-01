# analysis_bci

Script MATLAB per l'analisi **offline** della pipeline BCI CVSA/MI: calibrazione dei classificatori (GMM + QDA), validazione "pseudo-online" replicando in MATLAB esattamente la logica dei nodi ROS (C++), e aggregazione/plot dei risultati multi-soggetto raccolti nei vari siti sperimentali (Graz, SMC).

Non è un package ROS (non viene compilato da catkin): è una cartella di script/funzioni MATLAB pensati per esecuzione interattiva (usano `uigetfile`, `clear`, aprono più figure).

## Dipendenze

* **MATLAB** con Statistics and Machine Learning Toolbox (`fitgmdist`, `posterior`, `evalclusters`, ...).
* [**BioSig**](http://biosig.sourceforge.net/) per `sload` (lettura file `.gdf`).
* [**yamlmatlab**](https://github.com/ewiger/yamlmatlab) per `ReadYaml` (lettura dei parametri delle run, vedi `utils/loadParameters.m`).
* **EEGLAB** (funzione `topoplot`) per gli script che generano topoplot (`results_graz/grand_avg_fisher_topo.m`, `results_SMC/topoplots.m`, `hilbert/hilbertPower_topoplot.m`).

Diversi script contengono `addpath(...)` con path **assoluti hardcoded** (es. `/home/paolo/Local/Matlab/yamlmatlab`, `/home/paolo/chanlocs39.mat`) — vanno adattati alla propria macchina prima dell'esecuzione.

## Struttura

| Cartella | Contenuto |
|---|---|
| `analysis_approach/` | Calibrazione da zero: a partire dai `.gdf` di calibrazione addestra il **GMM** (intentional/non-intentional control) e prepara il dataset per il **QDA**, per CVSA e per MI (varianti antneuro a 39 canali e gtec a 16 canali). |
| `equal_ros/` | Funzioni "gemelle" del codice C++ online (processing, CAR/CSD + Hilbert, artifact rejection, applicazione GMM/QDA, integrazione temporale). Usate da tutti gli altri script per garantire che l'analisi offline rispecchi esattamente la pipeline ROS. |
| `simulation_ROS/` | `evaluation_pipeline_*`: replay di registrazioni reali (`.gdf` + `.yaml` di configurazione + modelli GMM/QDA salvati) per confrontare il metodo **"my"** (GMM+QDA) col metodo **"traditional"** (solo QDA), calcolando accuratezza, timeout, AUC del GMM, SNR di controllo e mappe soglia-vs-performance. |
| `results_graz/` | Aggregazione e plot multi-soggetto dei risultati del sito Graz (calibrazione GMM, risultati online QDA, Fisher score/topoplot). |
| `results_SMC/` | Stesso tipo di aggregazione/plot per il sito SMC (`evaluation_pipeline_mi_gtec.m`, `grandAverage_plots.m`, `topoplots.m`). |
| `hilbert/` | Analisi esplorative di band power via trasformata di Hilbert (per canale, per ROI, topoplot). |
| `icnic/` | Confronto tra GMM e k-means come metodo per separare campioni "intentional control" (IC) da quelli "non intentional control" (NIC). |
| `utils/` | Funzioni di supporto comuni: metriche (`computeMetrics*`), R² punto-biseriale (`calc_r2_from_data`), ricerca del picco alfa (`analyze_alpha_peak`), caricamento parametri YAML (`loadParameters`), plot ERD (`plot_erd_timecourse`), correzione eventi GDF (`modify_events_gdf`). |

## Dati richiesti

Gli script si aspettano dati organizzati come in `recordings/record_cvsa/<soggetto>/.../` e `recordings/record_mi/<soggetto>/.../`, con:

* **`calibration/gdf/*.gdf`** — registrazioni di calibrazione, accoppiate a **`calibration/parameters/*.yaml`** (config della run, letta da `loadParameters.m`).
* **`evaluation/my/`** e **`evaluation/trad/`** — registrazioni di valutazione online per il metodo GMM+QDA ("my") e per il solo QDA ("trad"), ciascuna con i modelli GMM/QDA effettivamente usati in quella run (`gmm_*.yaml`, `qda_*.yaml`) accanto al `.gdf`/`.mat`.
* I file **`.gdf`** devono contenere gli eventi standard del protocollo (cue classi, `786` fissazione, `781` inizio continuous feedback, `897/898/899` hit/miss/timeout) e un numero di canali coerente con l'headset usato (16 per g.tec, 39 per antneuro — vedi sotto).
* Gli script di aggregazione (`results_graz/`, `results_SMC/`) leggono `.mat` già prodotti da una run precedente di `simulation_ROS`/`analysis_approach` (es. `GMM_Validation_*.mat`, `Metrics_*.mat`, `online_*.mat`, `results_*.mat`), non i `.gdf` direttamente.

## Avvertenze

* Molti path sono **assoluti e hardcoded** per la macchina di sviluppo (`/home/paolo/cvsa/ic_cvsa_ws/...`): da aggiornare se si esegue altrove.
* Gli script sono pensati per esecuzione manuale/interattiva, non come funzioni riutilizzabili in pipeline automatizzate.
