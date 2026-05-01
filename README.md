# C. elegans Volume Analysis

Analysis of volume, speed, and thrashing-rate changes over time in *C. elegans* across multiple NaCl concentrations using WormLab tracking data and reproducible notebooks.

## Repository Layout

```
celegans-volume-analysis/
├── Notebooks/
│   ├── 01_compile_genotype_data.ipynb      # Interactive data compilation
│   ├── 02_celegans_volume_analysis.ipynb   # Primary analysis workflow
│   └── 03_archiving.ipynb                  # Archive complete analysis
├── scripts/
│   ├── compile_genotype_data.py            # Builds compiled CSVs per metric
│   └── run_volume_analysis.py              # Convenience runner
├── data/
│   ├── raw/                                # Raw tracking data by date / condition / video
│   │   └── <date>/
│   │       └── <condition>/               # e.g. 150, 450, 1050
│   │           ├── 0001/                  # Video folder
│   │           │   ├── area.csv
│   │           │   ├── length.csv
│   │           │   ├── position.csv
│   │           │   ├── fit.csv
│   │           │   └── curvature.csv      # or CurvatureMap_Track_N.csv
│   │           └── 0002/
│   ├── metadata/                          # Centralised YAML files, one per condition
│   │   └── <condition>.yaml               # e.g. 150.yaml, 450.yaml
│   ├── processed/                         # Generated compiled CSVs
│   ├── results/                           # Analysis outputs (plots, summary CSVs)
│   └── archive/                           # Timestamped analysis snapshots
├── sample.yaml                            # Template metadata YAML
├── requirements.txt
└── README.md
```

---

## Notebook 01 – Compile Genotype Data

`Notebooks/01_compile_genotype_data.ipynb` builds per-condition compiled CSVs from raw tracking exports.

1. Recursively discovers condition folders that contain numbered video subdirectories.
2. Loads metadata from `data/metadata/<condition>.yaml` (falls back to a YAML inside the condition folder for legacy data).
3. Compiles all video tracks for each metric (area, length, fit, curvature) side-by-side using an **inner merge on `Frame`**, so only frames present in every file are kept. Videos differing by ≤ 2 frames at the end are silently trimmed; larger differences print a warning.
4. Writes `<condition>_compiled_<metric>.csv` files to `data/processed/`.
5. Creates `compilation_summary.csv` (file paths) and `genotype_metadata.csv` (consolidated metadata).

**Supported filename patterns:** legacy (`area.csv`, `curvature.csv`), WormLab batch exports (`BatchExport.csv_Area.csv`), and curvature aliases (`CurvatureMap_Track_N.csv`, `spinecurvature`, etc.).

---

## Notebook 02 – Volume & Motility Analysis Workflow

`Notebooks/02_celegans_volume_analysis.ipynb` runs the full multi-modal analysis. Cells run top-to-bottom.

### Volume analysis

| Cell | What it does |
|------|-------------|
| Multi-genotype pipeline | Loads compiled area/length/fit CSVs, computes volume ($V = \pi \cdot \text{Area}^2 / (4 \cdot \text{Length})$), applies quality filters, normalizes to %V0, smooths, and plots mean ± SE curves for each concentration. |
| Per-track direct measurements | Computes V0 (first 28 frames), V_final (last 28 frames), max ΔV, and time to 50% volume change without curve fitting. |
| Statistical analysis – time to 50% | Bar chart + individual points, one-way ANOVA with significance annotation. |
| Statistical analysis – max ΔV | Same layout for maximum volume change. |

**Volume track filters:** fit ≥ 0.9 · initial volume ≥ 1.5 µm³ · duration ≥ 400 s valid measurements.

### Speed analysis

Speed is computed from WormLab position tracks (14 fps):

- **Lookback window:** 28 frames (2 s) — displacement from frame *i − 28* to frame *i* divided by elapsed time.
- The first 28 frames are NaN by construction, preventing start-of-recording spikes.
- Frames with fit quality < 0.9 are set to NaN.
- No additional smoothing; the 2-second lookback provides inherent temporal averaging.
- **V0 per track:** median speed during the first 10 s (robust to outlier frames).
- **Normalized speed:** % of per-track V0.

Outputs: per-track long-form CSV, per-condition summary CSV, absolute speed plot (µm/s), and normalized speed plot (% V0).

### Thrashing-rate analysis

Thrashing rate is estimated from compiled curvature maps:

- Midbody segments (30–70% of body) are averaged to form a single bending trace per track.
- Fit masking (< 0.9), 14-frame rolling median smooth, and baseline subtraction are applied.
- Bending events (peaks and troughs) are detected with `scipy.signal.find_peaks`.
- **Rolling rate** is computed in a 5-second window: `bend count / (2 × window)` in Hz.
- **Time to thrashing suppression:** first time the rolling rate stays ≤ 25% of the track's initial rate (first 1 s baseline) for ≥ 5 s.

Outputs: per-track long-form CSV, per-condition summary CSV, composite mean ± SE plot, multipanel plot, and suppression summary table.

### Cross-modal timing comparison

Compares when each modality reaches suppression across volume loss, motility loss, and thrashing loss:

- Two suppression thresholds: **50%** and **25%** of initial value.
- Sustained criterion: signal must remain at or below threshold for ≥ 5 s.
- Motility and thrashing baselines use the first **10 s** of the recording.
- Per-track event tables are outer-merged and pairwise Δt values computed.
- Bootstrap resampling (within concentration) produces 95% CIs on earliest-event probabilities.

---

## Notebook 03 – Archiving

`Notebooks/03_archiving.ipynb` creates timestamped archives of complete analyses:

1. Packages raw data, processed data, notebooks, scripts, and configuration files.
2. Creates `data/archive/celegans_analysis_YYYYMMDD_HHMMSS/`.
3. Generates a `MANIFEST.json` with metadata for all archived files.

---

## Key Outputs

| Path | Contents |
|------|----------|
| `data/processed/<condition>_compiled_<metric>.csv` | Compiled tracking data per condition and metric |
| `data/processed/compilation_summary.csv` | Index of all compiled files |
| `data/processed/genotype_metadata.csv` | Consolidated metadata from all YAML files |
| `data/results/<condition>_speed_tracks_long.csv` | Per-track speed time series |
| `data/results/<condition>_speed_summary.csv` | Mean ± SE speed per condition |
| `data/results/<condition>_thrashing_tracks_long.csv` | Per-track thrashing-rate time series |
| `data/results/<condition>_thrashing_summary.csv` | Mean ± SE thrashing rate per condition |
| `data/results/speed_mean_se.png` | Absolute speed plot |
| `data/results/speed_normalised_mean_se.png` | Normalised speed plot |
| `data/results/thrashing_rate_composite_mean_se.png` | Thrashing rate composite plot |
| `data/archive/celegans_analysis_YYYYMMDD_HHMMSS/` | Timestamped complete analysis archive |

---

## Getting Started

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
2. **Organise raw data** under `data/raw/<date>/<condition>/0001/`, `0002/`, etc. Add a `data/metadata/<condition>.yaml` for each condition (see `sample.yaml` for the template).
3. **Compile:** run `Notebooks/01_compile_genotype_data.ipynb` (or `python scripts/compile_genotype_data.py`) to generate compiled CSVs in `data/processed/`.
4. **Analyse:** run `Notebooks/02_celegans_volume_analysis.ipynb` top-to-bottom. Results are written to `data/results/`.

---

## Requirements

- Python 3.9+
- pandas, numpy, matplotlib, seaborn, scipy, pyyaml, jupyter (see `requirements.txt`)

## Citation

If you use this workflow in your research, please cite appropriately.
