# C. elegans Volume and Locomotion Analysis

Reproducible analysis of osmotic-stress responses in *C. elegans* using WormLab tracking data.

The current workflow focuses on three related analyses:

1. **Body-volume response**
2. **Locomotor-speed response**
3. **Volume–speed response kinetics**, including normalized response progress and T50 timing

Experimental conditions are defined through YAML metadata and ordered using `condition_order`, rather than hard-coded condition lists.

---

## Repository Layout

```text
celegans-volume-analysis/

├── Notebooks/
│   ├── 01_compile_genotype_data.ipynb
│   ├── 02_celegans_volume_analysis.ipynb
│   └── 03_archiving.ipynb
│
├── data/
│   ├── raw/
│   │   └── <experiment_id>/
│   │       └── wormlab_exports/
│   │           └── <strain>/
│   │               └── <condition>/
│   │                   ├── 0001/
│   │                   ├── 0002/
│   │                   └── ...
│   │
│   ├── metadata/
│   │   └── *.yaml
│   │
│   ├── processed/
│   │   ├── genotype_metadata.csv
│   │   ├── compilation_summary.csv
│   │   └── ...
│   │
│   ├── results/
│   │   └── analysis figures and derived tables
│   │
│   └── archive/
│       └── timestamped analysis snapshots
│
├── sample.yaml
├── requirements.txt
└── README.md
```

For example:

```text
data/raw/
└── 260311/
    └── wormlab_exports/
        └── N2/
            ├── 21/
            │   ├── 0001/
            │   ├── 0002/
            │   └── ...
            ├── 150/
            ├── 250/
            ├── 350/
            ├── 450/
            ├── 550/
            ├── 650/
            ├── 750/
            ├── 850/
            ├── 950/
            └── 1050/
```

Local Python environments such as `.venv/` or `venv/` are development environments and are not part of the analysis data structure.

---

# Analysis Workflow

The workflow is organized into three notebooks that should normally be run in order:

```text
01_compile_genotype_data.ipynb
        ↓
02_celegans_volume_analysis.ipynb
        ↓
03_archiving.ipynb
```

---

# Notebook 01 – Compile Experimental Data

`Notebooks/01_compile_genotype_data.ipynb`

This notebook compiles WormLab exports and experimental metadata into standardized files for downstream analysis.

## Metadata

Experimental metadata are stored as YAML files under:

```text
data/metadata/
```

Each metadata file describes an experimental condition and can contain:

- experiment ID
- assay and analysis dates
- genotype
- strain
- genetic background
- animal age
- treatment information
- NaCl concentration
- assay conditions
- arena information
- imaging parameters
- WormLab configuration
- experimental notes
- optional per-video metadata

A template is provided in:

```text
sample.yaml
```

## Condition ordering and control definition

Experimental groups are explicitly ordered using:

```yaml
condition_order: <integer>
```

The condition with:

```yaml
condition_order: 0
```

is defined as the experimental control.

Positive integer values define the order of the experimental conditions.

This allows downstream analysis and plotting to derive condition order and control identity directly from metadata rather than hard-coding a particular NaCl concentration.

Notebook 01 validates `condition_order` so missing, non-numeric, non-integer, non-finite, or duplicate values are detected before downstream analysis.

## Compilation outputs

Processed files are written under:

```text
data/processed/
```

including:

```text
genotype_metadata.csv
compilation_summary.csv
```

along with compiled WormLab measurements used by the volume analysis.

`genotype_metadata.csv` serves as the authoritative compiled metadata table for downstream analysis.

---

# Notebook 02 – Volume and Locomotion Analysis

`Notebooks/02_celegans_volume_analysis.ipynb`

This is the primary analysis notebook.

The current workflow contains three major analysis sections:

1. Volume response
2. Locomotor-speed response
3. Normalized Volume–Speed kinetic comparison

Cells are intended to run from top to bottom.

A single global figure-style section controls plotting conventions throughout the notebook.

---

# Volume Analysis

## Volume estimation

Worm volume is estimated from WormLab body area and length measurements using:

\[
V = \frac{\pi A^2}{4L}
\]

where:

- \(A\) = measured body area
- \(L\) = measured body length
- \(V\) = estimated body volume

Quality-control criteria are applied before downstream volume measurements are calculated.

---

## Volume time course

Volume trajectories are normalized to each worm's initial state and summarized within each osmotic condition.

Time-course figures show the concentration- and time-dependent volume response following osmotic stress.

Experimental conditions are displayed according to their metadata-defined `condition_order`.

---

## Final volume endpoint

The primary volume endpoint is calculated from each worm's final 20 s of valid measurements.

Endpoint figures display:

- individual worms
- group mean
- SEM
- sample size
- statistical comparisons with control

Each worm is treated as an individual biological observation.

For comparisons with the metadata-defined control:

- Welch's two-sample t-test is used
- each experimental condition is compared with control
- Holm correction is applied across the control comparisons

---

# Volume Transition Kinetics

For worms producing a measurable volume response, transition kinetics are estimated from the individual trajectory.

The analysis identifies:

- initial plateau
- final plateau
- response amplitude
- V20
- V50
- V80
- T20
- T50
- T80
- transition duration
- transition rate

A control-derived response-amplitude threshold is used to distinguish measurable responses from fluctuations too small to support reliable kinetic normalization.

Persistent crossing criteria are used so kinetic measurements are not determined by a single noisy observation.

## Primary timing metric

The principal timing phenotype is:

```text
T50
```

T50 represents the time at which an individual worm reaches 50% of its measured response.

Nonresponders and tracks without a valid crossing remain missing rather than being assigned artificial timing values.

---

# Locomotor Speed Analysis

Locomotor speed is calculated directly from WormLab X–Y position data.

## Raw-data discovery

Speed measurements are obtained from the raw WormLab hierarchy:

```text
data/raw/
└── <experiment_id>/
    └── wormlab_exports/
        └── <strain>/
            └── <condition>/
                └── <video>/
```

This hierarchy provides technical identifiers including:

- experiment ID
- strain
- condition
- video ID

Condition order and control identity are obtained from:

```text
data/processed/genotype_metadata.csv
```

The speed analysis does not require per-video YAML metadata matching.

---

## Primary speed definition

The primary speed parameters are:

```text
Frame rate:                    14 fps
Lookback:                      14 frames
Nominal lookback duration:     1 s
Fit threshold:                 0.9
Minimum valid-window fraction: 0.80
```

Speed is calculated from endpoint displacement:

\[
\text{Speed}
=
\frac{\sqrt{(X_t-X_{t-14})^2+(Y_t-Y_{t-14})^2}}
{\Delta t}
\]

where \(\Delta t\) is the actual elapsed time between endpoints.

## Speed validity requirements

A speed observation requires:

- finite X and Y coordinates at both endpoints
- Fit ≥ 0.9 at both endpoints
- an endpoint frame difference of exactly 14 frames
- positive elapsed time
- sufficient valid measurements within the expected inclusive 15-frame window

At least 80% of the expected measurements within the window must satisfy the measurement-validity criteria.

Missing positions and poor-fit endpoints are not interpolated or replaced.

---

## Speed time course

The primary 1-s speed measurements are retained without additional analytical smoothing.

For visualization and time-course summarization:

```text
1-s speed measurements
        ↓
mean within each worm × 5-s time bin
        ↓
condition mean ± SEM across worms
```

This reduces frame-level visual noise while preserving the worm as the statistical unit.

The 5-s binning is a visualization/summarization step and does not redefine the underlying speed measurement.

---

## Speed sensitivity analysis

Important speed-analysis parameters are evaluated through sensitivity analysis.

Alternative lookback windows include:

```text
7 frames   = 0.5 s
14 frames  = 1 s
28 frames  = 2 s
56 frames  = 4 s
```

Alternative minimum valid-window fractions include:

```text
0.60
0.80
1.00
```

The primary analysis uses:

```text
14-frame lookback
0.80 valid-window fraction
```

Sensitivity analysis is used to determine whether major biological conclusions depend strongly on these choices.

---

## Final speed endpoint

The primary speed endpoint is the mean absolute speed during each worm's final 20 s of valid measurements.

Endpoint figures display:

- individual worms
- group mean ± SEM
- sample size
- statistical significance relative to control

Statistics parallel the volume endpoint analysis:

- Welch's two-sample t-test
- each experimental condition versus the metadata-defined control
- Holm correction across comparisons

Absolute speed in µm/s is the primary locomotor measurement.

---

# Volume–Speed Response Kinetics

Volume and speed have different units and response amplitudes.

To compare their **response timing** directly, each phenotype is independently transformed onto a common response-progress scale.

## Response-progress normalization

For each individual worm and phenotype:

\[
\text{Response progress}(t)
=
100
\times
\frac{X_{\mathrm{initial}}-X(t)}
{X_{\mathrm{initial}}-X_{\mathrm{final}}}
\]

Therefore:

```text
0%   = initial state
50%  = half of the measured response completed
100% = final measured response
```

This is not percent-of-baseline normalization.

For example, a volume trajectory changing from:

```text
1.00 → 0.70
```

is transformed from:

```text
0% → 100% response progress
```

even though its absolute volume decrease is only 30%.

Similarly, a speed trajectory changing from:

```text
250 µm/s → 5 µm/s
```

is independently transformed from:

```text
0% → 100% response progress
```

The resulting trajectories therefore compare **response kinetics rather than absolute response magnitude**.

---

# Individual-First Normalization

Response-progress normalization is performed before group averaging.

The analysis order is:

```text
individual trajectory
        ↓
initial and final values
        ↓
response eligibility
        ↓
individual response-progress normalization
        ↓
5-s time binning
        ↓
condition mean ± SEM
```

Group-average trajectories are never normalized after averaging.

Normalized values are not analytically clipped to 0–100%, allowing transient excursions outside the expected range to remain visible for QC.

---

# Response Eligibility

Response-progress normalization becomes unstable when the initial-to-final response amplitude is very small.

Eligibility is therefore determined before normalized trajectories are summarized.

Eligibility is based on explicit response-amplitude and data-quality criteria rather than excluding tracks because their normalized trajectories appear noisy.

## Volume eligibility

Volume eligibility uses the control-derived response-amplitude threshold established by the volume transition-kinetics analysis.

## Speed eligibility

Speed eligibility is derived independently from control speed responses because volume and speed have different units and noise characteristics.

Technical eligibility also requires adequate valid measurements in the initial and final reference windows.

Tracks failing eligibility remain documented in audit outputs but do not contribute to normalized response-progress summaries.

---

# Conditions Used for Volume–Speed Comparison

The normalized kinetic comparison focuses on osmotic conditions with measurable responses:

```text
350
450
550
650
750
850
950
1050 mM NaCl
```

Lower-stress conditions:

```text
21
150
250 mM NaCl
```

remain part of the absolute volume and speed analyses but are not included in the 0–100% response-progress comparison.

This prevents small low-stress fluctuations from being artificially expanded onto a full-response scale.

---

# Volume–Speed Co-plots

For each included NaCl concentration, normalized Volume and Speed response-progress trajectories are plotted together.

Each panel displays:

- Volume response progress
- Speed response progress
- condition mean ± SEM
- 0% and 100% response references
- eligible sample sizes

These plots provide a direct visualization of the temporal relationship between locomotor inhibition and body-volume loss.

---

# T50 Analysis

For both Volume and Speed:

```text
T50 = time at which response progress reaches 50%
```

Persistent crossing and interpolation logic are used so T50 is not determined by a single noisy observation.

Invalid crossings and non-crossers remain missing.

Volume and Speed T50 values can then be compared across NaCl concentrations.

---

# ΔT50 Analysis

When Volume and Speed measurements can be reliably paired for the same individual worm:

\[
\Delta T50
=
T50_{\mathrm{Volume}}
-
T50_{\mathrm{Speed}}
\]

Interpretation:

```text
ΔT50 > 0
Speed reaches 50% response before Volume.

ΔT50 ≈ 0
Volume and Speed reach 50% response at similar times.

ΔT50 < 0
Volume reaches 50% response before Speed.
```

This provides a direct measure of the temporal separation between locomotor inhibition and volume loss.

---

# Data Structure and Statistical Units

The workflow preserves experimental and technical identifiers wherever available, including:

- experiment ID
- strain
- condition
- condition order
- video ID
- track ID
- worm ID

Individual worms are the biological observations used for endpoint and kinetic summaries.

Video identity is retained so potential within-video clustering or technical effects can be evaluated.

---

# Figures and Prism Exports

Notebook 02 uses a single global plotting-style configuration.

Where appropriate, figures are saved using the shared:

```python
save_figure_with_prism(...)
```

workflow.

This produces figure files together with Prism-ready CSV tables constructed directly from the underlying analysis variables.

Scientific data are not reconstructed from Matplotlib artists.

Typical outputs include:

- PDF figures
- PNG figures
- Prism-ready CSV files
- endpoint statistics
- response-eligibility audit tables
- individual T50 measurements
- paired ΔT50 measurements

Results are written under:

```text
data/results/
```

---

# Quality Control

The workflow retains QC information rather than silently removing problematic observations.

Current QC includes:

- Fit-based measurement validity
- speed valid-window coverage
- initial/final reference-window coverage
- response-amplitude eligibility
- response-normalization audit tables
- valid T50 crossing status
- failed-crossing reasons
- sample sizes by condition
- retained-observation fractions
- speed lookback sensitivity
- speed valid-fraction sensitivity

Excluded observations remain identifiable through audit outputs whenever possible.

---

# Notebook 03 – Archiving

`Notebooks/03_archiving.ipynb`

The archiving notebook creates timestamped snapshots of the analysis so the data, metadata, notebooks, and outputs associated with a particular analysis state can be preserved.

Archives are written under:

```text
data/archive/
```

The archive workflow can preserve:

- raw data
- processed data
- metadata
- notebooks
- analysis outputs
- configuration files
- file manifests

---

# Key Project Locations

| Location | Contents |
|---|---|
| `data/raw/` | Original WormLab exports |
| `data/metadata/` | Experimental YAML metadata |
| `data/processed/` | Compiled and standardized data |
| `data/processed/genotype_metadata.csv` | Consolidated experimental metadata |
| `data/processed/compilation_summary.csv` | Compilation records |
| `data/results/` | Figures, statistics, QC, and derived analysis tables |
| `data/archive/` | Timestamped analysis archives |
| `Notebooks/01_compile_genotype_data.ipynb` | Data and metadata compilation |
| `Notebooks/02_celegans_volume_analysis.ipynb` | Primary scientific analysis |
| `Notebooks/03_archiving.ipynb` | Analysis archiving |
| `sample.yaml` | Metadata template |
| `requirements.txt` | Python dependencies |

---

# Getting Started

## 1. Create a Python environment

For example:

```bash
python -m venv .venv
```

Activate the environment using the appropriate command for your operating system.

## 2. Install dependencies

```bash
pip install -r requirements.txt
```

## 3. Add experimental metadata

Create YAML metadata files under:

```text
data/metadata/
```

Use:

```text
sample.yaml
```

as the template.

Each experimental condition must have a valid `condition_order`.

Exactly one control condition within an experiment should have:

```yaml
condition_order: 0
```

## 4. Organize raw data

Raw WormLab exports should follow:

```text
data/raw/
└── <experiment_id>/
    └── wormlab_exports/
        └── <strain>/
            └── <condition>/
                ├── 0001/
                ├── 0002/
                └── ...
```

For example:

```text
data/raw/
└── 260311/
    └── wormlab_exports/
        └── N2/
            ├── 21/
            ├── 150/
            ├── 250/
            ├── 350/
            ├── 450/
            ├── 550/
            ├── 650/
            ├── 750/
            ├── 850/
            ├── 950/
            └── 1050/
```

## 5. Compile

Run:

```text
Notebooks/01_compile_genotype_data.ipynb
```

from top to bottom.

This generates the processed data and consolidated metadata required by Notebook 02.

## 6. Analyze

Run:

```text
Notebooks/02_celegans_volume_analysis.ipynb
```

from top to bottom.

Analysis outputs are written under:

```text
data/results/
```

## 7. Archive

When an analysis state should be preserved, run:

```text
Notebooks/03_archiving.ipynb
```

---

# Reproducibility Principles

This repository follows several analysis principles:

- raw measurements are kept separate from derived measurements
- metadata define experimental identity, ordering, and control status
- analysis parameters and thresholds are explicit
- individual worms are processed before group averaging
- missing or invalid kinetic measurements are not replaced with artificial values
- QC exclusions are documented rather than silently applied
- sensitivity analyses evaluate important methodological choices
- video identity is retained for potential cluster-aware analysis
- plotting and Prism exports derive from the same source data
- validated calculations are preserved when additional analyses are introduced

---

# Requirements

The workflow is Python-based and uses:

- pandas
- NumPy
- Matplotlib
- Seaborn
- SciPy
- PyYAML
- Jupyter
- IPython kernel support

See:

```text
requirements.txt
```

for the current environment dependencies.

---

# Citation

If this workflow or code is used in research, cite the associated publication or repository as appropriate.