# C. elegans Volume Analysis

Analysis of volume changes over time in *C. elegans* across multiple genotypes using Wormlab tracking data and reproducible notebooks.

## Repository Layout

```
celegans-volume-analysis/
├── Notebooks/
│   ├── 01_compile_genotype_data.ipynb      # Interactive data compilation
│   ├── 02_celegans_volume_analysis.ipynb   # Primary analysis workflow
│   └── 03_archiving.ipynb                  # Archive complete analysis
├── scripts/
│   ├── compile_genotype_data.py            # Builds compiled CSVs per metric
│   └── run_volume_analysis.py              # Convenience runner (stub)
├── data/
│   ├── raw/                                # Raw tracking data by genotype/video
│   │   └── <genotype>/
│   │       ├── config.yaml                 # Metadata (strain, NaCl, frame rate, etc.)
│   │       ├── 0001/                       # Video folder
│   │       │   ├── area.csv
│   │       │   ├── length.csv
│   │       │   ├── position.csv
│   │       │   └── fit.csv
│   │       └── 0002/                       # Another video
│   ├── processed/                          # Generated compiled CSVs
│   └── archive/                            # Timestamped analysis snapshots
├── sample.yaml                             # Template config.yaml file
├── requirements.txt
└── README.md
```

## Notebook 01 – Compile Genotype Data

Use `Notebooks/01_compile_genotype_data.ipynb` to build compiled CSVs from raw tracking data. The notebook:

1. Loads metadata from `config.yaml` files in each genotype folder (containing strain, genotype, NaCl concentration, frame rate, pixel resolution, and other experimental parameters).
2. Discovers video folders (0001, 0002, etc.) containing area.csv, length.csv, position.csv, and fit.csv files.
3. Aligns `Frame`/`Time` columns across metrics and compiles tracks side-by-side.
4. Writes `<genotype>_compiled_<metric>.csv` files to `data/processed/` for each genotype.
5. Creates summary CSVs: `compilation_summary.csv` (file paths) and `genotype_metadata.csv` (all config.yaml metadata).

## Notebook 02 – Volume Analysis Workflow

Cell order in `Notebooks/02_celegans_volume_analysis.ipynb` reflects the updated analysis strategy:

1. **N2 preprocessing** – cleans column headers, matches area/length/fit tracks, and writes `N2_compiled_volume.csv`.
2. **N2 visualization** – plots raw per-track volumes.
3. **N2 normalization** – normalizes to the first 5 s, applies a 28-frame rolling mean, and plots mean ± SEM.
4. **Multi-genotype pipeline** – iterates over every compiled genotype, applies shared filters, produces both smoothed raw normalized curves (56-frame window) and per-track exponential-fit curves, then plots genotype overlays for each representation.
5. **Fit summary** – aggregates the per-track exponential parameters (A, k, half-life) into genotype-level tables for downstream comparisons.
6. **Per-track decay statistics** – plots per-track exponential k values as bar charts with individual points, reports ANOVA plus Holm-adjusted pairwise Welch tests, and annotates group sizes to highlight significant genotype differences.

## Notebook 03 – Archiving

Use `Notebooks/03_archiving.ipynb` to create timestamped archives of complete analyses:

1. Packages all raw data (including config.yaml files), processed data, notebooks, scripts, and configuration files.
2. Creates a timestamped archive in `data/archive/celegans_analysis_YYYYMMDD_HHMMSS/`.
3. Generates a `MANIFEST.json` with metadata about all archived files.
4. Provides a self-contained, reproducible analysis package for sharing or long-term storage.

## Volume, Filtering, and Normalization

- **Volume formula:** $V = \frac{\pi \cdot \text{Area}^2}{4 \cdot \text{Length}}$
- **Track filters:**
    - Fit score ≥ 0.9
    - Initial volume ≥ 1.5 µm³
    - Duration ≥ 400 s of valid measurements
- **Normalization:** Each track is scaled to the mean volume within the first 5 s (V0) and reported as %V0.
- **Smoothing:**
    - Single-genotype plots use a 28-frame rolling mean.
    - Multiyour data:
    - Create genotype folders in `data/raw/<genotype>/`
    - Add a `config.yaml` file in each genotype folder (see `sample.yaml` for template)
    - Place video folders (0001, 0002, etc.) with tracking CSVs inside each genotype folder
3. Compile the data:
    - Open `Notebooks/01_compile_genotype_data.ipynb` and run all cells to generate compiled CSVs in `data/processed/`
4. Analyze:
    - Laprocessed/<genotype>_compiled_<metric>.csv` – compiled tracking data per genotype and metric
- `data/processed/compilation_summary.csv` – index of all compiled files
- `data/processed/genotype_metadata.csv` – consolidated metadata from all config.yaml files
- `data/processed/*_compiled_volume.csv` – per-genotype volume time series after filtering
- `data/archive/celegans_analysis_YYYYMMDD_HHMMSS/` – timestamped complete analysis archives
- Notebook-rendered plots:
    - Smoothed raw normalized mean ± SEM across genotypes
    - Exponential-fit mean ± SEM across genotypes
    - Per-track exponential decay rates with individual track points, group means, and ±SE bars
- Fit parameter tables showing `A`, `k`, and half-life per genotype
    ```
2. Prepare compiled inputs (pick one of the options):
    - **Interactive:** open `Notebooks/01_compile_genotype_data.ipynb`, set the raw-data root if needed, and run all cells to regenerate the compiled area/length/fit CSVs.
    - **Batch script:** run `python scripts/compile_genotype_data.py` for the same result without a notebook.
3. Launch `Notebooks/02_celegans_volume_analysis.ipynb` in Jupyter (VS Code, JupyterLab, etc.) and run the cells top to bottom. The notebook writes per-genotype volume CSVs, plots both curve types, and saves summary tables in `celegans-volume-analysis/`.

## Key Outputs

- `data/compiled/*_compiled_volume.csv` – per-genotype volume time series after filtering.
- `summary_*.csv` – aggregated statistics per genotype and across genotypes.
- Notebook-rendered plots:
    - Smoothed raw normalized mean ± SEM across genotypes.
    - Exponential-fit mean ± SEM across genotypes.
    - Per-track exponential decay rates with individual track points, group means, and ±SE bars.
- Fit parameter tables showing `A`, `k`, and half-life per genotype.

## Notes and Assumptions, pyyaml (see `requirements.txt`)

- Wormlab exports must share aligned `Frame`/`Time` series per genotype; misaligned tracks are logged and aligned by frame when possible.
- All volumes are expressed in µm³; normalized plots show %V0.
- Update `FIT_THRESHOLD`, `INITIAL_VOLUME_MIN`, `MIN_DURATION`, `V0_WINDOW`, `SMOOTH_WINDOW`, and `FIT_WINDOW` within the notebook to explore alternative criteria.

## Requirements

- Python 3.9+
- pandas, numpy, matplotlib, seaborn, scipy, jupyter (see `requirements.txt`).

## Citation

If you use this workflow in your research, please cite appropriately.
