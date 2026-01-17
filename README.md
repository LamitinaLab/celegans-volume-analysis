# C. elegans Volume Analysis

Analysis of volume changes over time in *C. elegans* across multiple genotypes using Wormlab tracking data and reproducible notebooks.

## Repository Layout

```
celegans-volume-analysis/
├── Notebooks/
│   ├── 01_compile_genotype_data.ipynb      # Interactive data compilation
│   └── 02_celegans_volume_analysis.ipynb   # Primary analysis workflow
├── scripts/
│   ├── compile_genotype_data.py            # Builds compiled CSVs per metric
│   └── run_volume_analysis.py              # Convenience runner (stub)
├── data/
│   ├── raw/                                # Wormlab exports (per genotype)
│   └── compiled/                           # Generated *_compiled_*.csv files
├── requirements.txt
└── README.md
```

## Notebook 01 – Compile Genotype Data

Use `Notebooks/01_compile_genotype_data.ipynb` when you need to rebuild the compiled CSVs from freshly exported Wormlab files. The notebook walks through:

1. Locating raw CSV folders under `data/raw/` and normalizing their headers.
2. Aligning `Frame`/`Time` columns across area/length/fit metrics.
3. Writing `*_compiled_area.csv`, `*_compiled_length.csv`, and `*_compiled_fit.csv` into `data/compiled/` for every genotype.

This notebook mirrors the logic in `scripts/compile_genotype_data.py`, but keeps the steps visible in Jupyter with intermediate diagnostics.

## Notebook 02 – Volume Analysis Workflow

Cell order in `Notebooks/02_celegans_volume_analysis.ipynb` reflects the updated analysis strategy:

1. **N2 preprocessing** – cleans column headers, matches area/length/fit tracks, and writes `N2_compiled_volume.csv`.
2. **N2 visualization** – plots raw per-track volumes.
3. **N2 normalization** – normalizes to the first 5 s, applies a 28-frame rolling mean, and plots mean ± SEM.
4. **Multi-genotype pipeline** – iterates over every compiled genotype, applies shared filters, produces both smoothed raw normalized curves (56-frame window) and per-track exponential-fit curves, then plots genotype overlays for each representation.
5. **Fit summary** – aggregates the per-track exponential parameters (A, k, half-life) into genotype-level tables for downstream comparisons.
6. **Per-track decay statistics** – plots per-track exponential k values as bar charts with individual points, reports ANOVA plus Holm-adjusted pairwise Welch tests, and annotates group sizes to highlight significant genotype differences.

## Volume, Filtering, and Normalization

- **Volume formula:** $V = \frac{\pi \cdot \text{Area}^2}{4 \cdot \text{Length}}$
- **Track filters:**
    - Fit score ≥ 0.9
    - Initial volume ≥ 1.5 µm³
    - Duration ≥ 400 s of valid measurements
- **Normalization:** Each track is scaled to the mean volume within the first 5 s (V0) and reported as %V0.
- **Smoothing:**
    - Single-genotype plots use a 28-frame rolling mean.
    - Multi-genotype raw overlays use a 56-frame rolling mean to reduce noise before averaging.
- **Exponential fitting:** Each normalized track is fit to $V = A e^{k t}$ within a configurable window (default 0–120 s). The mean ± SEM of these fitted curves is plotted separately from the raw smoothed curves.

## Running the Analysis

1. Install dependencies (in a virtual environment if desired):
    ```bash
    pip install -r requirements.txt
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

## Notes and Assumptions

- Wormlab exports must share aligned `Frame`/`Time` series per genotype; misaligned tracks are logged and aligned by frame when possible.
- All volumes are expressed in µm³; normalized plots show %V0.
- Update `FIT_THRESHOLD`, `INITIAL_VOLUME_MIN`, `MIN_DURATION`, `V0_WINDOW`, `SMOOTH_WINDOW`, and `FIT_WINDOW` within the notebook to explore alternative criteria.

## Requirements

- Python 3.9+
- pandas, numpy, matplotlib, seaborn, scipy, jupyter (see `requirements.txt`).

## Citation

If you use this workflow in your research, please cite appropriately.
