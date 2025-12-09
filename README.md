# C. elegans Volume Analysis

Analysis of volume changes over time in *C. elegans* across different genotypes using Wormlab tracking data.

## Overview

This repository contains a Jupyter notebook for analyzing volume changes in *C. elegans* worms over time. The analysis compares four genotypes:
- **WT** (wildtype)
- **dr170**
- **dr180**
- **dr170 dr180** (double mutant)

## Volume Calculation

Volume is calculated using the formula:

```
V ≈ (π × Area²) / (4 × Length)
```

Where:
- **Area**: Projected area from Wormlab (μm²)
- **Length**: Body length from Wormlab (μm)
- **Volume**: Calculated volume (μm³)

## Features

- Load area, length, and fit data from Wormlab CSV files
- Calculate volume for each tracked worm
- Filter tracks by fit quality (default threshold: 0.9)
- Normalize volume to baseline (V0 = average volume over first 5 seconds)
- Apply 5-second rolling average smoothing
- Generate plots with mean ± SEM for each genotype
- Export summary statistics to CSV files

## Data Structure

Expected directory structure:
```
├── WT compiled_data/
│   ├── WT compiled_area.csv
│   ├── WT compiled_length.csv
│   └── WT compiled_fit.csv
├── dr170 compiled_data (1)/
│   ├── dr170 compiled_area.csv
│   ├── dr170 compiled_length.csv
│   └── dr170 compiled_fit.csv
├── dr180 compiled_data/
│   ├── dr180 compiled_area.csv
│   ├── dr180 compiled_length.csv
│   └── dr180 compiled_fit.csv
└── dr170 dr180 compiled_data (1)/
    ├── dr170 dr180 compiled_area.csv
    ├── dr170 dr180 compiled_length.csv
    └── dr170 dr180 compiled_fit.csv
```

## Usage

### Installation

1. Clone this repository
2. Install required packages:
```bash
pip install pandas numpy matplotlib seaborn scipy
```

Or run the first cell in the notebook:
```python
%pip install pandas numpy matplotlib seaborn scipy
```

### Running the Analysis

Open `celegans_volume_analysis.ipynb` in Jupyter and run cells sequentially:

1. **Step 0**: Install packages (if needed)
2. **Step 1**: Import libraries
3. **Step 2**: Define parameters
4. **Step 3-8**: Process and plot WT data (example walkthrough)
5. **Step 9**: Process all remaining genotypes
6. **Step 10**: Generate combined plot with all genotypes
7. **Step 11-13**: Summary statistics and data export

## Parameters

Customize analysis parameters in Step 2:

```python
FIT_THRESHOLD = 0.9          # Minimum fit quality (0-1)
V0_TIME_WINDOW = 5.0         # Baseline window (seconds)
SMOOTHING_WINDOW = 5.0       # Smoothing window (seconds)
```

## Color Scheme

- **WT**: Black line with grey SEM envelope
- **dr170**: Red line with pink SEM envelope
- **dr180**: Blue line with light blue SEM envelope
- **dr170 dr180**: Orange line with light orange SEM envelope

## Output Files

- `volume_analysis_all_genotypes.png` - Combined plot
- `summary_all_genotypes.csv` - Combined summary data
- `summary_WT.csv`, `summary_dr170.csv`, etc. - Individual genotype summaries

## Requirements

- Python 3.7+
- pandas
- numpy
- matplotlib
- seaborn
- scipy

## License

This analysis code is provided as-is for scientific research purposes.

## Citation

If you use this analysis in your research, please cite appropriately.
