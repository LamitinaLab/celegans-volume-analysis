#!/usr/bin/env python3
"""Add a new 21 mM plot cell to the notebook."""

import json
from pathlib import Path

nb_path = Path('Notebooks/02_celegans_volume_analysis.ipynb')

with open(nb_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

# Create the new cell
new_cell = {
    'cell_type': 'code',
    'metadata': {'language': 'python'},
    'source': [
        '# Individual %deltaV/V0 tracks for 21 mM - all tracks overlaid\n',
        'import matplotlib.pyplot as plt\n',
        'import numpy as np\n',
        'import pandas as pd\n',
        'import seaborn as sns\n',
        '\n',
        "TARGET_CONDITION = '21'\n",
        "SMOOTH_W = int(globals().get('SMOOTH_WINDOW', 5))\n",
        '\n',
        '# Get data from volume_long\n',
        "if 'volume_long' in globals() and isinstance(volume_long, pd.DataFrame):\n",
        "    df = volume_long[volume_long['Condition_mM'].astype(str) == TARGET_CONDITION].copy()\n",
        'else:\n',
        "    raise ValueError('volume_long not found. Run volume analysis pipeline first.')\n",
        '\n',
        'if df.empty:\n',
        "    raise ValueError(f'No data found for {TARGET_CONDITION} mM.')\n",
        '\n',
        '# Prepare data\n',
        "df = df[['Condition_mM', 'Track', 'Time', 'Volume_norm_pct']].copy()\n",
        "df['Time'] = pd.to_numeric(df['Time'], errors='coerce')\n",
        "df['Volume_norm_pct'] = pd.to_numeric(df['Volume_norm_pct'], errors='coerce')\n",
        "df = df.dropna(subset=['Time', 'Volume_norm_pct'])\n",
        '\n',
        '# Calculate %deltaV/V0\n',
        "df['DeltaV_over_V0_pct'] = df['Volume_norm_pct'] - 100.0\n",
        '\n',
        '# Smooth each track\n',
        "df = df.sort_values(['Track', 'Time']).copy()\n",
        "df['Delta_smooth'] = (\n",
        "    df.groupby('Track', sort=False)['DeltaV_over_V0_pct']\n",
        "      .transform(lambda s: s.rolling(window=SMOOTH_W, min_periods=1, center=True).mean())\n",
        ')\n',
        '\n',
        '# Create overlay plot\n',
        "sns.set_style('whitegrid')\n",
        'fig, ax = plt.subplots(figsize=(12, 6))\n',
        '\n',
        "tracks = sorted(df['Track'].dropna().unique().tolist())\n",
        'colors = plt.cm.tab20(np.linspace(0, 1, len(tracks)))\n',
        '\n',
        'for i, tr in enumerate(tracks):\n',
        "    dtr = df[df['Track'] == tr].sort_values('Time')\n",
        '    ax.plot(\n',
        "        dtr['Time'], dtr['Delta_smooth'],\n",
        '        linewidth=1.5, alpha=0.7,\n',
        '        label=tr, color=colors[i]\n',
        '    )\n',
        '\n',
        "ax.axhline(0.0, linestyle='--', linewidth=1.0, color='black', alpha=0.5)\n",
        "ax.set_xlabel('Time (s)', fontsize=12)\n",
        "ax.set_ylabel('%deltaV/V0', fontsize=12)\n",
        "ax.set_title(f'{TARGET_CONDITION} mM: Individual %deltaV/V0 Trajectories (All Tracks Overlaid)', fontsize=14)\n",
        'ax.grid(alpha=0.3)\n',
        "ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)\n",
        '\n',
        'plt.tight_layout()\n',
        '\n',
        "if 'RESULTS_DIR' in globals():\n",
        "    out_png = RESULTS_DIR / f'individual_deltaV_over_V0_{TARGET_CONDITION}mM_overlay.png'\n",
        "    out_pdf = RESULTS_DIR / f'individual_deltaV_over_V0_{TARGET_CONDITION}mM_overlay.pdf'\n",
        "    fig.savefig(out_png, dpi=300, bbox_inches='tight')\n",
        "    fig.savefig(out_pdf, bbox_inches='tight')\n",
        "    print(f'Saved: {out_png}')\n",
        "    print(f'Saved: {out_pdf}')\n",
        '\n',
        'plt.show()'
    ]
}

# Add the new cell
nb['cells'].append(new_cell)

# Save the notebook
with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print('New cell added successfully!')
