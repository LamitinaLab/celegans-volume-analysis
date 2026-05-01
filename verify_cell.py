#!/usr/bin/env python3
import json
from pathlib import Path

nb_path = Path('Notebooks/02_celegans_volume_analysis.ipynb')
with open(nb_path) as f:
    nb = json.load(f)

print(f'Total cells: {len(nb["cells"])}')
last_cell = nb['cells'][-1]
print(f'Last cell type: {last_cell["cell_type"]}')
print(f'First line: {last_cell["source"][0][:80] if last_cell["source"] else "empty"}')
