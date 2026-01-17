from pathlib import Path
import pandas as pd
import re

DATA_ROOT = Path(r"e:\Todd\OSR\251128\OSR compiled data all genotypes\data\raw")
OUT_DIR = Path(r"e:\Todd\OSR\251128\OSR compiled data all genotypes\data\compiled")
OUT_DIR.mkdir(parents=True, exist_ok=True)

METRICS = ['area', 'length', 'fit']

def sanitize_prefix(s: str) -> str:
    s = re.sub(r"[^0-9A-Za-z_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip('_')
    return s


def collect_files_for_genotype(geno_dir: Path):
    files_by_metric = {m: [] for m in METRICS}
    for p in geno_dir.rglob('*.csv'):
        name = p.stem.lower()
        for m in METRICS:
            if m in name:
                files_by_metric[m].append(p)
                break
    return files_by_metric


def read_with_header(path: Path):
    try:
        df = pd.read_csv(path, skiprows=1)
    except Exception as e:
        print(f"ERROR reading {path}: {e}")
        raise
    return df


def compile_metric_side_by_side(genotype: str, files: list, out_dir: Path, metric: str):
    if not files:
        print(f"No files for {genotype} / {metric}")
        return None

    base_df = None
    out_df = None
    for f in sorted(files):
        df = read_with_header(f)
        if 'Frame' not in df.columns or 'Time' not in df.columns:
            print(f"Skipping {f} (missing Frame/Time)")
            continue

        # Rename track columns
        cols = [c for c in df.columns if c not in ['Frame', 'Time']]
        prefix = sanitize_prefix(f.stem)
        rename_map = {c: f"{prefix}_{c}" for c in cols}
        df = df.rename(columns=rename_map)

        if base_df is None:
            base_df = df[['Frame', 'Time']].copy()
            out_df = base_df.copy()

        else:
            # Validate Frame/Time alignment
            if not (df['Frame'].equals(base_df['Frame']) and df['Time'].equals(base_df['Time'])):
                print(f"Warning: Frame/Time mismatch for {f} (attempting index alignment)")
                # Attempt to align by Frame if possible
                df = df.set_index('Frame')
                out_df = out_df.set_index('Frame')

                # Prepare tracks and ensure unique column names before joining
                df_tracks = df.drop(columns=['Time']).copy()
                for col in list(df_tracks.columns):
                    new_col = col
                    i = 1
                    while new_col in out_df.columns:
                        new_col = f"{col}__dup{i}"
                        i += 1
                    if new_col != col:
                        df_tracks = df_tracks.rename(columns={col: new_col})

                out_df = out_df.join(df_tracks, how='outer')
                out_df = out_df.reset_index()
                if 'Time' not in out_df.columns and 'Time' in df.columns:
                    out_df['Time'] = df['Time'].values
                continue

        # Append track columns
        track_cols = [c for c in df.columns if c not in ['Frame', 'Time']]
        out_df = pd.concat([out_df.reset_index(drop=True), df[track_cols].reset_index(drop=True)], axis=1)

    if out_df is None:
        print(f"No valid data for {genotype}/{metric}")
        return None

    out_path = out_dir / f"{genotype}_compiled_{metric}.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Wrote {out_path} (shape={out_df.shape})")
    return out_path


def main():
    if not DATA_ROOT.exists():
        print(f"DATA_ROOT not found: {DATA_ROOT}")
        return

    for geno_dir in sorted(p for p in DATA_ROOT.iterdir() if p.is_dir()):
        genotype = geno_dir.name.replace(' ', '_')
        print(f"\nProcessing genotype: {genotype}")
        files_by_metric = collect_files_for_genotype(geno_dir)
        for metric in METRICS:
            compile_metric_side_by_side(genotype, files_by_metric.get(metric, []), OUT_DIR, metric)

    print('\nCompilation complete.')

if __name__ == '__main__':
    main()
