import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Load per-track means generated earlier
csv_fp = 'per_track_means.csv'
df = pd.read_csv(csv_fp)

# Genotypes
groups = sorted(df['Genotype'].unique())
if 'WT' not in groups:
    raise SystemExit('WT not found in per_track_means.csv')

wt_vals = df.loc[df['Genotype'] == 'WT', 'Mean500_600'].dropna().values
results = []
paired = []

for g in groups:
    if g == 'WT':
        continue
    vals = df.loc[df['Genotype'] == g, 'Mean500_600'].dropna().values
    n_wt = len(wt_vals)
    n_g = len(vals)
    if n_wt < 2 or n_g < 2:
        results.append({'Genotype': g, 'n_WT': n_wt, 'n_group': n_g, 'U': np.nan, 'p_uncorrected': np.nan, 'p_adj': np.nan, 'rank_biserial': np.nan, 'note': 'insufficient_n'})
        continue
    # Mann-Whitney U test (two-sided)
    try:
        u_stat, pval = stats.mannwhitneyu(wt_vals, vals, alternative='two-sided')
    except Exception as e:
        results.append({'Genotype': g, 'n_WT': n_wt, 'n_group': n_g, 'U': np.nan, 'p_uncorrected': np.nan, 'p_adj': np.nan, 'rank_biserial': np.nan, 'note': f'error:{e}'})
        continue
    # rank-biserial effect size: r = 1 - (2*U)/(n1*n2)
    r_rb = 1.0 - (2.0 * u_stat) / (n_wt * n_g)
    results.append({'Genotype': g, 'n_WT': n_wt, 'n_group': n_g, 'U': float(u_stat), 'p_uncorrected': float(pval), 'p_adj': np.nan, 'rank_biserial': float(r_rb), 'note': ''})
    paired.append(pval)

# Adjust p-values (Holm)
if len(paired) > 0:
    reject, p_adj, _, _ = multipletests(paired, method='holm')
    # write adjusted back
    k = 0
    for i, row in enumerate(results):
        if not np.isnan(row['p_uncorrected']):
            results[i]['p_adj'] = float(p_adj[k])
            results[i]['significant'] = bool(reject[k])
            k += 1
        else:
            results[i]['significant'] = False

res_df = pd.DataFrame(results)
res_df = res_df[['Genotype', 'n_WT', 'n_group', 'U', 'rank_biserial', 'p_uncorrected', 'p_adj', 'significant', 'note']]
res_df.to_csv('pairwise_WT_results.csv', index=False)
print('Saved pairwise_WT_results.csv')
print('\nSummary:')
print(res_df.to_string(index=False))
