import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = Path(r"e:\Todd\OSR\251128\OSR compiled data all genotypes")
genotypes = ['WT', 'dr170', 'dr180', 'dr170 dr180', 'osm6p811']
FIT_THRESHOLD = 0.9
V0_TIME_WINDOW = 5.0
T_LOW, T_HIGH = 500.0, 600.0

out_rows = []
per_track_rows = []
group_values = {}

for g in genotypes:
    area_fp = DATA_DIR / f"{g} compiled_area.csv"
    length_fp = DATA_DIR / f"{g} compiled_length.csv"
    fit_fp = DATA_DIR / f"{g} compiled_fit.csv"
    if not (area_fp.exists() and length_fp.exists() and fit_fp.exists()):
        print(f"Missing files for genotype {g}:\n  {area_fp}\n  {length_fp}\n  {fit_fp}")
        group_values[g] = np.array([])
        continue

    area = pd.read_csv(area_fp)
    length = pd.read_csv(length_fp)
    fit = pd.read_csv(fit_fp)

    track_cols = [c for c in area.columns if c not in ['Frame', 'Time']]
    # calculate volume
    vol = area.copy()
    for col in track_cols:
        # avoid division by zero
        vol[col] = (np.pi * (area[col].astype(float) ** 2)) / (4.0 * length[col].astype(float))

    # apply fit filter
    for col in track_cols:
        if col in fit.columns:
            mask = fit[col] < FIT_THRESHOLD
            vol.loc[mask, col] = np.nan

    # normalize to V0
    norm = vol.copy()
    v0_mask = vol['Time'] <= V0_TIME_WINDOW
    for col in track_cols:
        try:
            v0 = vol.loc[v0_mask, col].mean()
        except Exception:
            v0 = np.nan
        if pd.notna(v0) and v0 > 0:
            norm[col] = (vol[col] / v0) * 100.0
        else:
            norm[col] = np.nan

    # compute per-track mean in 500-600s window
    window_mask = (norm['Time'] >= T_LOW) & (norm['Time'] <= T_HIGH)
    per_track_means = norm.loc[window_mask, track_cols].mean(axis=0, skipna=True)
    per_track_means = per_track_means.dropna()

    group_values[g] = per_track_means.values

    # store per-track rows
    for track, val in per_track_means.items():
        per_track_rows.append({'Genotype': g, 'Track': track, 'Mean500_600': val})

# Save per-track means
per_track_df = pd.DataFrame(per_track_rows)
per_track_df.to_csv('per_track_means.csv', index=False)
print('Saved per-track means to per_track_means.csv')

# Prepare groups for testing (only groups with at least 2 tracks)
groups = {g: vals for g, vals in group_values.items() if len(vals) >= 2}
if len(groups) < 2:
    print('Not enough groups with data to compare. Exiting.')
    sys.exit(1)

# Normality tests (Shapiro) and report
normalities = {}
for g, vals in groups.items():
    # Shapiro requires 3<=n<=5000; if outside range, we'll skip or use normal approximation
    if len(vals) >= 3 and len(vals) <= 5000:
        stat, p = stats.shapiro(vals)
        normalities[g] = (stat, p)
    else:
        normalities[g] = (np.nan, np.nan)

print('\nShapiro-Wilk normality tests (statistic, p-value):')
for g, (stat, p) in normalities.items():
    print(f"  {g}: {stat:.4f}, p={p:.4g}")

# Levene test for equal variances
levene_stat, levene_p = stats.levene(*groups.values(), center='median')
print(f"\nLevene test for equal variances: stat={levene_stat:.4f}, p={levene_p:.4g}")

# Decide test
all_normal = all((np.isnan(v[1]) or v[1] > 0.05) for v in normalities.values())
homogeneous = levene_p > 0.05

# Try to import advanced posthoc packages
have_statsmodels = True
have_posthocs = True
try:
    import statsmodels.api as sm
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from statsmodels.stats.multitest import multipletests
except Exception:
    have_statsmodels = False

try:
    import scikit_posthocs as sp
except Exception:
    have_posthocs = False

results = []

if all_normal and homogeneous:
    print('\nAll groups approx normal and variances homogeneous -> running one-way ANOVA and Tukey HSD')
    fstat, fp = stats.f_oneway(*groups.values())
    print(f"ANOVA F={fstat:.4f}, p={fp:.4g}")
    results.append({'test': 'ANOVA', 'F': fstat, 'p': fp})

    if have_statsmodels:
        # assemble data for Tukey
        data = per_track_df.copy()
        data = data[data['Genotype'].isin(groups.keys())]
        tuk = pairwise_tukeyhsd(endog=data['Mean500_600'], groups=data['Genotype'], alpha=0.05)
        print('\nTukey HSD results:')
        print(tuk)
        tuk_df = pd.DataFrame(data=tuk._results_table.data[1:], columns=tuk._results_table.data[0])
        tuk_df.to_csv('tukey_hsd_results.csv', index=False)
        print('Saved Tukey HSD results to tukey_hsd_results.csv')
    else:
        print('statsmodels not available; cannot run Tukey HSD. Consider installing statsmodels.')

elif not all_normal:
    print('\nNot all groups pass normality -> using Kruskal-Wallis (non-parametric)')
    hstat, hp = stats.kruskal(*groups.values())
    print(f"Kruskal-Wallis H={hstat:.4f}, p={hp:.4g}")
    results.append({'test': 'Kruskal-Wallis', 'H': hstat, 'p': hp})

    if have_posthocs:
        data = per_track_df.copy()
        data = data[data['Genotype'].isin(groups.keys())]
        dunn = sp.posthoc_dunn(data, val_col='Mean500_600', group_col='Genotype', p_adjust='holm')
        dunn.to_csv('dunn_posthoc_results.csv')
        print('Saved Dunn posthoc results to dunn_posthoc_results.csv')
    else:
        print('scikit-posthocs not available; cannot run Dunn post-hoc. Consider installing scikit-posthocs.')

else:
    # normal but heteroscedastic
    print('\nData approx normal but variances unequal -> using pairwise Welch t-tests with p-value correction')
    gen_names = list(groups.keys())
    pairwise = []
    pvals = []
    comparisons = []
    for i in range(len(gen_names)):
        for j in range(i+1, len(gen_names)):
            a = groups[gen_names[i]]
            b = groups[gen_names[j]]
            tstat, pval = stats.ttest_ind(a, b, equal_var=False, nan_policy='omit')
            comparisons.append((gen_names[i], gen_names[j]))
            pvals.append(pval)
            pairwise.append({'group1': gen_names[i], 'group2': gen_names[j], 't': tstat, 'p': pval})
    # adjust p-values
    if have_statsmodels:
        reject, pvals_corrected, _, _ = multipletests(pvals, method='holm')
        for k, pw in enumerate(pairwise):
            pw['p_adj'] = pvals_corrected[k]
        pd.DataFrame(pairwise).to_csv('pairwise_welch_results.csv', index=False)
        print('Saved pairwise Welch test results to pairwise_welch_results.csv')
    else:
        print('statsmodels not available; saving unadjusted pairwise p-values')
        pd.DataFrame(pairwise).to_csv('pairwise_welch_results_unadj.csv', index=False)

# Save results summary
res_df = pd.DataFrame(results)
res_df.to_csv('stats_results_summary.csv', index=False)
print('\nSaved summary to stats_results_summary.csv')
print('Done')
