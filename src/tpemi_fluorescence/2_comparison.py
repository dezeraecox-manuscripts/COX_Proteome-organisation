import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import statsmodels.api as sa
import statsmodels.formula.api as sfa

from statsmodels import stats

from loguru import logger

logger.info('Import OK')

input_path = 'results/tpemi_fluorescence/initial_cleanup/compiled.csv'
output_folder = 'results/tpemi_fluorescence/comparison/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def apply_oneway_anova(df, xcol, group_col):
    ols_model = f'{xcol} ~ C({group_col})'
    lm = sfa.ols(ols_model, data=df).fit()
    anova = sa.stats.anova_lm(lm)

    # complete post-hoc test with bonferroni correction
    pairs = lm.t_test_pairwise(f"C({group_col})", method='bonferroni')

    return anova, pairs.result_frame

# Read in summarised data
compiled = pd.read_csv(input_path)
compiled.drop([col for col in compiled.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Test difference to 1 using one-samp ttest
ttests = {
    treatment: list(ttest_1samp(df['norm_fluorescence'], popmean=1))
    for treatment, df in compiled.groupby('treatment')
}
ttests = pd.DataFrame(ttests).T.reset_index()
ttests.columns = ['treatment', 'tstat', 'pval']
ttests['pval_adjusted'] = stats.multitest.multipletests(
    ttests['pval'], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]
ttests['sign'] = ['***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else np.nan)) for pval in ttests['pval']]
ttests.to_csv(f'{output_folder}ttest_result.csv')

# Test difference using ANOVA with post-hoc
anova, multi_comparison = apply_oneway_anova(
    compiled, xcol='norm_fluorescence', group_col='treatment')
anova.to_csv(f'{output_folder}anova.csv')
multi_comparison.to_csv(f'{output_folder}multi_comparison.csv')
