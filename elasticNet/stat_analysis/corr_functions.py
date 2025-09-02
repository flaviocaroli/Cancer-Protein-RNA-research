import pandas as pd 
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr


def compute_correlation_for_protein(transcriptomics, proteomics, protein_name, p_threshold=0.05, bootstrap=True, n_bootstraps=100):
    """
    Computes Pearson correlations between the given protein and all genes.
    
    For each gene:
      - If bootstrap=False, uses scipy.stats.pearsonr and returns the correlation if p-value < p_threshold, otherwise 0.
      - If bootstrap=True, performs bootstrapping to compute a 95% confidence interval and the median correlation.
        The median bootstrap correlation is returned if the confidence interval excludes 0; otherwise 0 is returned.
    
    Parameters:
      transcriptomics: DataFrame (samples x genes)
      proteomics: DataFrame (samples x proteins)
      protein_name: str, name of the protein in proteomics DataFrame
      p_threshold: float, significance threshold for the p-value (used when bootstrap=False)
      bootstrap: bool, if True, use bootstrapping to assess significance instead of the parametric p-value
      n_bootstraps: int, number of bootstrap iterations if bootstrap=True
    
    Returns:
      pd.Series where index = genes, and value = correlation (or 0 if not significant)
    """
    if protein_name not in proteomics.columns:
        raise ValueError(f"Protein '{protein_name}' not found in proteomics.")
    
    # Prepare a Series to store correlation values for each gene
    correlations = pd.Series(index=transcriptomics.columns, dtype=float)

    # Iterate over all genes
    for gene in transcriptomics.columns:
        # Identify samples with non-null values for both the protein and the gene
        common_samples = proteomics[protein_name].dropna().index.intersection(
            transcriptomics[gene].dropna().index
        )
        # Pearson correlation requires at least 2 samples
        if len(common_samples) > 1:
            data_prot = proteomics.loc[common_samples, protein_name].values
            data_gene = transcriptomics.loc[common_samples, gene].values

            if not bootstrap:
                # Compute using the standard Pearson correlation test
                corr_val, p_val = pearsonr(data_prot, data_gene)
                correlations[gene] = corr_val if p_val < p_threshold else 0.0
            else:
                # Compute the observed correlation (for reference)
                observed_corr, _ = pearsonr(data_prot, data_gene)
                boot_corrs = []
                # Bootstrap resampling: sample with replacement and compute correlations
                for _ in range(n_bootstraps):
                    indices = np.random.choice(len(common_samples), size=len(common_samples), replace=True)
                    boot_corr, _ = pearsonr(data_prot[indices], data_gene[indices])
                    boot_corrs.append(boot_corr)
                # Compute 95% confidence interval from bootstrap distribution
                ci_lower, ci_upper = np.percentile(boot_corrs, [2.5, 97.5])
                # Compute the median of the bootstrap correlations
                median_boot = np.median(boot_corrs)
                # Keep the median correlation only if the confidence interval excludes 0
                if ci_lower > 0 or ci_upper < 0:
                    correlations[gene] = median_boot
                else:
                    correlations[gene] = 0.0
        else:
            correlations[gene] = 0.0

    return correlations
