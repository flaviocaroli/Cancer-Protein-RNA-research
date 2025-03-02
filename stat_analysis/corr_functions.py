import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def compute_correlation_for_protein(transcriptomics, proteomics, protein_name, p_threshold=0.05):
    """
    Computes Pearson correlations between the given protein and all genes.
    Returns a pd.Series where index = genes, value = correlation (or 0 if p >= p_threshold).
    
    transcriptomics: DataFrame (samples x genes)
    proteomics: DataFrame (samples x proteins)
    protein_name: str
    p_threshold: float, significance threshold for p-value
    """
    if protein_name not in proteomics.columns:
        raise ValueError(f"Protein '{protein_name}' not found in proteomics.")

    # Prepare a Series to store correlation values for each gene
    correlations = pd.Series(index=transcriptomics.columns, dtype=float)

    # Iterate over all genes
    for gene in transcriptomics.columns:
        # Overlapping samples that have non-null values for both the protein and this gene
        common_samples = proteomics[protein_name].dropna().index.intersection(
            transcriptomics[gene].dropna().index
        )
        # Pearson correlation requires at least 2 samples
        if len(common_samples) > 1:
            corr_val, p_val = pearsonr(
                proteomics.loc[common_samples, protein_name],
                transcriptomics.loc[common_samples, gene]
            )
            # Keep correlation if p < p_threshold; otherwise 0
            correlations[gene] = corr_val if p_val < p_threshold else 0.0
        else:
            correlations[gene] = 0.0

    return correlations
