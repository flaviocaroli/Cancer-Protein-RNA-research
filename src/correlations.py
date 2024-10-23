# Imports

import pandas as pd
import notebooks.standardised_pipeline_utils as standardised_pipeline_utils

# Processing and matching

def initial_processing(transcriptomics: pd.DataFrame, proteomics: pd.DataFrame):

    transcriptomics = None
    proteomics = None
    sample_info = None

    # transcriptomics Updahya processing 
    transcriptomics = transcriptomics.set_index('Unnamed: 0')
    transcriptomics = transcriptomics.rename(index = dict(zip(sample_info['DepMap_ID'],
                                                          sample_info['CCLE_Name'])),
                                        columns = lambda x : str(x).split(' ')[0])
    
    transcriptomics = transcriptomics.transpose()
    assert len(transcriptomics.columns[transcriptomics.columns.duplicated()]) == 0, "columns contain duplicates"
    # drop rows (genes) with more than 20% NaNs (only removed 2 genes)
    # rows with same index are collapsed together into one row with avg values
    # drop rows that are indexed by datetime object
    # drop rows containing ":"
    transcriptomics_processed = standardised_pipeline_utils.process(transcriptomics)
    transcriptomics_processed = transcriptomics_processed.transpose()
    
    ### proteomics Updahya processing 
    proteomics.set_index('Gene_Symbol', inplace=True)
    proteomics = proteomics.loc[:, proteomics.columns.str.contains('_TenPx')]
    # Checking for cell lines repeated in >1 Ten-plexes
    proteomics.filter(regex='SW948_LARGE_INTESTINE|CAL120_BREAST|HCT15_LARGE_INTESTINE').columns   
    # Eliminating the cell lines that do not correlate well with transcriptomics data as mentioned in the paper 
    proteomics.drop(columns=['SW948_LARGE_INTESTINE_TenPx11', 'CAL120_BREAST_TenPx02', 'HCT15_LARGE_INTESTINE_TenPx30'], 
                        inplace=True)
    proteomics = proteomics.rename(columns = lambda x : str(x).split('_TenPx')[0]) 
    assert len(proteomics.columns[proteomics.columns.duplicated()]) == 0, "columns contain duplicates"  
    proteomics_processed = standardised_pipeline_utils.process(proteomics)
    proteomics_processed = proteomics_processed.transpose()

    ### Matching cell lines

    transcriptomics = transcriptomics_processed
    proteomics = proteomics_processed
    matched_cell_lines = [line for line in proteomics.index if line in transcriptomics.index]
    transcriptomics_matched = transcriptomics.reindex(matched_cell_lines)
    proteomics_matched = proteomics.reindex(matched_cell_lines)
    # only keeping the 369 cell lines present both in transcriptomics and proteomics data
    return transcriptomics_matched, proteomics_matched

def correlations(transcriptomics_matched, proteomics_matched):
    transcriptomics_corr, proteomics_corr = standardised_pipeline_utils.match_proteins_samples(transcriptomics_matched.T,
                                                                                           proteomics_matched.T)
    #Pearson
    pearson_correlations = transcriptomics_corr.corrwith(proteomics_corr,
                                                      axis=1,
                                                      method='pearson')
    #Spearman
    spearman_correlations = transcriptomics_corr.corrwith(proteomics_corr,
                                                      axis=1,
                                                      method='spearman')
    
    #Kendall
    kendall_correlations = transcriptomics_corr.corrwith(proteomics_corr,
                                                      axis=1,
                                                      method='kendall')
    correlations = pd.DataFrame({'pearson': pearson_correlations,
                             'spearman': spearman_correlations,
                             'kendall': kendall_correlations})
    
    return correlations