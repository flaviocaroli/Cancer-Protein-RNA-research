# Imports 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import notebooks.standardised_pipeline_utils as standardised_pipeline_utils

def second_processing(transcriptomics, proteomics, transcriptomics_matched):
    # removed row with nan's. no more nan's in RNA data.
    transcriptomics.drop(index="TT_OESOPHAGUS", inplace=True)

    proteomics = proteomics.dropna(axis=1)
    ### Re-match no-nan data    
    indices = list(transcriptomics_matched.index)

    transcriptomics = transcriptomics.reindex(indices)
    proteomics = proteomics.reindex(indices) 

    # Removing unexpressed genes from transcriptomics
    # Removing genes with summed expression below 0.5
    transcriptomics = transcriptomics.loc[:, transcriptomics.sum() >= 0.5]   
    
    #LogLog data
    # Applying another log-transformation
    transcriptomics = np.log2(1000 * transcriptomics + 1)

    return transcriptomics, proteomics