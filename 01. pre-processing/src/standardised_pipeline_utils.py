import datetime 
import numpy as np 
import pandas as pd 

get_common_proteins = lambda df1, df2: df1.index[df1.index.isin(df2.index)]
get_common_samples = lambda df1, df2: np.intersect1d(df1.columns, df2.columns)

def match_proteins_samples(df1, df2):
    common_proteins = get_common_proteins(df1, df2)
    common_samples = get_common_samples(df1, df2)
    print("Number of common proteins: ", len(common_proteins))
    print("Number of common samples: ", len(common_samples))
    
    return df1.reindex(common_proteins)[common_samples], df2.reindex(common_proteins)[common_samples]

def correlate_genewise(df1, df2, cname, method='spearman'):
    df1_subset, df2_subset = df1.copy(), df2.copy()
    correlation = df1_subset.corrwith(df2_subset, axis=1, method=method)
    print("Median", method.title(), "Correlation: ", round(correlation.median(), 4))    
    return correlation.to_frame(cname)

def dropna(dataframe, non_null_threshold=0.8, replace_zero=False):
    if(replace_zero):
        subset = dataframe.copy(deep=True).replace(0, np.nan)
    else: 
        subset = dataframe.copy(deep=True)
    non_null_columns = len(dataframe.columns)*non_null_threshold
    subset.dropna(thresh=non_null_columns, inplace=True) 
    return subset
   
#eliminate rows with >20% of null values and compute mean value for the protein isoforms
def process(dataframe):
    # If the data does not contain NA values (this is especially seen in the older studies transcriptomic data), 
    # replace the zero with NA values and then drop the rows with >20% NA values
    replace_zero = True if (dataframe.isnull().sum().sum() == 0) else False
    dataframe_processed = dropna(dataframe, replace_zero=replace_zero)
    dataframe_processed = dataframe_processed.groupby(dataframe_processed.index).mean()
    dataframe_processed.drop(index=[index for index in dataframe_processed.index if type(index) is datetime.datetime], inplace=True)
    dataframe_processed.drop(index=[index for index in dataframe_processed.index if ':' in index], inplace=True)
    print("Dimensions: ", dataframe_processed.shape)
    return dataframe_processed