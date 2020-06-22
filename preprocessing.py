import pandas as pd 
import numpy as np


def preprocess(df):
    """
    Does all preprocessing steps including log transform and quantile normalization.
    """
    log_transformed_data = log_transform(df)
    quantile_normalized_data = quantile_norm(log_transformed_data)
    return quantile_normalized_data


def log_transform(df):
    """
    Compute log transform to initially normalize data and decrease skewness.
    """
    return np.log2(df + 1)


def quantile_norm(df):
    """
    Compute quantile normalization. 
    First, get means of each row from matrix with sorted columns.
    For loop: 
        1. Create list of tuples (index, value) from original matrix.
        2. Sort each tuple by value while preserving original index.
        3. Replace value of each tuple with rank means (in order).
        4. Reorder tuples by index (to get original order).
        5. Replace original col with new quantile-normalized col.
    If the original matrix has repeated values, I preserved the original
    order of the values.
    """
    rank = mean_rank(df)
    for col in df:
        d = list(df[col].to_dict().items())
        d = sorted(d, key=lambda x:x[1])
        d = [(d[index][0],rank[index]) for index in range(len(d))]
        d = sorted(d, key=lambda x:x[0])
        df[col] = [val for index,val in d]
    return df


def mean_rank(df):
    """
    Find the mean rank of the matrix. In other words,
    after sorting each column in the matrix, find the 
    mean of each row.
    """
    sorted_dict = {}
    for col in df: 
        sorted_dict.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(sorted_dict)
    return sorted_df.mean(axis=1).tolist()