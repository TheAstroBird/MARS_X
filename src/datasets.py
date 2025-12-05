import pandas as pd
import warnings as wr

def refine(df: pd.DataFrame,
           inertia: str="inertia",
           Love_number: str="Love_number"):
    """Refines the DataFrame from records, where values of MOI or Love number does not fit to observational data

    :param df: DataFrame used
    :type df: pd.DataFrame
    :param inertia: The name of inertia value in the DataFrame df, defaults to 'inertia'
    :type inertia: str
    :param Love_number: The name of Love number value in the DataFrame df, defaults to 'Love_number'
    :type Love_number: str"""
    df_new = df[(0.3634 <= df[inertia]) & (df[inertia] <= 0.3646) &
                (0.166 <= df[Love_number].apply(lambda x: complex(x).real)) &
                (df[Love_number].apply(lambda x: complex(x).real) <= 0.182)]
    if df_new.shape[0] == df.shape[0]:
        wr.warn(message="DataFrame was not changed!!", category=Warning)
    return df_new