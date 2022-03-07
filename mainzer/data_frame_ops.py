import pandas as pd


def round_df(
    df: pd.DataFrame,
    columns2float: dict={},
    columns2int: dict=[]
):
    for col, decimals in columns2float.items():
        df[col] = df[col].round(decimals=decimals)
    for col in columns2int:
        df[col] = df[col].astype(int)
    return df