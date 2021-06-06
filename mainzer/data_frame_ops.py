def round_df(df, columns2float={}, columns2int=[]):
    for col, decimals in columns2float.items():
        df[col] = df[col].round(decimals=decimals)
    for col in columns2int:
        df[col] = df[col].astype(int)
    return df