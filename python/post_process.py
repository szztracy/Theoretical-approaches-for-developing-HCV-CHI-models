import pandas as pd

def concat(row):
    df = pd.read_csv(row.output_file)
    df['ga_iter'] = row.ga_iter
    df['params_id'] = row.param_iter
    df['run'] = row.run
    return df

def create_topn_file(result_file, top_n, top_n_file):
    header = ['param_score', 'ga_iter', 'param_iter', 'run', 'mse', 'params', 'output_file']
    df = pd.read_csv(result_file, delimiter='|', names=header)
    df = df.sort_values(['param_score'])
    df = df.reset_index()
    top10 = top_n * 10
    dfs = df.iloc[0: top10].apply(concat, axis=1)
    pd.concat(dfs.to_list(), axis=0).to_csv(top_n_file, index=False)
