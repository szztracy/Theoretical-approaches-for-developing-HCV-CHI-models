import pandas as pd
import json
import numpy as np


class EclipseMinMaxConstraint :

    def __init__(self, ga_params, min_name, max_name):
        for i, param in enumerate(ga_params):
            if param.name == min_name:
                self.min_idx = i
                self.min_param = param
            elif param.name == max_name:
                self.max_idx = i
                self.max_param = param

    def check_constraint(self, individual):
        min_val = individual[self.min_idx]
        max_val = individual[self.max_idx]
        return max_val - min_val >= 1.01

    def redraw(self, individual):
        """ Redraws the constrained params until satisfied
        """
        while not self.check_constraint(individual):
            individual[self.min_idx] = self.min_param.randomDraw()
            individual[self.max_idx] = self.max_param.randomDraw()

    def remutate(self, original_individual, mutated_individual, mu, indpb):
        """ Update the mutated individual using the values from the original
            individual until constraint is satisfied
        """

        while not self.check_constraint(mutated_individual):
            mutated_individual[self.min_idx] = self.min_param.mutate(original_individual[self.min_idx], mu=mu, 
                indpb=indpb)
            mutated_individual[self.max_idx] = self.max_param.mutate(original_individual[self.max_idx], mu=mu, 
                indpb=indpb)


def calculate_j_v2(viral_loads_file, results_file, na_log_cutoff):
    print("Calc for V2")
    header = ['susceptible','eclipse','infected','infecVirusLoad','vp1','vp2']
    rdf = pd.read_csv(results_file, names=header)
    rdf = rdf.replace("-Infinity", -np.inf)
    rdf = rdf.astype({'infecVirusLoad' : "float64"})
    rdf['hour'] = np.arange(len(rdf)) + 1

    vl = pd.read_csv(viral_loads_file)
    vl['hour'] =vl['Day'] * 24
    vl['result'] = rdf[rdf['hour'].isin(vl['hour'].values)]['infecVirusLoad'].values
    # where vl['viral_load'] is -9999 (-9999 is flag for NA) if vl['result] is < 3.6 then return 
    # vl['result'] - vl['viral_load'] = 0, else treat vl['viral_load'] and do as normal
    vl['diff'] = 0
    vl.loc[vl['viral_load'] != -9999, 'diff'] = vl['result'] - vl['viral_load']
    vl.loc[(vl['result'] > na_log_cutoff) & (vl['viral_load'] == -9999), 'diff' ] = vl['result'] - na_log_cutoff
    
    j = np.sqrt(np.sum(np.square(vl['diff'])) / (vl.shape[0] * np.square(0.5)))
    if np.isposinf(j) or np.isneginf(j):
        j = 999999
    return j


def calculate_j_acute_mice(viral_loads_file, rdf):
    vl = pd.read_csv(viral_loads_file)
    # vl = vl[vl.hour >= 1.0]
    vl['result'] = rdf[rdf['time'].isin(vl['hour'].values)]['HCV'].values
    vl['diff'] = vl['result'] - vl['viral_load']

    j = np.sqrt(np.sum(np.square(vl['diff'])) / (vl.shape[0] * np.square(0.5)))
    if np.isposinf(j) or np.isneginf(j):
        j = 999999
    return j


def calculate_j(viral_loads_file, results_file, na_log_cutoff):
    rdf = pd.read_csv(results_file)
    rdf['log_iv_load'] = 0
    idxs = rdf.index[rdf['infecVirusLoad'] > 1]
    c1 = rdf.columns.get_loc('infecVirusLoad')
    c2 = rdf.columns.get_loc('log_iv_load')
    rdf.iloc[idxs, c2] = np.log10(rdf.iloc[idxs, c1])
    # first row is end of first hour
    rdf['hour'] = np.arange(len(rdf)) + 1

    vl = pd.read_csv(viral_loads_file)
    vl['hour'] =vl['Day'] * 24
    vl['result'] = rdf[rdf['hour'].isin(vl['hour'].values)]['log_iv_load'].values
    # where vl['viral_load'] is -9999 (-9999 is flag for NA) if vl['result] is < 3.6 then return 
    # vl['result'] - vl['viral_load'] = 0, else treat vl['viral_load'] and do as normal
    vl['diff'] = 0
    vl.loc[vl['viral_load'] != -9999, 'diff'] = vl['result'] - vl['viral_load']
    vl.loc[(vl['result'] > na_log_cutoff) & (vl['viral_load'] == -9999), 'diff' ] = vl['result'] - na_log_cutoff
    
    return np.sqrt(np.sum(np.square(vl['diff'])) / (vl.shape[0] * np.square(0.5)))


def json_to_input_file(json_string, output_file):
    params = ["gridWidth","gridHeight", "initialHepatocytes","initialViralLoad", "eclipsePhase_MIN",
                "eclipsePhase_MAX", "viralProdSteadyState", "viralProdSteepnessGrowthCurve", "viralProdCycleMidPoint",
                "viralProdCycleConstant", "viralProdCycleDecayExponent", "viralClearanceRate", "viralInfectionRate",
                "infectiousVirusProportion", "infectiousVirusProportion_firstday", "chg_infectiousVirusProportion",
                "virusProdAmplificationFactor", "timeAtInfVirusPropChanges", "timeAtViralProdDrops"]
    
    j = json.loads(json_string)
    values = []
    with open(output_file, 'w') as f_out:
        f_out.write('Param,Value\n')
        for p in params:
            f_out.write('{},{}\n'.format(p, j[p]))
            values.append(str(j[p]))
    return values


def concat_files(fs, ts=None):
    dfs = []
    for f in fs:
        df = pd.read_csv(f)
        if ts is not None:
            df = df[df['time'].isin(ts)]
        dfs.append(df)

    df = pd.concat(dfs)
    return df


def within_range(es_files, is_files, es_ts, is_ts, lb_targets, ub_targets):
    df = concat_files(es_files, es_ts)
    ess = df.groupby('time').mean()
    ess['group'] = ['e{}_ni'.format(int(x / 24)) for x in ess.index]
    ess = ess.set_index('group')
    ess.rename(columns={'HCV': 'val'}, inplace=True)


    df = concat_files(is_files, is_ts)
    iss = df.groupby('time').mean()
    iss['group'] = ['i{}_ni'.format(int(x / 24)) for x in iss.index]
    iss = iss.set_index('group')
    iss.rename(columns={'percentage': 'val'}, inplace=True)
    ess_iss = pd.concat([ess, iss])

    lbd = json.loads(lb_targets)
    lbd = {k: v for k, v in lbd.items() if k.endswith('_ni')}
    lb = pd.Series(lbd, name="current_lower_bounds")
    ubd = json.loads(ub_targets)
    ubd = {k: v for k, v in ubd.items() if k.endswith('_ni')}
    ub = pd.Series(ubd, name="current_upper_bounds")
    bounds = pd.concat((lb, ub), axis=1)
    df = bounds.assign(val=ess_iss)

    # print(df, flush=True)
    res = df.val.between(df.current_lower_bounds, df.current_upper_bounds, inclusive='both')
    # print(res, flush=True)
    results = [str(x) for x in df.val]
    # print(results)
    return (res.all(), results)
