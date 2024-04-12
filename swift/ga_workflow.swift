import files;
import string;
import sys;
import io;
import stats;
import python;
import math;
import location;
import assert;
import R;

import EQPy;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

string resident_work_ranks = getenv("RESIDENT_WORK_RANKS");
string r_ranks[] = split(resident_work_ranks,",");

file model_sh = input(emews_root+"/scripts/acute_mouse.sh");
string model_x = argv("model_x");
int trials = string2int(argv("trials", "10"));
string expected_viral_load_file = argv("viral_load_file");
string deap_cfg_file = argv("deap_cfg");

string json_to_input_line = """
import json
params = json.loads('%s')
line_items = [f'--{k} {v}' for k, v in params.items()]
input_line = ' '.join(line_items)
""";


string result_code_template = """
import pandas as pd
import shutil
import utils

fs = [%s]
df = utils.concat_files(fs)
res = df.groupby('time').mean()
# make 'time' a column rather than index
res = res.reset_index()
res.to_csv('%s', index=False)
j = utils.calculate_j_acute_mice('%s', res)
shutil.rmtree('%s')
""";


app (file out, file err) run(string input_line, string result_file, string instance_dir, int seed)
{
    "bash" model_sh input_line result_file emews_root model_x instance_dir seed @stdout=out @stderr=err;
}


(string score, string log_string) obj(string params, int ga_iter, int param_iter) {
    string results[];
    string input_line_code = json_to_input_line % params;
    string input_line = python_persist(input_line_code, "input_line");
    string instance_id = "instance_%i_%i" % (ga_iter, param_iter);
    string instance = "%s/instances/%s" % (turbine_output, instance_id);
    string mean_result_file = "%s/results/mean_result_%i_%i.csv" % (turbine_output, ga_iter, param_iter);
    make_dir(instance) => {
      foreach i in [0:trials-1:1] {
            file out <instance + "/" + ("%d_out.txt" % i)>;
            file err <instance + "/" + ("%d_err.txt" % i)>;
            string output_file = "%s/result_%d.csv" % (instance, i);
            (out,err) = run(input_line, output_file, instance, i) =>
            results[i] = "'" + output_file + "'";
        }
    }

    string fs = string_join(results, ",");
    string code = result_code_template % (fs, mean_result_file, expected_viral_load_file, instance);
    score = python_persist(code, "str(j)");
    log_string = "%s|%s|%s|%s" % (instance_id, score, params, mean_result_file);
}

(void v) loop (location ME) {
    for (boolean b = true, int i = 1;
       b;
       b=c, i = i + 1)
  {
    // gets the model parameters from the python algorithm
    string params =  EQPy_get(ME);
    boolean c;
    if (params == "DONE")
    {
        string finals =  EQPy_get(ME);
        // TODO if appropriate
        // split finals string and join with "\\n"
        // e.g. finals is a ";" separated string and we want each
        // element on its own line:
        // multi_line_finals = join(split(finals, ";"), "\\n");
        string fname = "%s/final_result" % (turbine_output);
        file results_file <fname> = write(finals) =>
        printf("Writing final result to %s", fname) =>
        // printf("Results: %s", finals) =>
        v = make_void() =>
        c = false;
    }
    else if (params == "EQPY_ABORT")
    {
        printf("EQPy Aborted");
        string why = EQPy_get(ME);
        // TODO handle the abort if necessary
        // e.g. write intermediate results ...
        printf("%s", why) =>
        v = propagate() =>
        c = false;
    }
    else
    {

        string param_array[] = split(params, ";");
        string results[];
        string log[];
        foreach p, j in param_array {
            results[j], log[j]  = obj(p, i, j);
        }

        string res = join(results, ";");
        //printf("passing %s", res);
        string fname = "%s/ga_result_%d.csv" % (turbine_output, i);
        file results_file <fname> = write(join(log, "\n") + "\n");
        EQPy_put(ME, res) => c = true;
    }
  }
}

(void o) start (int ME_rank) {
  location deap_loc = locationFromRank(ME_rank);
    EQPy_init_package(deap_loc,"deap_ga") =>
    EQPy_get(deap_loc) =>
    EQPy_put(deap_loc, deap_cfg_file) =>
      loop(deap_loc) => {
        EQPy_stop(deap_loc);
        o = propagate();
    }
}

// deletes the specified directory
app (void o) rm_dir(string dirname) {
  "rm" "-rf" dirname;
}

// call this to create any required directories
app (void o) make_dir(string dirname) {
  "mkdir" "-p" dirname;
}

// anything that need to be done prior to a model runs
// (e.g. file creation) can be done here
//app (void o) run_prerequisites() {
//
//}


main() {
  int rank = string2int(r_ranks[0]);
  start(rank);
}
