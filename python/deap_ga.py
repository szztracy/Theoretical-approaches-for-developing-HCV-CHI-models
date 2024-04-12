import yaml
import random
import time
import math
import csv
import json
import time
import os
import pickle

import numpy as np

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import eqpy, ga_utils
import utils

# list of ga_utils parameter objects
ga_params = None
# ec_constraint = None

def eaSimple(population, toolbox, cxpb, mutpb, ngen, stats=None,
             halloffame=None, verbose=True, pkl_file=None):
    """This algorithm reproduce the simplest evolutionary algorithm as
    presented in chapter 7 of [Back2000]_.

    :param population: A list of individuals.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The probability of mating two individuals.
    :param mutpb: The probability of mutating an individual.
    :param ngen: The number of generation.
    :param stats: A :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: Whether or not to log the statistics.
    :returns: The final population
    :returns: A class:`~deap.tools.Logbook` with the statistics of the
              evolution

    The algorithm takes in a population and evolves it in place using the
    :meth:`varAnd` method. It returns the optimized population and a
    :class:`~deap.tools.Logbook` with the statistics of the evolution. The
    logbook will contain the generation number, the number of evaluations for
    each generation and the statistics if a :class:`~deap.tools.Statistics` is
    given as argument. The *cxpb* and *mutpb* arguments are passed to the
    :func:`varAnd` function. The pseudocode goes as follow ::

        evaluate(population)
        for g in range(ngen):
            population = select(population, len(population))
            offspring = varAnd(population, toolbox, cxpb, mutpb)
            evaluate(offspring)
            population = offspring

    As stated in the pseudocode above, the algorithm goes as follow. First, it
    evaluates the individuals with an invalid fitness. Second, it enters the
    generational loop where the selection procedure is applied to entirely
    replace the parental population. The 1:1 replacement ratio of this
    algorithm **requires** the selection procedure to be stochastic and to
    select multiple times the same individual, for example,
    :func:`~deap.tools.selTournament` and :func:`~deap.tools.selRoulette`.
    Third, it applies the :func:`varAnd` function to produce the next
    generation population. Fourth, it evaluates the new individuals and
    compute the statistics on this population. Finally, when *ngen*
    generations are done, the algorithm returns a tuple with the final
    population and a :class:`~deap.tools.Logbook` of the evolution.

    .. note::

        Using a non-stochastic selection method will result in no selection as
        the operator selects *n* individuals from a pool of *n*.

    This function expects the :meth:`toolbox.mate`, :meth:`toolbox.mutate`,
    :meth:`toolbox.select` and :meth:`toolbox.evaluate` aliases to be
    registered in the toolbox.

    .. [Back2000] Back, Fogel and Michalewicz, "Evolutionary Computation 1 :
       Basic Algorithms and Operators", 2000.
    """
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream, flush=True)

    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))

        # Vary the pool of individuals
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        # Replace the current population by the offspring
        population[:] = offspring
        if pkl_file is not None:
            pf = f'{pkl_file}_{gen}.pkl'
            with open(pf, 'wb') as fout:
                pickle.dump([population, halloffame], fout)

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream, flush=True)

    return population, logbook


def obj_func(x):
    return 0

# {"batch_size":512,"epochs":51,"activation":"softsign",
#"dense":"2000 1000 1000 500 100 50","optimizer":"adagrad","drop":0.1378,
#"learning_rate":0.0301,"conv":"25 25 25 25 25 1"}
def create_list_of_json_strings(list_of_lists, super_delim=";"):
    # create string of ; separated jsonified maps
    res = []
    global ga_params
    for l in list_of_lists:
        jmap = {}
        for i,p in enumerate(ga_params):
            jmap[p.name] = l[i]

        jstring = json.dumps(jmap)
        res.append(jstring)

    return (super_delim.join(res))

def create_fitnesses(params_string):
    """return equivalent length tuple list
    :type params_string: str
    """
    params = params_string.split(";")
    # get length
    res = [(i,) for i in range(len(params))]
    return (res)

def queue_map(obj_func, pops):
    # Note that the obj_func is not used
    # sending data that looks like:
    # [[a,b,c,d],[e,f,g,h],...]
    if not pops:
        return []
    # print("OUT: population")
    eqpy.OUT_put(create_list_of_json_strings(pops))
    # print("WAITING for IN_get()")
    result = eqpy.IN_get()
    # print("IN_get complete")
    split_result = result.split(';')
    # TODO determine if max'ing or min'ing and use -9999999 or 99999999
    return [(float(x),) if not math.isnan(float(x)) else (float(99999999),) for x in split_result]
    #return [(float(x),) for x in split_result]

def make_random_params():
    """
    Performs initial random draw on each parameter
    """
    global ga_params
    individual = [p.randomDraw() for p in ga_params]
    # redraws until the constraint is satisfied
    # ec_constraint.redraw(individual)
    return individual

def parse_init_params(params_file):
    init_params = []
    with open(params_file) as f_in:
        reader = csv.reader(f_in)
        header = next(reader)
        for row in reader:
            init_params.append(dict(zip(header, row)))
    return init_params

def update_init_pop(pop, params_file, fraction_of_pop):
    global ga_params
    if fraction_of_pop > 0:
        print("Reading initial population from {}".format(params_file))
        init_params = parse_init_params(params_file)
        count = int(len(pop) * fraction_of_pop)
        if count > len(init_params):
            raise ValueError("Not enough initial params to set the population: size of init params < population size * fraction_of_pop")

        print("Replacing {} individuals in random population with individuals from {}".format(count, params_file))
        for i, indiv in enumerate(pop[:count]):
            for j, param in enumerate(ga_params):
                indiv[j] = param.parse(init_params[i][param.name])
            
# keep as reference for log type
# def mutGaussian_log(x, mu, sigma, mi, mx, indpb):
#     if random.random() < indpb:
#         logx = math.log10(x)
#         logx += random.gauss(mu, sigma)
#         logx = max(mi, min(mx, logx))
#         x = math.pow(10, logx)
#     return x

# Returns a tuple of one individual
def custom_mutate(individual, indpb):
    """
    Mutates the values in list individual with probability indpb
    """

    # Note, if we had some aggregate constraint on the individual
    # (e.g. individual[1] * individual[2] < 10), we could copy
    # individual into a temporary list and mutate though until the
    # constraint was satisfied
    original_values = list(individual)
    for i, param in enumerate(ga_params):
        individual[i] = param.mutate(original_values[i], mu=0, indpb=indpb)
    # ec_constraint.remutate(original_values, individual, mu=0, indpb=indpb)

    return individual,

def cxUniform(ind1, ind2, indpb):
    c1, c2 = tools.cxUniform(ind1, ind2, indpb)
    # while (not ec_constraint.check_constraint(c1)) or (not ec_constraint.check_constraint(c2)):
    #    c1, c2 = tools.cxUniform(ind1, ind2, indpb)
    return (c1, c2)

def timestamp(scores):
    return str(time.time())

def run():
    eqpy.OUT_put("Params")
    param_file = eqpy.IN_get()
    with open(param_file) as f_in:
        params = yaml.safe_load(f_in)
    
    num_iter = params['num_iter']
    num_pop = params['num_pop']
    seed = params['seed']
    strategy = params['strategy']
    mut_prob = params['mut_prob']
    ga_params_file = params['ga_params_file']
    # init_pop_file = params['init_pop_file']
    # init_pop_replacement_fraction = params['init_pop_replacement_fraction']
    
    random.seed(seed)
    global ga_params
    ga_params = ga_utils.create_parameters(ga_params_file)
    # global ec_constraint
    # ec_constraint = utils.EclipseMinMaxConstraint(ga_params, 'epmi', 'epmx')
   
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)
    toolbox = base.Toolbox()
    toolbox.register("individual", tools.initIterate, creator.Individual,
                     make_random_params)

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", obj_func)
    toolbox.register("mate", cxUniform, indpb=0.5)
    mutate_indpb = mut_prob
    toolbox.register("mutate", custom_mutate, indpb=mutate_indpb)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("map", queue_map)

    pop = toolbox.population(n=num_pop)
    # update_init_pop(pop, init_pop_file, init_pop_replacement_fraction)
   
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)
    stats.register("ts", timestamp)

    # num_iter-1 generations since the initial population is evaluated once first
    mutpb = mut_prob
    start_time = time.time()
    if strategy == 'simple':
        pkl_file = f'{os.getenv("TURBINE_OUTPUT")}/GA'
        print(pkl_file, flush=True)
        pop, log = eaSimple(pop, toolbox, cxpb=0.5, mutpb=mutpb, ngen=num_iter - 1,
                            stats=stats, halloffame=hof, verbose=True, pkl_file=pkl_file)
    elif strategy == 'mu_plus_lambda':
        mu = int(math.floor(float(num_pop) * 0.5))
        lam = int(math.floor(float(num_pop) * 0.5))
        if mu + lam < num_pop:
            mu += num_pop - (mu + lam)

        pop, log = algorithms.eaMuPlusLambda(pop, toolbox, mu=mu, lambda_=lam,
                                             cxpb=0.5, mutpb=mutpb, ngen=num_iter - 1,
                                             stats=stats, halloffame=hof, verbose=True)
    else:
        raise NameError('invalid strategy: {}'.format(strategy))

    end_time = time.time()

    fitnesses = [str(p.fitness.values[0]) for p in pop]

    eqpy.OUT_put("DONE")
    # return the final population
    eqpy.OUT_put("{}\n{}\n{}\n{}\n{}".format(create_list_of_json_strings(pop), ';'.join(fitnesses),
        start_time, log, end_time))
