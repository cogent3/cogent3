#!/usr/bin/env python
"""Algorthims for function optimisation

   great_deluge() is a hillclimbing algorithm based on:    
      Gunter Dueck: New Optimization Heuristics, The Great Deluge Algorithm
      and the Record-to-Record Travel. Journal of Computational Physics, Vol.
      104, 1993, pp. 86 - 92

   ga_evolve() is a basic genetic algorithm in which all internal functions can
      be overridden

   NOTE: both optimisation functions are generators.
"""

from numpy.random import normal

__author__ = "Daniel McDonald and Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Production"

def _simple_breed(best, num, mutation_rate, random_f):
    """Returns num copies of parent with mutation_rate changes"""
    result = []
    score, parent = best
    for child_number in range(num):
        if random_f() <= mutation_rate:
            child = parent.mutate()
            result.append(child)
        else:
            result.append(parent)
    return result

def _simple_score(child, target):
    """Returns the childs score as defined by the childs scoring function"""
    return child.score(target)
    
def _simple_init(parent, num):
    """Creates a list parent copies"""
    return [parent.copy() for i in range(num)]

def _simple_select(population, scores):
    """Returns a tuple: (best_score, best_child)"""
    scored = zip(scores, population)
    scored.sort()
    return scored[0]

def great_deluge(a, step_factor=500, max_iter=100, max_total_iters=1000):
    """This generator makes random variations of the object a to minimize cost.

       Yields are performed at the end of each iteration and a tuple containing
       ((iter_count, total_iters), a) is returned. iter_count is used to
       kill the while loop in the event that no new objects are found with a
       better cost. iter_count gets reset each time an object with a better
       cost is found. total_iters will kill the while loop when the total
       number of iterations through the loop reaches max_total_iters

       Object a must implement methods cost() and perturb() for evaluating
       the score and making mutations respectively. Usually, you'll want to
       write a wrapper that passes these through to methods of an internal
       data object, or functions acting on that object.
    """
    water_level = curr_cost = a.cost() # can't be worse than initial guess
    step_size = abs(water_level)/step_factor
    iter_count = 0
    total_iters = 0
    while iter_count < max_iter and total_iters < max_total_iters:
        new = a.perturb()
        new_cost = new.cost()
        if new_cost < water_level:
            if new_cost < curr_cost:
                water_level = max(curr_cost, water_level - step_size)
                iter_count = 0      # WARNING: iter_count is reset here!
            curr_cost = new_cost
            a = new
        else:
            iter_count += 1
        yield ((iter_count, total_iters), a)
        total_iters += 1

def ga_evolve(parent, target, num, mutation_rate=0.01, score_f=_simple_score, 
              breed_f=_simple_breed, select_f=_simple_select, 
              init_f=_simple_init, random_f=normal, max_generations=1000):
    """Evolves a population based on the parent to the target

       Parent must implement methods copy(), mutate(), and score(target) to be
       used with the simple default functions.
    
       Yields are performed at the end of each iteration and contain the tuple
       (generation, best). The default functions return the tuple
       (generation, (best_score, best_obj)).

       Arguments:
         parent:          Object to create initial population from.
         target:          The goal of the evolution.
         num:             Population size.
         mutation_rate:   Rate at which objects in the population are mutated.
         score_f:         Function to score the object against the target.
         breed_f:         Function to create new population with mutations
         select_f:        Function to select best object(s) from the population
         random_f:        Function to be used in breed_f
         max_generations: Kills while loop if max_generations is reached

       Overload default functions:
         score_f:   Must take an object and a target score. Returns objects 
                    score.
         breed_f:   Must take a tuple containing (scores, objects), the size of 
                    population, a mutation rate and random function to use. 
                    Returns a list containing the initial population. Default 
                    function takes only the best object, but this may not be 
                    desired behavior.
         select_f:  Must take a population and scores. Returns a tuple 
                    containing the best scores and objects in the population. 
                    Default function returns only the best score and object.
         init_f:    Must take an object and the size of the population. Returns
                    a list containing the starting population
    """
    generation = 0
    population = init_f(parent, num)
    while generation < max_generations:
        scores = [score_f(child, target) for child in population]
        best = select_f(population, scores)
        population = breed_f(best, num, mutation_rate, random_f)
        yield (generation, best)
        generation += 1
