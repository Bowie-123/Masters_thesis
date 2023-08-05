import random
import numpy as np
from math import exp
import pandas as pd
import os

class Daphnia:
    def __init__(self, resistance):
        self.resistance = resistance
        self.fitness = 1
        self.parasite = None

    def update_fitness(self, gamma, parasite_virulence):
        self.fitness = max(0, 1 - parasite_virulence * (1 - self.resistance) - gamma * self.resistance)

class Parasite:
    def __init__(self, gene_A, gene_B):
        self.gene_A = gene_A
        self.gene_B = gene_B
        self.fitness = 1
        self.virulence = 0 
        self.host = None

    def calculate_virulence(self, scarcity):
        return 1 / (1 + exp(self.gene_A * (-scarcity + self.gene_B)))

    def update_fitness(self, h, virulence, k, host_resistance, scarcity, c):
        if self.host is not None:
            self.fitness = max(0, 1 + h * virulence - k * host_resistance + (1-scarcity) + scarcity * virulence - c * virulence)


def scarcityfunc(number_gens, alpha=None, Vspeed=None):
    scarcitylist = []
    
    S_t = 0.5
    for i in range(number_gens):
        if i % Vspeed == 0:
            Urandom = random.random()
            S_t = S_t * alpha + ((1 - alpha) * Urandom)
        scarcitylist.append(S_t)
    return scarcitylist[:number_gens]

def create_population(num_daphnia, num_parasites):
    daphnia_pop = [Daphnia(np.random.normal(0.5, 0.1)) for _ in range(num_daphnia)]
    parasite_pop = [Parasite(np.random.normal(0.0, 0.1), np.random.normal(0.5, 0.1)) for _ in range(num_parasites)]
    return daphnia_pop, parasite_pop 

def couple_symbionts(daphnia_pop, parasite_pop):
    paired_parasites = random.sample(parasite_pop, len(daphnia_pop))
    free_parasites = [parasite for parasite in parasite_pop if parasite not in paired_parasites]

    for i, daphnia in enumerate(daphnia_pop):
        daphnia.parasite = paired_parasites[i]
        paired_parasites[i].host = daphnia

    return daphnia_pop, paired_parasites, free_parasites

def run_interaction(daphnia_pop, paired_parasites, free_parasites, gamma, k, h, c, scarcity):
    
    for parasite in paired_parasites + free_parasites:
        parasite.virulence = parasite.calculate_virulence(scarcity)
                
    avg_virulence = np.mean([parasite.virulence for parasite in paired_parasites + free_parasites])
 
    for daphnia in daphnia_pop:
        if daphnia.parasite is not None:
            daphnia.update_fitness(gamma, daphnia.parasite.virulence)
    
    for parasite in paired_parasites:
        if parasite.host is not None:
            parasite.update_fitness(h, parasite.virulence, k, parasite.host.resistance, scarcity, c)
    
    avg_daph_fitness = np.mean([daphnia.fitness for daphnia in daphnia_pop])
    avg_parasite_fitness = np.mean([parasite.fitness for parasite in paired_parasites + free_parasites])
    
    return daphnia_pop, paired_parasites, free_parasites, avg_virulence, avg_daph_fitness, avg_parasite_fitness

def make_cumulfitlist(population):
    pop = [org.fitness for org in population]
    cumul_fitlist = np.cumsum(pop)
    return cumul_fitlist

def get_parent(cumul_fitlist):
    target = random.uniform(0, cumul_fitlist[-1])
    for i, cumul_fit in enumerate(cumul_fitlist):
        if cumul_fit > target:
            return i
    return len(cumul_fitlist) - 1

def mutate_gene(gene, organism_type=None):
    mutation_rate = 0.01

    if random.random() < mutation_rate:
        
        mutation_step = np.random.normal(0, 0.1)   # Gaussian mutation with mean 0 and stdev 0.1
        # mutation_step = random.choice([-0.05, 0.05])
        gene += mutation_step

        # Ensure the resistance gene R of the daphnia does not go above 1 or below 0
        if organism_type == "daphnia":
            if gene > 1:
                gene = 1
            elif gene < 0:
                gene = 0

    return gene

def reproduce(population, organism_type):
    new_pop = []
    cumul_fitlist = make_cumulfitlist(population)
    
    for _ in range(len(population)):
        parent_index = get_parent(cumul_fitlist)
        parent = population[parent_index]

        if organism_type == "daphnia":
            new_resistance = mutate_gene(parent.resistance, organism_type)
            offspring = Daphnia(new_resistance)
        elif organism_type == "parasite":
            new_gene_A = mutate_gene(parent.gene_A)
            new_gene_B = mutate_gene(parent.gene_B)
            offspring = Parasite(new_gene_A, new_gene_B)

        new_pop.append(offspring)

    return new_pop

def main(params):
    all_results = []
     
    for run in range(params["number_runs"]):
        
        total_generations = params["number_gens"]               
        scarcity_list = params["scarcity_function"](params["number_gens"], params["alpha"], params["Vspeed"])
        daphnia_pop, parasite_pop = create_population(params["num_daphnia"], params["num_parasites"])
        daphnia_pop, paired_parasites, free_parasites = couple_symbionts(daphnia_pop, parasite_pop)

        generation_results = []
        
        for generation in range(total_generations):
                        
            scarcity = scarcity_list[generation]
            daphnia_pop, paired_parasites, free_parasites, avg_virulence, avg_daph_fitness, avg_parasite_fitness = run_interaction(daphnia_pop, paired_parasites, free_parasites, params["gamma"], params["k"], params["h"], params["c"], scarcity)

            avg_R = np.mean([daphnia.resistance for daphnia in daphnia_pop])
            avg_A = np.mean([parasite.gene_A for parasite in paired_parasites + free_parasites])
            avg_B = np.mean([parasite.gene_B for parasite in paired_parasites + free_parasites])

            generation_results.append({
                "avg_R": avg_R,
                "avg_virulence": avg_virulence,
                "avg_A": avg_A,
                "avg_B": avg_B,
                "scarcity": scarcity,
                "avg_daph_fitness": avg_daph_fitness,
                "avg_parasite_fitness": avg_parasite_fitness
            })

            parasite_pop = paired_parasites + free_parasites
            daphnia_pop = reproduce(daphnia_pop, "daphnia")
            parasite_pop = reproduce(parasite_pop, "parasite")
        
            daphnia_pop, paired_parasites, free_parasites = couple_symbionts(daphnia_pop, parasite_pop)

        all_results.append(generation_results)
              
    return all_results
