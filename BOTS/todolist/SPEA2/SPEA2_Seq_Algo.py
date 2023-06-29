import math
import operator
import random
from Bio.Data import CodonTable
from Bio.SeqUtils import seq3
from todolist.SPEA2.Fitness_Functions import eval_gc_content, eval_hairpins, eval_homopolymers, eval_host, eval_repeats, eval_restriction_sites, eval_splice_sites, eval_start_sites
from todolist.SPEA2.Sequence_Container import SequenceContainer

max_generations = 3


def n_dimensional_euclidian_distance(individual1, individual2):
    assert isinstance(individual2, SequenceContainer) and isinstance(individual1, SequenceContainer)
    gc = (getattr(individual2, "gc_fitness") - getattr(individual1, "gc_fitness")) ** 2
    homo = (getattr(individual2, "homo_fitness") - getattr(individual1, "homo_fitness")) ** 2
    host = (getattr(individual2, "host_fitness") - getattr(individual1, "host_fitness")) ** 2
    repeat = (getattr(individual2, "repeats_fitness") - getattr(individual1, "repeats_fitness")) ** 2
    restriction = (getattr(individual2, "restriction_fitness") - getattr(individual1, "restriction_fitness")) ** 2
    splice = (getattr(individual2, "splice_fitness") - getattr(individual1, "splice_fitness")) ** 2
    start = (getattr(individual2, "start_fitness") - getattr(individual1, "start_fitness")) ** 2
    hairpin = (getattr(individual2, "hairpins_fitness") - getattr(individual1, "hairpins_fitness")) ** 2
    return math.sqrt(gc + homo + host + repeat + restriction + splice + start + hairpin)


def crossover(parents, probability_crossover):
    """number_of_crossovers = len(parents) * probability_crossover
    if number_of_crossovers is 0 or len(parents) < 2 or parents[0] is not SequenceContainer:
        return
    for current_num_crossovers in range(number_of_crossovers):
        parent2 = parent1 = random.randint(1, len(parents))
        while parent1 is parent2:
            parent2 = random.randint(1, len(parents))
        crossover_point = random.randint(1, len(parents-1))
        parents.extend(getattr(parents[parent1], "sequence").tomutable()[:crossover_point].extend(getattr(parents[parent2], "sequence").tomutable()[crossover_point-1:]))
    """
    return parents


def crossover_and_mutate(parents, probability_crossover, probability_mutation, codon_use_table):
    population = []
    for individual in parents:
        population.append(SequenceContainer(
                          mutate_sequence(individual, codon_use_table, probability_mutation))
                          )
    population = crossover(population, probability_crossover)
    return population


def get_sum_fitness(elem):
    return getattr(elem, "gc_fitness") + getattr(elem, "homo_fitness") + getattr(elem,"host_fitness") + getattr(elem, "repeats_fitness") + getattr(elem, "restriction_fitness") + getattr(elem,"splice_fitness") + getattr(elem, "start_fitness") + getattr(elem, "hairpins_fitness")


def get_parents(archive):
    selected = sorted(archive, key=get_sum_fitness, reverse=False)
    return selected[:math.ceil(0.2 * len(archive))]


def get_density(elem):
    return getattr(elem, "density")


def truncate_similar_individuals(archive, archive_size):
    sorted_archive = sorted(archive, key=get_density, reverse=False)
    return sorted_archive[archive_size:]


def populate_archive_with_remaining_best(archive, archive_size, dominated):
    num_individuals_to_add = archive_size - len(archive)
    dominated_set = sorted(dominated, key=lambda sequence:
                           getattr(sequence, "gc_fitness") +
                           getattr(sequence, "homo_fitness") +
                           getattr(sequence, "host_fitness") +
                           getattr(sequence, "repeats_fitness") +
                           getattr(sequence, "restriction_fitness") +
                           getattr(sequence, "splice_fitness") +
                           getattr(sequence, "start_fitness") +
                           getattr(sequence, "hairpins_fitness"), reverse=False)
    dominated_set = dominated_set[:num_individuals_to_add]
    archive.extend(dominated_set)
    return archive


def get_non_dominated_solutions(union):
    """
    Determine pareto dominance, that is
    dominant: solution that has the best scores without betraying the other scores
    non-dominant: solution that betrays the scores, compared to another one
    its a two way relationship
    :param union:
    :return:
    """
    non_dominated_set = []
    dominated_set = []
    pareto_dominance_graph = [[True for x in range(len(union))] for y in range(len(union))]
    # if any vector is larger it is dominated, else it is not dominated
    for idx in range(len(union)):
        for jdx in range(len(union)):
            if getattr(union[idx], "gc_fitness") > getattr(union[jdx], "gc_fitness") or \
                    getattr(union[idx], "homo_fitness") > getattr(union[jdx], "homo_fitness") or \
                    getattr(union[idx], "host_fitness") > getattr(union[jdx], "host_fitness") or \
                    getattr(union[idx], "repeats_fitness") > getattr(union[jdx], "repeats_fitness") or \
                    getattr(union[idx], "restriction_fitness") > getattr(union[jdx], "restriction_fitness") or \
                    getattr(union[idx], "splice_fitness") > getattr(union[jdx], "splice_fitness") or \
                    getattr(union[idx], "start_fitness") > getattr(union[jdx], "start_fitness") or \
                    getattr(union[idx], "hairpins_fitness") > getattr(union[jdx], "hairpins_fitness"):
                pareto_dominance_graph[idx][jdx] = False
    for idx in range(len(union)):
        add = True
        for jdx in range(len(union)):
            if pareto_dominance_graph[idx][jdx] is False:
                add = False
                break
        if add is True:
            non_dominated_set.append(union[idx])
        else:
            dominated_set.append(union[idx])
    return non_dominated_set, dominated_set


def calculate_solution_density(union, nearest_neighbours):
    """
    Evaluate distance between kth nearest neighbours, where k is sqrt(len(union))
    :param nearest_neighbours:
    :param union:
    :return:
    """
    kth = math.ceil(math.sqrt(len(union)))
    for idx in range(len(union)):
        assert isinstance(union[idx], SequenceContainer)
        setattr(union[idx], "density", sorted(nearest_neighbours[idx])[kth])


def calculate_raw_fitness(individual):
    """
    Determine the fitness relative to the others based on how dominant they are
    :param individual:
    :return:
    """
    pass


def calculate_fitness(population, gc_parameters, ancestor_sequence, restriction_sites):
    for individual in population:
        assert isinstance(individual, SequenceContainer)
        setattr(individual, "gc_fitness", eval_gc_content(individual, gc_parameters))
        setattr(individual, "homo_fitness", eval_homopolymers(individual))
        setattr(individual, "host_fitness", eval_host(individual, ancestor_sequence))
        setattr(individual, "repeats_fitness", eval_repeats(individual))
        setattr(individual, "restriction_fitness", eval_restriction_sites(individual, restriction_sites))
        setattr(individual, "splice_fitness", eval_splice_sites(individual))
        setattr(individual, "start_fitness", eval_start_sites(individual))
        setattr(individual, "hairpins_fitness", eval_hairpins(individual))


def initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation, codon_use_table):
    """
    Initializes a population of sequences based on a (hopefully) codon optimized
    sequence
    :param codon_use_table:
    :param probability_mutation: float determines probability of mutation
    :param population_size: int determines minimum population size
    :param problem_size: int helps determines minimum population size
    :param ancestor_sequence: Bio.seq.seq codon optimized sequence
    :return:
    """
    population = []
    mutable_seq = ancestor_sequence.tomutable()
    if problem_size < population_size:
        size = population_size
    else:
        size = problem_size
    for idx in range(size):
        population.append(SequenceContainer(mutable_seq.toseq()))
        mutate_sequence(population[-1], codon_use_table, probability_mutation)


    return population


def mutate_sequence(individual, codon_use_table, mutation_probability=0.05, offset=0):
    """
    Takes a single sequence and gives it a random number of mutations.
    :param individual:
    :param codon_use_table:
    :param mutation_probability:
    :param offset:
    :return:
    """
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()
    num_codons = len(mutable_seq) // 3
    num_mutations = math.ceil(num_codons * mutation_probability)
    for _ in range(num_mutations):
        position = 3 * random.randrange(0, len(mutable_seq) // 3)
        codon_idx = slice(offset + position, (offset + 3) + position)
        new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
        mutable_seq[codon_idx] = new_codon
    return mutable_seq.toseq()


def mutate_codon(codon_in, codon_use_table):
    """Select a synonymous codon in accordance with the frequency of use
    in the host organism.
    Args:
    codon_in (Bio.Seq.Seq): A single codon.
    Returns:
        Bio.Seq.Seq: A new codon.
    """
    amino_acid = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()
    synonymous_codons, codon_use_freq = codon_use_table[amino_acid]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    return codon_out


def determine_neighbours(union):
    neighbours = [[0 for x in range(len(union))] for y in range(len(union))]
    for idx in range(len(union)):
        for jdx in range(len(union)):
            if idx != jdx:
                distance = n_dimensional_euclidian_distance(union[idx], union[jdx])
                neighbours[idx][jdx] = neighbours[jdx][idx] = distance
    return neighbours


def optimize_with_strength_pareto_evolutionary_algorithm(population_size, archive_size, problem_size,
                                                         probability_crossover, probability_mutation,
                                                         ancestor_sequence, codon_use_table, gc_parameters,
                                                         restriction_sites):
    population = initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation,
                                       codon_use_table)
    archive = []
    for generation in range(0, max_generations):
        calculate_fitness(population, gc_parameters, ancestor_sequence, restriction_sites)
        union = population.copy()
        union.extend(archive)
        nearest_neighbours = determine_neighbours(union)
        calculate_solution_density(union, nearest_neighbours)
        for individual in union:
            calculate_raw_fitness(individual)
        archive, dominated = get_non_dominated_solutions(union)
        if len(archive) < archive_size:
            archive = populate_archive_with_remaining_best(archive, archive_size, dominated)
        elif len(archive) > archive_size:
            archive = truncate_similar_individuals(archive, archive_size)
        parents = get_parents(archive)
        population = crossover_and_mutate(parents, probability_crossover, probability_mutation, codon_use_table)

    for item in archive:
        print(getattr(item, "sequence"))
    return archive





