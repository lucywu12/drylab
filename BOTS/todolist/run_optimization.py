import threading
import time
from Bio.Alphabet import IUPAC
from todolist.SPEA2.Codon_Use import *
from todolist.SPEA2.Bio_Structures import GC_content, RestrictionEnzymes
from todolist.SPEA2.SPEA2_Seq_Algo import optimize_with_strength_pareto_evolutionary_algorithm


class ThreadingExample(object):
    """ Threading example class
    The run() method will be started and it will run in the background
    until the application exits.
    """

    def __init__(self, population_size, archive_size, problem_size,
                 probability_crossover, probability_mutation,
                 ancestor_sequence, codon_use_table, gc_parameters,
                 restriction_sites, interval=1):
        print("threading example")
        self.interval = interval

        thread = threading.Thread(target=self.run, args=(population_size, archive_size, problem_size,
                                                         probability_crossover, probability_mutation,
                                                         ancestor_sequence, codon_use_table, gc_parameters,
                                                         restriction_sites))
        thread.daemon = True  # Daemonize thread
        thread.start()  # Start the execution

    def run(self, population_size, archive_size, problem_size,
            probability_crossover, probability_mutation,
            ancestor_sequence, codon_use_table, gc_parameters,
            restriction_sites):
        print("run")
        optimize_with_strength_pareto_evolutionary_algorithm(population_size, archive_size, problem_size,
                                                             probability_crossover, probability_mutation,
                                                             ancestor_sequence, codon_use_table, gc_parameters,
                                                             restriction_sites)
        print('completed')
