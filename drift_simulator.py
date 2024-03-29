import numpy as np
import argparse

import pdb


class StoreAsArray(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run drift simulator with parameters\
        intitial frequencies or step (step default 0.1), total_pop (default 96),\
        no_males (default 48), no_females (default 48), cycles (default 1), sims (default 1),\
        if plots should be created (default False) and output file name (default drift_simulation.csv)"
    )
    parser.add_argument(
        "-i",
        "--initial",
        action=StoreAsArray,
        type=float,
        nargs="+",
        help="Initial allele frequency/ies",
    )
    parser.add_argument(
        "-t",
        "--total_pop",
        type=int,
        default=96,
        help="Population size in each generation (number of progeny)",
    )
    parser.add_argument(
        "-m",
        "--no_males",
        type=int,
        default=48,
        help="Number of males contributing gametes",
    )
    parser.add_argument(
        "-f",
        "--no_females",
        type=int,
        default=48,
        help="Number of females contributing gametes",
    )
    parser.add_argument(
        "-c", "--cycles", default=1, help="Number of generations of random mating"
    )
    parser.add_argument(
        "-p",
        "--plot",
        default=False,
        help="True/False, create a plot of allele frequencies for each simulation",
    )
    parser.add_argument(
        "-s", "--sims", default=1, type=int, help="Number of simulations to perform"
    )
    parser.add_argument(
        "-e",
        "--step",
        type=float,
        help="Step between allele frequences. E.g., the default of 0.1 results in simulating drift for starting allele frequencies of 0, 0.1, 0.2,...1",
    )
    parser.add_argument(
        "-q",
        "--quantiles",
        action=StoreAsArray,
        type=float,
        nargs="+",
        help="Array of quantiles to calculate, defaults to .1, .9 ... .000001, .999999 (12 total)",
        default=np.array(
            [
                0.1,
                0.9,
                0.01,
                0.99,
                0.001,
                0.999,
                0.0001,
                0.9999,
                0.00001,
                0.99999,
                0.000001,
                0.999999,
            ]
        ),
    )
    parser.add_argument(
        "-r",
        "--run_quantiles",
        type=bool,
        default=True,
        help="Whether or not to calcualte quantiles as given by quantiles.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        default="drift_simulation",
        help="Prefix for files to write output to.",
    )

    return parser


def drift_neutral(
    initial: float,
    total_pop: int,
    no_males: int,
    no_females: int,
    cycles: int,
    sims: int,
) -> np.array:
    """
    Python implementation of DriftSimulator (R)
    by Timothy M Beissinger, 2012-01-25.
    Simulates drift for a neutral locus.
    Args:
        initial (numeric, 1>=n>=0):
            Assumed allele frequency in the initial population
        total_pop (int): population size in each generation (number of progeny)
        no_males (int): Number of males contributing gametes
        no_females (int): Number of females contributing gametes
        cycles (int): Number of generations of random mating
        sims (int): Number of simulations to perform
    Returns:
        A vector with one element for each simulation performed.
        Each element of the returned vector represents the final
        allele frequency for one simulation.
    """
    from math import trunc

    float_ones = initial * total_pop * 2
    trunc_ones = trunc(float_ones)

    initial1 = np.ones(trunc_ones)
    initial0 = np.zeros((2 * total_pop) - trunc_ones)
    # array to hold current alleles in population
    freq_array = np.concatenate((initial1, initial0))

    # final_frequencies stores only the final cycle's
    # allele frequencies. This makes sense when considering
    # > 1 simulation. Otherwise, it is identical to the later-declared
    # gen_freq
    final_frequencies = np.zeros(sims)

    for s in range(sims):
        # create an array to hold the allele frequency
        # at the end of each cycle
        gen_freq = np.zeros(cycles + 1)
        # cycle 0 is simply the initial allele frequency
        gen_freq[0] = initial

        for c in range(cycles):
            m_prog = np.random.choice(freq_array, no_males)
            m_prog = np.resize(m_prog, total_pop)
            f_prog = np.random.choice(freq_array, no_females)
            f_prog = np.resize(f_prog, total_pop)
            gen_freq[c + 1] = (np.sum(m_prog) + np.sum(f_prog)) / (2 * total_pop)
        final_frequencies[s] = gen_freq[c + 1]

    return final_frequencies


def calculate_quantiles(
    frequency_array: np.array, quantile_array: np.array
) -> np.array:
    # once upgrade numpy, add in method option
    return np.quantile(
        frequency_array,
        quantile_array,
        # method="median_unbiased",
    )


def write_matrix(prefix: str, suffix: str, matrix: np.array):
    file_name = "{}_{}.txt.gz".format(prefix, suffix)
    np.savetxt(file_name, matrix, fmt="%f")


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.initial is None and not args.step:
        raise ValueError(
            "Either intial frequencies or an interval step need to be provided."
        )

    if args.step:
        num = (1 / args.step) + 1
        initial = np.linspace(0, 1, num)
    else:
        initial = args.initial

    # set dimensions for results matrix
    # the first column of the matrix will be the intial frequencies
    # largely for clarity.
    mat_n = initial.shape[0]
    mat_m = args.sims + 1
    results_matrix = np.zeros((mat_n, mat_m))
    results_matrix[:, 0] = initial

    qs = args.quantiles.shape[0]
    quantile_matrix = np.zeros((mat_n + 1, qs))
    quantile_matrix[0] = args.quantiles

    row = 0
    for freq in np.nditer(initial):
        final_sim_frequencies = drift_neutral(
            freq, args.total_pop, args.no_males, args.no_females, args.cycles, args.sims
        )
        results_matrix[row, 1:] = final_sim_frequencies
        if args.run_quantiles:
            qts = calculate_quantiles(final_sim_frequencies, args.quantiles)
            quantile_matrix[row + 1] = qts
        row += 1

    write_matrix(args.outfile, "sim_frequencies", results_matrix)

    if args.run_quantiles:
        write_matrix(args.outfile, "quantiles", quantile_matrix)

# Several things to add:
# 1) should un-randomize, make sure get same values returned
# given the same frequency (do for multiple frequencies so can feel safer)
# 2) Should start working on unit tests. Should factor out function.
# 3) Add in plotting abilities (nice to have)
