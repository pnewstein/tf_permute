"""
Simulates a neuron picking transcription factors to express at random
"""
from pathlib import Path
import pickle

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

HERE = Path(__file__).parent

def simulate(
    n_cells: int, n_genes: int, prob: float, seed: int, resolution: int
) -> tuple[np.ndarray, float]:
    """
    simulates n_cells turning on n_genes at random reslution times uses a
    random seed and a probability (between 0 and 1) of independently turning on
    each gene returns the number (value) each number of intersections [index]
    occured and the average number of cells with no genes turned on
    """
    rng = np.random.default_rng(seed)
    out = np.zeros(n_cells).astype(np.int64)
    n_no_expressed = 0
    "the number of cells with 0hdtfs"
    for _ in range(resolution):
        rand_bools: np.ndarray = (
            rng.uniform(size=n_genes * n_cells).reshape((n_cells, n_genes)) < prob
        )
        n_no_expressed += sum(rand_bools.sum(axis=1) == 0)
        unique_combos = set((tuple(gene_list) for gene_list in rand_bools))
        n_intersections = rand_bools.shape[0] - len(unique_combos)
        out[n_intersections] += 1
    return out, n_no_expressed / resolution


def plot_simulation_results_on_axs(simulation_results: np.ndarray, ax: plt.Axes):
    """
    makes a figure showing the frequency of each number of intersections
    """
    # x_labels = [f"{i} intersections" if i % 10 == 0 else " " for i in range(len(simulation_results))]
    x_labels = range(len(simulation_results))

    # zoom in on plot so that you are only plotting where data exists
    # cut off begining
    full_columns = simulation_results == 0
    edge_inds = np.nonzero(np.diff(full_columns))[0]
    start_ind = 0 if simulation_results[0] != 0 else edge_inds[0]
    end_ind = -1 if simulation_results[-1] != 0 else edge_inds[-1] + 1
        
    # simulation_results = simulation_results[start_ind: end_ind]
    # x_labels = x_labels[start_ind: end_ind]
    ax.set_xlim(start_ind, end_ind)
    ax.bar(x_labels, simulation_results)


def main():
    # Simulate c. elegans
    c_elegans_path = HERE / "worm_hdtf_simulation_results.pickle"
    if c_elegans_path.exists():
        sim_results, p_cent_0 = pickle.loads(c_elegans_path.read_bytes())
    else:
        sim_results, p_cent_0 = simulate(118, 68, 0.10269192422731803, 2132, 100_000)
        c_elegans_path.write_bytes(pickle.dumps((sim_results, p_cent_0)))
    fig, ax = plt.subplots()
    plot_simulation_results_on_axs(sim_results, ax)
    plt.show()


if __name__ == '__main__':
    main()
