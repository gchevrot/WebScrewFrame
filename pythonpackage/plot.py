import matplotlib.pyplot as plt

def plot_distances(data, figsize=(8,6), font_size=16):
    """
    Plot a histogram of the distances between atoms
    """
    chains = set(data.index.values)

    fig, ax = plt.subplots(figsize = figsize)
    ax.set_xlabel(r"Distances [$\AA$]", size=font_size)
    ax.set_ylabel("Frequency", size=font_size)
    # Ticks - change the font size
    ax.tick_params(labelsize = font_size)

    for chain in chains:
        ax.hist(data[chain], bins=25, alpha=0.4, edgecolor='lightgrey')

    # Legend
    ax.legend(chains, frameon=False,
              title = 'Chain id',
              loc='best', prop={'size': font_size-4});
