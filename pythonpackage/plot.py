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

def plot_msd(data, percentage=1, legend=None):
    """
    Plot the mean square displacement

    :param data: list of msd
    :type data: pd.Series
    :param percentage: percentage of the data that will be plotted,
                       value between 0 and 1.
    :type percentage: float
    :param legend: data legend
    :type legend: list of string
    """
    fig, ax = plt.subplots(figsize=(10,7))

    xmax, ymax = 0, 0
    # plot data
    for msd in data:
        ax.plot(msd.values)
        ax.set_xlabel('#amino acids', size=16)
        ax.set_ylabel(r'$MSD [nm^2]$', size=16)
        if int(len(msd)*percentage)-1 > xmax:
            xmax = int(len(msd)*percentage)-1
        if msd.iloc[int(len(msd)*percentage)-1] > ymax:
            ymax = msd.iloc[int(len(msd)*percentage)-1]
        ax.set_xlim([0, xmax])
        ax.set_ylim([0, ymax])

    # labels
    ax.tick_params(labelsize=16)

    # Legend
    if legend:
        ax.legend(legend,
                  loc='best',
                  frameon=True,
                  shadow=True,
                  facecolor='#FFFFFF',
                  framealpha=0.9,
                  fontsize=14)

    # Remove spines and ticks
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none');
