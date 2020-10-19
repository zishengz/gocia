import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
matplotlib.use('Agg')

def heatmap(matrix, baseName):
    dimension = len(matrix)
    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap=cm.Blues)
    fig.colorbar(im)
    ax.tick_params(top=True, bottom=False,
                    labeltop=True, labelbottom=False)
    ax.set_xticks(np.arange(dimension))
    ax.set_yticks(np.arange(dimension))
    # ... and label them with the respective list entries
    ax.set_xticklabels([str(i) for i in np.arange(dimension)], fontsize='8')
    ax.set_yticklabels([str(i) for i in np.arange(dimension)], fontsize='8')
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right",
            rotation_mode="anchor")
    fig.tight_layout()
    plt.savefig(
        '%s_heat.png'%(baseName),
        dpi=200
        )