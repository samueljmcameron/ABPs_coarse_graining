import matplotlib.pyplot as plt

def JupyterPlots():

    import matplotlib.pyplot as plt
    import numpy as np
    
    gold = (np.sqrt(5)-1)/2
    figsize=[3.73,3.73*gold]
    params = {'backend': 'pdf',
              'axes.labelsize': 14,
              'font.size': 14,
              'legend.fontsize': 10,
              'xtick.labelsize': 11,
              'ytick.labelsize': 11,
              'text.usetex': True,
              'figure.figsize': figsize,
              'figure.dpi': 200,
              'text.latex.preamble': [r"\usepackage{amstext}"]}
    plt.rcParams.update(params)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{bm}')

    return figsize


