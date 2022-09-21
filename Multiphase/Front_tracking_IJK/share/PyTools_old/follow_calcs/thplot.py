#!/bin/python
import glob, sys, os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# regle un pb de mpl pour enregistrer des figures trop lourdes
matplotlib.rcParams['agg.path.chunksize'] = 100000


# ce script va definir des fonctions a utiliser dans la console et peut s'executer avec son main

def read_temperature_out(head, tdeb=0., tfin=10**6):
    """Read the temperature_*.out files and stock them separately into pandas dataframe.

    Args:
        head: name of the data file
        tdeb: starting time for statistics
        tfin: end time for statistics

    Returns:
        list of pandas DataFrame containing ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u'] for each temperature_*.out
         file.
    """

    fict = "%s_source_temperature*.out" % head
    print("columns : ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u']")
    listt = glob.glob(fict)
    listt.sort()
    print('file list to read : ', listt)

    matt_list = []
    for i, fict in enumerate(listt):
        print(fict)
        matt = pd.read_csv(fict, sep=' ', names=['tstep',  'time', 'theta_adim', 'nu', 'rho_cp_u'])
        matt = matt.query('time > @tdeb')
        matt_list.append(matt)
    return matt_list


def plot_temperature_stat(matt_list, fig_name='DNS', savefig=True, show=True):
    """Creates figures from the temperature stat from ``matt_list``.

    Args:
        matt_list: list of pandas DataFrame containing ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u']
        fig_name: the prefix name for the figures saved
        savefig (bool): flag to save the figures
        show (bool): flag to show figures
    """

    ftemp = []
    for i, matt in enumerate(matt_list):
        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['theta_adim'], '-', label='Theta_adim_moy')
        plt.title('temperature %d' % i)
        plt.legend()

        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['nu'], '-', label='Nu')
        plt.title('temperature %d' % i)
        plt.legend()

        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['rho_cp_u'], '-', label='rho_cp_u')
        plt.title('temperature %d' % i)
        plt.legend()

    if savefig:
        for i, f in enumerate(ftemp):
            f.gca().grid(True)
            f.savefig(fig_name+'_T%d.png' % i)

    if show:
        plt.show()


if __name__ == '__main__':
    # ici on lance les methodes qui lisent les fichiers temperature.out puis on trace les statistiques selon les
    # arguments qui sont passes en entree.
    print(os.getcwd())
    jdd = glob.glob('*.data')[0]
    jdd = jdd.split('.')[0]
    parser = argparse.ArgumentParser(description="Ce programme peut enregistrer des images de l'évolution de la "
                                                 "thermique du cacul")
    parser.add_argument("--data", "-d", help="name for the header of _temperature_*.out", default=jdd)
    parser.add_argument("--fig_name", help="the prefix name of the figures", default="fig")
    parser.add_argument("--no_plot", action="store_true", help="plot the temperature statistics")
    parser.add_argument("--tdeb", "-s", type=float, help="stating time used for the plot", default=0.)
    parser.add_argument("--tfin", "-e", type=float, help="end time used for the plot", default=5.e6)
    args = parser.parse_args()

    print("provided header for name_fig is %s." % args.fig_name)

    matt_list = read_temperature_out(args.data, args.tdeb, args.tfin)
    plot_temperature_stat(matt_list, fig_name=args.fig_name, show=not args.no_plot)

