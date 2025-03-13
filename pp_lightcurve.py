#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
import argparse
import warnings
import os
from pathlib import Path
warnings.filterwarnings("ignore")

def lister(path: Path, name_pattern: str, return_type='name', object_type="file") -> list:
    """
    Lists all objects that match the naming pattern by the given path
    :param path: path to the directory
    :param name_pattern: name and format of the files: (e.g. "sim_*.dat")
    :param return_type: "name" - returns only the names of the files in the directory
                        "path" - returns full paths to the files in the directory
    :return: list of paths (or names) of the files that satisfy given conditions
    """
    # get paths to every object with specified name pattern
    objects = sorted(list(Path(path).glob('{}'.format(name_pattern))))
    # check if the object is a file or a directory
    if object_type == "file":
        objects = [obj_path for obj_path in objects if os.path.isfile(obj_path)]
    elif object_type == "dir":
        objects = [obj_path for obj_path in objects if os.path.isdir(obj_path)]
    # get only the names of the objects
    if return_type == 'name':
        objects = sorted([obj.name for obj in objects])
    return objects


def lightcurve_plots(data, flagged_data=None, target=None, path: str=None):
    """build lightcurve plots for targets"""
    fig, ax = plt.subplots()
    ax.set_title(target)
    ax.set_xlabel('Minutes after {:s} UT'.format(Time(data['julian_date'].min(), format='jd').to_value('iso', subfmt='date_hm')))
    ax.set_ylabel('Magnitude')
    ax.errorbar(data['time_from_start']*1440,
                data['mag'],
                yerr=data['sig'],
                linestyle='', color='red',
                marker='o', capsize=3)
    # plot flagged data
    if flagged_data is not  None:
        ax.errorbar(flagged_data['time_from_start']*1440,
                flagged_data['mag'],
                yerr=flagged_data['sig'],
                linestyle='', color='orange',
                marker='o', capsize=3)
    ax.set_ylim([ax.get_ylim()[1], ax.get_ylim()[0]])
    ax.set_xticklabels = [Time(t, format='jd').to_value('iso',
                                                        subfmt='date_hm')
                          for t in plt.xticks()[0]]
    ax.grid()
    fig.savefig(path, format='png', dpi=200)


def plot_lightcurve(data_path, plot_flagged=True, target_name: str=None, save_name: str='lightcurve.csv'):
    """"""""
    # check if file is csv
    if str(data_path).endswith('.csv'):
        # read csv file with lightcurve data
        df = pd.read_csv(data_path)
        zerotime = df['julian_date'].min()
        df['time_from_start'] = df['julian_date'] - zerotime
        good_data = df[(~df['rejected']) & (df['sextractor_flags'] == 0)]
        flagged_data = df[(~df['rejected']) & (df['sextractor_flags'] > 1)]

    # plot data
    if not plot_flagged:
        flagged_data = None
    lightcurve_plots(data=good_data, flagged_data=flagged_data, target=target_name, path=save_name)
    print('Done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='creates the plot of the lightcurve')
    parser.add_argument('-file_path', help='which file to use for photometry data',
                        default=None)
    parser.add_argument('-plot_flagged', help='select if to plot data flagged by sextractor',
                        default=False, action='store_true')
    parser.add_argument('-target_name', help='name of the object to use as a header',
                        default=None)
    parser.add_argument('-save_name', help='name to use for the resulting file',
                        default=None)

    args = parser.parse_args()


    file_path = args.file_path
    plot_flagged = args.plot_flagged
    target_name = args.target_name
    save_name = args.save_name
    print(save_name)

    if file_path is None:
        file_path = lister(os.getcwd(), '*.csv', 'path', 'file')[0]

    if save_name is None:
        save_name = 'lightcurve.png'
    if target_name is not None:
        save_name = f'{target_name}_lightcurve.png'

    plot_lightcurve(data_path=file_path, plot_flagged=plot_flagged, target_name=target_name, save_name=save_name)


