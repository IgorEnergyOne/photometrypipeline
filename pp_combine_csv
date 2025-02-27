#!/usr/bin/env python3
# pp_combine_csv.py

import os
import argparse
import pandas as pd
from pathlib import Path


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


def combine_csv_files(path_core_dir: str,
                      dir_name_pattern: str,
                      file_name_pattern: str,
                      out_file_path: str,
                      append=False):
    """
    Combines all csv files in the given directory with the given name pattern
    """
    # get paths to every directory with specified name pattern
    dirs = lister(Path(path_core_dir), name_pattern=dir_name_pattern, object_type="dir", return_type="path")
    # get csv files in every directory
    target_data = []
    control_data = []
    for dir in dirs:
        # get control star_data for each directory
        control_star_path = lister(dir, name_pattern="photometry_Control_Star.csv",
                                   object_type='file', return_type='path')
        control_data.append(pd.read_csv(control_star_path[0]))
        # get target data for each directory
        file = lister(dir, name_pattern=file_name_pattern, object_type='file', return_type='path')
        target_data.append(pd.read_csv(file[0]))
    # combine data
    control_data = pd.concat([*control_data], ignore_index=True)
    control_data.rename(columns={"mag": "mag_control", "sig": "sig_control",
                                 "source_ra": "control_ra", "source_dec": "control_dec",
                                 "inst_mag": "control_mag_inst", "inst_sig": "control_sig_inst"}, inplace=True)
    target_data = pd.concat([*target_data], ignore_index=True)
    combined_data = target_data.join(control_data[["mag_control", "sig_control", "control_ra", "control_dec",
                                                   "control_mag_inst", "control_sig_inst"]],
                                     how='left', rsuffix='_control')
    combined_data['reduc_julian_date'] = combined_data['julian_date'] - combined_data.iloc[0]['julian_date']
    # calculate relative instrumental magnitude (source - control star)
    combined_data['rel_mag'] = combined_data['inst_mag'] - combined_data['control_mag_inst']
    # calculate instrumental magnitude error
    combined_data['rel_sig'] = (combined_data['inst_sig']**2 + combined_data['control_sig_inst']**2)**0.5
    if append:
        combined_data.to_csv(out_file_path, mode='a', header=False, index=False)
    else:
        combined_data.to_csv(out_file_path, index=False)


if __name__ == '__main__':
    # command line arguments
    parser = argparse.ArgumentParser(description='combine results of photometry from different directories into one file')
    parser.add_argument('-path_core_dir',
                        help='path to directory containing directories with partial results', default=None)
    parser.add_argument('-dirs_pattern',
                        help='name pattern for directories',
                        default="*")
    parser.add_argument('-files_pattern',
                        help='name pattern for files containing data',
                        default="photometry*_.csv")
    parser.add_argument('-out_path',
                        help='path where to output the combined results',
                        default=None)
    parser.add_argument('-append',
                        help='append data to the existing csv file',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    core_path = str(args.path_core_dir) if args.path_core_dir is not None else os.getcwd()
    dir_name_pattern = str(args.dirs_pattern)
    file_name_pattern = str(args.files_pattern)
    out_path = args.out_path
    append = args.append
    if out_path is None:
        out_path = os.path.join(core_path, 'combined_results.csv')

    if append:
        print('Appending data to {}'.format(out_path))
        combine_csv_files(core_path, dir_name_pattern, file_name_pattern, out_path, append=True)
    else:
        combine_csv_files(core_path, dir_name_pattern, file_name_pattern, out_path)
        print('Done! Results are saved in {}'.format(out_path))

