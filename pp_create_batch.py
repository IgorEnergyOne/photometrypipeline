#!/usr/bin/env python3

import shutil
import os
import argparse
from pathlib import Path


def lister(path: Path, name_pattern: str, return_type='name') -> list:
    """
    Lists all objects that match the naming pattern by the given path
    :param path: path to the directory
    :param name_pattern: name and format of the files: (e.g. "sim_*.dat")
    :param return_type: "name" - returns only the names of the files in the directory
                        "path" - returns full paths to the files in the directory
    :return: list of paths (or names) of the files that satisfy given conditions
    """
    if return_type == 'name':
        objects = list(Path(path).glob('{}'.format(name_pattern)))
        objects = sorted([obj.name for obj in objects])
    elif return_type == 'path':
        objects = list(Path(path).glob('{}'.format(name_pattern)))
    else:
        raise Exception('Specify the return type ("name" or "path")')
    return objects


def create_batches(pp_command: str, base_dir: Path, subdir_name: str, batch_fname: str='queue.txt',
                   batch_size: int = 3, overlap: int = 1):
    """
    Splits the observation series into separate batches of a given size, and creates a command file for each batch
    to be run with pp_run_batch command
    :param pp_command: pp_run command with parameters to be executed
    :param base_dir: base directory which contains the observation files
    :param subdir_name: name of the subdirectory to store batches
    :param batch_fname: name of the file with the list of commands for pp_run_batch
    :param batch_size: number of files in each batch
    :param overlap: number of files to overlap between batches
    """
    # create a list of commands for each batch for a text file
    batch_cmd = []
    fpaths = lister(base_dir,
                    name_pattern='*.fit*', return_type='path')
    fpaths = sorted(fpaths)
    base_subdir = base_dir / subdir_name
    # create a subdirectory for each batch and copy the files there
    os.mkdir(base_subdir)
    # get the base directory
    for i, file_idx in enumerate(range(0, len(fpaths), batch_size - overlap)):
        # path for every file in the batch
        batch = fpaths[file_idx: file_idx + batch_size]
        # create a subdirectory for each batch
        batch_dir_name = f"batch_{i+1:02d}"
        batch_subdir = base_subdir / batch_dir_name
        os.mkdir(batch_subdir)
        line = pp_command + " " + batch_dir_name + "\n"
        batch_cmd.append(line)
        # copy the files to the subdirectory
        [shutil.copy(file_path, batch_subdir / Path(file_path).name) for file_path in batch]

    with open(base_subdir / batch_fname, 'w') as f:
        f.writelines(batch_cmd)

    print('Done!')

if __name__ == '__main__':
    core_path = os.getcwd()
    parser = argparse.ArgumentParser(description='splitting observation into batches for batch processing')
    parser.add_argument('-pp_command',
                        help='pp_run command to be executed for each batch (e.g. pp_run -auto)',
                        default="pp_run -auto")
    parser.add_argument('-base_dir',
                        help='directory with the fits files to be split',
                        default=f"{core_path}")
    parser.add_argument('-subdir_name',
                        help='name of the directory to contain batches',
                        default="batches")
    parser.add_argument('-cmd_fname',
                        help='name of the file with the list of commands to be executed',
                        default="queue.txt")
    parser.add_argument('-batch_size',
                        help='number of files in each batch',
                        default=3)
    parser.add_argument('-overlap',
                        help='number of files to overlap between batches',
                        default=1)

    args = parser.parse_args()

    pp_command = str(args.pp_command)
    base_dir = Path(str(args.base_dir))
    subdir_name = str(args.subdir_name)
    cmd_filename = str(args.cmd_fname)
    batch_size = int(args.batch_size)
    overlap = int(args.overlap)

    create_batches(pp_command=pp_command, base_dir=base_dir, subdir_name=subdir_name,
                   batch_fname=cmd_filename, batch_size=batch_size, overlap=overlap)