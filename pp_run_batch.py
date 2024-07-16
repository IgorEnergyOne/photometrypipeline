#!/usr/bin/env python3

import os
import argparse
import subprocess
from pathlib import Path
from tqdm import tqdm

def pipeline_batch(fname: str):
    """
    wrapper around pp_run for processing of multiple directories in single run
    :param fname: file name with the list of commands to be executed
    :return:
    """
    with open(fname) as f:
        lines = f.readlines()
    try:
        for line in tqdm(lines, desc='Batch processing', unit='dirs'):
            # parse the command to get the path to directory
            path = line.rsplit(" ")[-1].strip('\n')

            # check if the path is absolute or not
            if Path(path).is_absolute():
                dir_path = path
            else:
                dir_path = Path(core_path) / path
            # change cwd to dir path
            os.chdir(dir_path)
            command = " ".join(line.rsplit(' ')[:-1]) + ' *.fit*'
            try:
                subprocess.call(['/bin/sh', '-i', '-c', command])
            except Exception as e:
                print(f'Error in {dir_path}, exception: {e}')
                continue
    except KeyboardInterrupt:
        print("The script was interrupted by the user. Aborting...")


if __name__ == '__main__':
    core_path = os.getcwd()
    parser = argparse.ArgumentParser(description='photometrypipeline batch processing')
    parser.add_argument('-file',
                        help='name of the file with the list of commands to be executed',
                        default="queue.txt")
    args = parser.parse_args()
    filename = str(args.file)

    pipeline_batch(filename)

