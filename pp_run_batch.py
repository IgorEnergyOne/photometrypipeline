#!/usr/bin/env python3

import os
import argparse
import subprocess
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from toolbox import lister


import sys

def pipeline_batch(fname: str, skip_pipeline: bool = False, by_filter: bool = False):
    """
    wrapper around pp_run for processing of multiple directories in single run
    :param fname: file name with the list of commands to be executed
    :param skip_pipeline: skip pipeline processing and perform only postprocessing
    :param by_filter: create separate csv and atlas files for each filter present
    :return:
    """
    has_error = False
    with open(fname) as f:
        lines = f.readlines()
        # remove empty lines
        lines = [line for line in lines if line.strip()]
    try:
        paths = []
        if skip_pipeline:
            print("Skipping pipeline processing")
        for line in tqdm(lines, desc='Batch processing', unit='dirs'):
            # parse the command to get the path to directory
            path = line.rsplit(" ")[-1].strip('\n')

            # check if the path is absolute or not
            if Path(path).is_absolute():
                dir_path = path
            else:
                dir_path = Path(core_path) / path
            # add path to the list
            paths.append(dir_path)
            # change cwd to dir path
            os.chdir(dir_path)
            if skip_pipeline:
                continue
            # check if cmd is marked as skip with '!': means skip pipeline processing for the directory
            elif line.startswith('!') and not skip_pipeline:
                print(f'Skipping {dir_path} for pipeline processing')
            else:
                # compile the command to run pipeline
                command = " ".join(line.rsplit(' ')[:-1]) + ' *.fit*'
                try:
                    # check_call raises CalledProcessError on non-zero exit code
                    subprocess.check_call(['/bin/sh', '-i', '-c', command])
                except Exception as e:
                    print(f'Error in {dir_path}, exception: {e}')
                    has_error = True
                    continue
        # change cwd one level up to the base directory
        prnt = Path(paths[0]).parent
        os.chdir(prnt)
        # if only one directory is processed, add empty string to the list
        if len(paths) == 1:
            paths = [paths[0], '']
        # try to combine data to one csv and atlas, build photometry curve
        try:
            subprocess.check_call(
                ['/bin/sh', '-i', '-c', f'pp_combine_csv -dirs_pattern {",".join(str(path) for path in paths)}'])
            atlas_cmd = f"pp_atlas -combine {' '.join(str(p) for p in paths)}"
            subprocess.check_call(['/bin/sh', '-i', '-c', atlas_cmd])
        except Exception as e:
            print(f"Error in post-processing combination: {e}")
            has_error = True

        if by_filter:
            # group directories by photometry filters
            filter_groups = {}
            for path in paths:
                if len(str(path)) == 0:
                    continue
                # get control star_data for each directory
                photo_file = lister(path, name_pattern="photometry_*_.csv",
                                    object_type='file', return_type='path')[0]
                if photo_file.exists():
                    photo_data = pd.read_csv(photo_file)
                    if not photo_data.empty and 'band' in photo_data.columns:
                        # Get the unique filter for this directory
                        filter_name = photo_data['band'].iloc[0]
                        if filter_name not in filter_groups:
                            filter_groups[filter_name] = []
                        filter_groups[filter_name].append(path)
            # Process each filter group separately
            for filter_name, filter_paths in filter_groups.items():
                print(f'\nProcessing filter: {filter_name}')
                if len(filter_paths) == 1:
                    filter_paths = [filter_paths[0], '']

                try:
                    csv_filename = f'combined_results_{filter_name}.csv'
                    # Combine CSV for this filter
                    subprocess.check_call(['/bin/sh', '-i', '-c',
                                     f'pp_combine_csv -dirs_pattern {",".join(str(p) for p in filter_paths)} -out_path {csv_filename}'])

                    # Create atlas file for this filter
                    if Path(csv_filename).exists():
                        atlas_cmd = f"pp_atlas -combine {' '.join(str(p) for p in filter_paths)} -fname_out combined_atlas_{filter_name}.ATL"
                        subprocess.check_call(['/bin/sh', '-i', '-c', atlas_cmd])
                        try:
                            target_name = photo_data["target"].iloc[0].replace(' ', '_')
                            # strip name of any special characters
                            target_name = target_name.replace('(', '').replace(')', '')
                        except KeyError:
                            target_name = f'Asteroid'
                        lightcurve_cmd = f'pp_lightcurve ' \
                                         f'-file_path {csv_filename} ' \
                                         f'-target_name "{target_name} ({filter_name} filter)" ' \
                                         f'-save_name lightcurve_{target_name}_{filter_name}.png ' \
                                         f'-plot_flagged'
                        subprocess.check_call(['/bin/sh', '-i', '-c', lightcurve_cmd])
                except Exception as e:
                    print(f"Error processing filter {filter_name}: {e}")
                    has_error = True

    except Exception as e:
        print(f'Error in {dir_path}, exception: {e}')
        if e is KeyboardInterrupt:
            print("The script was interrupted by the user. Aborting...")
        has_error = True
    
    if has_error:
        print("Batch processing completed with errors.")
        sys.exit(1)


def meta_batch(fname: str):
    """
    Process a list of batch queues.
    :param fname: path to the meta-queue file containing list of queue files and optional args
    """
    with open(fname) as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    # remove lines starting with '!' ( means skip pipeline processing for the directory)
    lines = [line for line in lines if not line.startswith('!')]


    print(f"Starting Meta-Batch processing with {len(lines)} queues.")

    successful_batches = 0
    failed_batches = []

    for line in tqdm(lines, desc='Meta Batch', unit='queue'):
        # format: /path/to/queue.txt [args...]
        parts = line.split()
        if not parts:
            continue
        
        queue_file = parts[0]
        extra_args = parts[1:]
        
        # change directory to where the queue file resides.
        
        qpath = Path(queue_file)
        if not qpath.is_absolute():
            # If relative path, assume it is CWD from queue batch file
            qpath = Path(os.getcwd()) / qpath

        # If the file exists, we use its parent as the working directory
        if not qpath.exists():
            print(f"Warning: Queue file {qpath} not found. Skipping.")
            failed_batches.append(f"{qpath} (File not found)")
            continue

        work_dir = qpath.parent
        queue_filename_only = qpath.name
        
        cmd = ['pp_run_batch', '-file', queue_filename_only] + extra_args
        
        print(f"\n>>> Launching batch: {queue_filename_only}")
        print(f"    In Directory: {work_dir}")
        print(f"    Command: {' '.join(cmd)}")
        
        try:
            # Run subprocess in the separate directory
            # check=True will raise CalledProcessError if return code != 0
            subprocess.check_call([str(c) for c in cmd], cwd=work_dir)
            successful_batches += 1
        except Exception as e:
            print(f"Error executing batch {qpath}: {e}")
            failed_batches.append(f"{qpath} ({str(e)})")
            if isinstance(e, KeyboardInterrupt):
                raise e

    print("\n" + "="*40)
    print("Meta-Batch Processing Summary")
    print("="*40)
    print(f"Total Batches: {len(lines)}")
    print(f"Successful:    {successful_batches}")
    print(f"Failed:        {len(failed_batches)}")
    
    if failed_batches:
        print("\nFailed Batches:")
        for fb in failed_batches:
            print(f"  - {fb}")
    print("="*40 + "\n")



if __name__ == '__main__':
    core_path = os.getcwd()
    parser = argparse.ArgumentParser(description='photometrypipeline batch processing')
    parser.add_argument('-file',
                        help='name of the file with the list of commands to be executed',
                        default="queue.txt")
    parser.add_argument('-skip_pipeline',
                        help='skip pipeline processing',
                        default=False, action='store_true')
    parser.add_argument('-by_filter',
                        help='perform creation of atlas and lightcurve files for each filter',
                        default=False, action='store_true')
    parser.add_argument('-meta',
                        help='run in meta-batch mode (file argument is a list of queues)',
                        default=False, action='store_true')
    args = parser.parse_args()
    filename = str(args.file)
    skip_pipeline = args.skip_pipeline
    by_filter = args.by_filter
    is_meta = args.meta

    if is_meta:
        meta_batch(filename)
    else:
        pipeline_batch(filename, skip_pipeline, by_filter)