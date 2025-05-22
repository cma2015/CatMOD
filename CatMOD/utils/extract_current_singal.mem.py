# -*- coding: UTF-8 -*-
# @Author: Minggui Song
# @Date: 2025-02-17
# @Description: This script is designed to extract current features based on genomic positions, 
#               and aggregate the results into files in npy format. 
# @Usage: Run the script with command line arguments to specify input file path, 
#         number of reads to process, number of threads, and output directory.
# @Example: python extract_current_singal.mem.py -f RRACH_motif_linear_Chr1D \
#                  -b  primary_mapped.Chr1D.bed -t 32 -o current_raw_Chr1D_datasets_base
# @Version: v0.0.1a

from argparse import ArgumentParser
from gc import collect
from pathlib import Path

import numpy as np
from rich.progress import Progress

def get_options():
    """Extract options from the command line."""
    parser = ArgumentParser(
        prog='extra_reference.py',
        usage='%(prog)s [options]',
        description='',
    )
    parser.add_argument('-l', '--list', type=str, dest='list',
                        required=True,
                        help='input bed list file path.')
    parser.add_argument('-ls', '--lines', type=int, dest='listlines',
                        required=False, default=0,
                        help='input bed list file lines.')
    parser.add_argument('-b', '--bed', type=str, dest='bed',
                        required=True,
                        help='input bed file path.')
    parser.add_argument('-o', '--output', type=str, dest='output',
                        required=True, help='output dir path.')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s v0.0.1a')
    options = parser.parse_args()
    return options


def read_bed(input_bed: str) -> dict:
    """Read the BED file and parse its contents into a dictionary structure.

    Args:
        input_bed (str): The path of file in BED format.

    Returns:
        dict: {chromosome: {strand: set((start, end), ...)}}
    """
    chr_strand_region_dict = {}
    with open(input_bed, 'r') as open_bed:
        for eachline in open_bed.readlines():
            sp = eachline.strip().split()
            chr_strand_region_dict.setdefault(sp[0], {}).setdefault(sp[5], set()).add((int(sp[1]), int(sp[2])))
    return chr_strand_region_dict


def read_list(input_list: str, chr_strand_region_dict: dict, list_lines: int = 0):
    """Read and process the files in the input list, 
       extract the characteristic information related to chromosomes, chains, and regions, 
       and save the results as.npy files.

    Args:
        input_list (str):  the list of file paths to be processed file paths.
        chr_strand_region_dict (dict): A dictionary containing chromosome, strand, and region information 
                                       to filter for qualifying features.
        list_lines (int, optional): Enter the number of rows in the list file for the progress bar display. Default is 0
    """
    region_info_dict, region_set = {}, set()
    with Progress() as progress:
        if list_lines:
            task = progress.add_task(
                '[green]INFO    [cyan]Reading ONT Current features...',
                total=list_lines)
        else:
            task = progress.add_task(
                '[green]INFO    [cyan]Reading ONT Current features...')
        with open(input_list, 'r') as open_list:
            for eachlist in open_list.readlines():
                with open(eachlist.strip(), 'r') as open_current:
                    for eachline in open_current.readlines():
                        sp = eachline.strip().split()
                        if (int(sp[1])+2, int(sp[2])-2) in chr_strand_region_dict.get(sp[0], {}).get(sp[5], set()):
                            region_string = f'{sp[0]}_{sp[5]}_{int(sp[2])-2}'
                            region_info_dict.setdefault(region_string, {}).setdefault('reads_id', []).append(eachlist.strip().split('/')[-1].split('.')[0])
                            region_info_dict.setdefault(region_string, {}).setdefault('reads_norm_mean', []).append([float(nm) for nm in sp[6].split(',')])
                            region_info_dict.setdefault(region_string, {}).setdefault('reads_norm_stdev', []).append([float(ns) for ns in sp[7].split(',')])
                            region_info_dict.setdefault(region_string, {}).setdefault('reads_current', []).append([float(cs) for bs in sp[8].split(';') for cs in bs.split(',')])
                            region_set.add(region_string)
                progress.advance(task)
    with Progress() as progress:
        task = progress.add_task(
            '[green]INFO    [cyan]Saving ONT Current features...',
            total=len(region_info_dict))
        for region_string in region_set:
            region_info = region_info_dict[region_string]
            np.save(f'{output_dir}/{region_string}.reads_id.npy', region_info['reads_id'])
            np.save(f'{output_dir}/{region_string}.reads_norm_mean.npy', region_info['reads_norm_mean'])
            np.save(f'{output_dir}/{region_string}.reads_norm_stdev.npy', region_info['reads_norm_stdev'])
            np.save(f'{output_dir}/{region_string}.reads_current.npy', region_info['reads_current'])
            # del region_info_dict[region_string]
            # collect()
            progress.advance(task)


def main():
    opt = get_options()
    global output_dir
    output_dir = opt.output
    output_path = Path(output_dir)
    if not output_path.is_dir():
        output_path.mkdir()
    chr_strand_region_dict = read_bed(opt.bed)
    read_list(opt.list, chr_strand_region_dict, opt.listlines)


if __name__ == '__main__':
    main()
