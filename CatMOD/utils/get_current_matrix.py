from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from argparse import ArgumentParser
from rich.progress import Progress

def get_options():
    """Extract options from the command line."""
    parser = ArgumentParser(
        prog='extra_sequence_matrix.py',
        usage='%(prog)s [options]',
        description='',
    )
    parser.add_argument('-f', '--files', type=str, dest='files',
                       required=False, help='input file list.')
    parser.add_argument('-d', '--dir', type=str, dest='dir',
                       required=False, help='input fast5 dir path.')
    parser.add_argument('-t', '--threads', type=int, dest='threads',
                        default=24, required=False,
                        help='number of threads to use. default [all]')
    parser.add_argument('-o', '--output', type=str, dest='output',
                        required=True, help='output dir path.')
    parser.add_argument('-n', '--name', type=str, dest='name',
                        required=True, help='output dir path.')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s v0.0.1a')
    options = parser.parse_args()
    return options

def process_line(eachline, current_raw_dir):
    sp = eachline.strip().split()
    chrom, strand, region = sp[0], sp[5], int(sp[2])
    region_string = f'{chrom}_{strand}_{region}'

    features_path = Path(f'{current_raw_dir}/current_raw_{chrom}_datasets_base/{region_string}.reads_id.npy')
    if features_path.is_file():
        reads_norm_mean = np.load(f'{current_raw_dir}/current_raw_{chrom}_datasets_base/{region_string}.reads_norm_mean.npy').astype(np.float32)
        reads_norm_stdev = np.load(f'{current_raw_dir}/current_raw_{chrom}_datasets_base/{region_string}.reads_norm_stdev.npy').astype(np.float32)
        reads_current = np.load(f'{current_raw_dir}/current_raw_{chrom}_datasets_base/{region_string}.reads_current.npy').astype(np.float32)

        if reads_norm_mean.shape[0] == reads_norm_stdev.shape[0] == reads_current.shape[0] > 0 and \
           reads_norm_mean.shape[1] == reads_norm_stdev.shape[1] == 5 and reads_current.shape[1] == 300:
            feature_vector = np.concatenate((reads_norm_mean.mean(0), reads_norm_stdev.mean(0), reads_current.mean(0)))
            return region_string, feature_vector
    return None, None

def extract_current_matrix(bed_file, current_raw_dir, output_dir, num_threads=4):
    current_dict = {'features': [], 'feature_names': []}

    with open(bed_file, 'r') as open_bed:
        lines = open_bed.readlines()

        with Progress() as progress:
            task = progress.add_task("[cyan]Processing BED file...", total=len(lines))

            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                futures = {executor.submit(process_line, eachline, current_raw_dir): eachline for eachline in lines}

                for future in as_completed(futures):
                    feature_name, result = future.result()
                    if result is not None:
                        current_dict['features'].append(result)
                        current_dict['feature_names'].append(feature_name)
                    progress.update(task, advance=1)

    np.save(f'{output_dir}/{matrix_name}.npy', np.array(current_dict['features'], dtype=np.float32))
    np.save(f'{output_dir}/{matrix_name}_names.npy', np.array(current_dict['feature_names'], dtype=np.object_))

if __name__ == '__main__':
    # iwgsc2.1.primary_mapped.all_sample2_depth20_no_overlap_with_Nanom6A_m6Anet_and_ngs.bed
    # m6Anet_Nanom6A_ngs_overlap_sites.bed
    opt = get_options()

    bed_file = opt.files
    threads = opt.threads
    current_raw_dir = '/home/smg/WorkDir/4.wheat/0.data/PAG40520/8.datasets'
    output_dir = opt.output
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    global matrix_name
    matrix_name = opt.name
    
    extract_current_matrix(bed_file, current_raw_dir, output_dir, num_threads=8)