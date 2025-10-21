import argparse
import copy
import os
import logging
import datetime
import pickle

import numpy as np

from utils import beta_rv, binomial_rv, normal_rv

def get_driver_loci(chr_size, freq, count=None):
    if count is None:
        successes = binomial_rv(n_trials=1, prob=freq, size=chr_size)
        loci = np.where(successes==1)[0]
    else:
        loci = np.random.randint(low=0, high=chr_size, size=count)
    return loci

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir', '-wd', required=True, help='path to base directory')
    parser.add_argument('--chr_file', '-c', required=True, help='path to file with chromosomes sizes')
    parser.add_argument('--driver_freq', '-d', required=False, default=5, help='frequency of driver loci or count per chromosome')
    parser.add_argument('--outprefix', '-o', required=False, help='string to add to outdir')
    args = parser.parse_args()

    chr_sizes = np.loadtxt(f'{args.working_dir}/{args.chr_file}', dtype=int)
    chr_file_str = '.'.join(args.chr_file.split('/')[-1].split('.')[:-1])
    cells_dir = f'{args.working_dir}/cells/'
    if not os.path.isdir(cells_dir):
        os.mkdir(cells_dir)
    cell_dir = f'{cells_dir}/{chr_file_str}_driver_loci_freq_{args.driver_freq}/'
    if not os.path.isdir(cell_dir):
        os.mkdir(cell_dir)

    day = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")

    if args.outprefix != None:
        cell_dir = f'{cell_dir}/{day}_{args.outprefix}/'
    else:
        cell_dir = f'{cell_dir}/{day}/'
    if not os.path.isdir(cell_dir):
        os.mkdir(cell_dir)
    seed = 42
    np.random.seed(42)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_handler = logging.FileHandler(f'{cell_dir}/output.log', 'w')
    logger.addHandler(logger_handler)
    logger.info(f'Creating cell with chrom sizes from file {args.chr_file}, drivers drawn as with frequency {args.driver_freq}')

    logger.info(f'Random seed {seed}')

    driver_freq = float(args.driver_freq)
    driver_count = int(driver_freq) if driver_freq > 1 else None
    drivers = {}
    for chr in range(len(chr_sizes)):
        loci = get_driver_loci(chr_sizes[chr], driver_freq, driver_count)
        drivers[chr] = loci
        logging.info(f'drivers in chr {chr}: {loci}')
        outprefix = f'{cell_dir}/chr{chr}_'
        pickle.dump(drivers, open(f'{outprefix}drivers.pickle', 'wb'))

    np.savetxt(f'{cell_dir}chr_sizes.txt', chr_sizes, fmt='%d')






if __name__ == '__main__':
    main()