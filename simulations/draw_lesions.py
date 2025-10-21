import argparse
import os
import logging
import pickle
import numpy as np

from utils import binomial_rv, beta_rv, TX_STRAND_DICT


def read_relavent_driver_info(chr_sizes, cell_dir):
    all_drivers = []
    for chr in range(len(chr_sizes)):
        drivers_file = f'chr{chr}_drivers.pickle'
        drivers = pickle.load(open(f'{cell_dir}/{drivers_file}', 'rb'))
        all_drivers.append(drivers)
    return drivers

def introduce_lesions_chrom_rate(chr_size, lesion_rate):
    chrom_lesion_success = binomial_rv(n_trials=1, prob=lesion_rate, size=chr_size)
    loci = np.where(chrom_lesion_success == 1)[0]
    strand1_successes = binomial_rv(n_trials=1, prob=0.5, size=len(loci))
    strand1_lesions = set(loci[strand1_successes == 1])
    strand2_lesions = set(loci[strand1_successes == 0])
    return strand1_lesions, strand2_lesions

def get_selection_coefficients(driver_s_first, driver_s_max, beta_a, beta_b, n_drivers):
    if driver_s_first is not None:
        s_list = [driver_s_first]
        if n_drivers > 1:
            s_rest = ((1+driver_s_max)/(1+driver_s_first))**(1/(n_drivers-1))-1
            s_list = s_list + [s_rest for d in range(n_drivers-1)]
        s_list = np.array(s_list)
    else:
        s_list = beta_rv(beta_a, beta_b, size=n_drivers)
    return s_list


def introduce_lesions_wrapper(chr_sizes, drivers_loci, lesion_rate, n_drivers, n_other, driver_s_first, driver_s_max, beta_a, beta_b, same_chrom, same_strand): #, strand_type):
    driver_lesions = {}
    other_lesions = {}
    driver_lesions_s = get_selection_coefficients(driver_s_first, driver_s_max, beta_a, beta_b, n_drivers)
    assigned_s = {}
    if lesion_rate is not None:
        raise NotImplementedError
        # for chrom_ind in range(len(chr_sizes)):
        #     chrom = chr_sizes[chrom_ind]
        #     chrom_lesions = introduce_lesions_chrom_rate(chrom, lesion_rate)
        #     homolog_lesions = introduce_lesions_chrom_rate(chrom, lesion_rate)
        #     driver_lesions[chrom_ind] = [chrom_lesions, homolog_lesions]
    else:
        if n_other != 0:
            raise NotImplementedError("Did not implement specific count of non-driver lesions")
            # chrom_other = np.random.randint(len(chr_sizes), size=n_other)
        # pick chr with drivers
        if same_chrom:
            chrom_with_drivers = np.random.randint(len(chr_sizes))
            possible_drivers = sorted(drivers_loci[chrom_with_drivers])
            if len(possible_drivers) < n_drivers:
                raise ValueError("number of required driver lesions is greater than drivers in selected chrom")
            drivers = np.random.choice(possible_drivers, size=n_drivers, replace=False)
            if same_strand:
                driver_lesions[chrom_with_drivers] = [[set(drivers), set()], [set(), set()]]

            else:

                strand1_ind = binomial_rv(n_trials=1, prob=0.5, size=n_drivers)
                strand1_drivers = drivers[strand1_ind == 1]
                strand2_drivers = drivers[strand1_ind == 0]
                driver_lesions[chrom_with_drivers] = [[set(strand1_drivers), set(strand2_drivers)], [set(), set()]]

            for d_ind in range(n_drivers):
                if chrom_with_drivers not in assigned_s:
                    assigned_s[chrom_with_drivers] = {}
                assigned_s[chrom_with_drivers][drivers[d_ind]] = driver_lesions_s[d_ind]
        else:
            chrom_with_drivers, counts = np.unique(np.random.randint(len(chr_sizes), size=n_drivers), return_counts=True)
            driver_s_ind = 0
            for chrom, count in zip(chrom_with_drivers, counts):
                possible_drivers = sorted(drivers_loci[chrom])
                drivers = np.random.choice(possible_drivers, size=count)

                # don't reassign s, since here drivers can be repeated on chorm/strand
                for d_ind in range(count):
                    if chrom not in assigned_s:
                        assigned_s[chrom] = {}
                    if drivers[d_ind] not in assigned_s[chrom]:
                        assigned_s[chrom][drivers[d_ind]] = driver_lesions_s[driver_s_ind]
                    driver_s_ind += 1


                chrom_ind = binomial_rv(n_trials=1, prob=0.5, size=count)
                chrom_drivers = drivers[chrom_ind == 1]
                homolog_drivers = drivers[chrom_ind == 0]

                chrom_strand1_ind = binomial_rv(n_trials=1, prob=0.5, size=len(chrom_drivers))
                chrom_1 = chrom_drivers[chrom_strand1_ind == 1]
                chrom_2 = chrom_drivers[chrom_strand1_ind == 0]

                hom_strand1_ind = binomial_rv(n_trials=1, prob=0.5, size=len(homolog_drivers))
                homolog_1 = homolog_drivers[hom_strand1_ind == 1]
                homolog_2 = homolog_drivers[hom_strand1_ind == 0]

                driver_lesions[chrom] = [[set(chrom_1), set(chrom_2)], [set(homolog_1), set(homolog_2)]]
    for chrom_ind in range(len(chr_sizes)):
        if chrom_ind not in driver_lesions:
            driver_lesions[chrom_ind] = [[set(), set()], [set(), set()]]
            assigned_s[chrom_ind] = {}

    for k, v in assigned_s.items():
         logging.info(f'{k}: {v}')
    return assigned_s, driver_lesions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir', '-wd', required=True, help='path to base directory')
    parser.add_argument('--cell_dir', '-c', required=True, help='path to file with cell driver info to use')
    parser.add_argument('--lesion_rate', '-lr', required=False, help='rate of lesion introduction, must be set if lesion counts are not set.')
    parser.add_argument('--lesion_counts', '-lc', required=False, help='tuple of # drivers, # other lesions to draw. Must be set if lesion rate not set.')
    parser.add_argument('--drivers_same_chrom', '-dc', required=False, default='False', help='whether drivers should be on same chrom, only relevant if using counts')
    parser.add_argument('--drivers_same_strand', '-ds', required=False, default='False',
                        help='whether drivers should be on same strand, only relevant if using counts and drivers same chrom set to True')
    parser.add_argument('--driver_s_first', '-s', required=False, default=None, help='selection coef to apply to first driver. Only use if drifer count specified'
                                                                                     'If multiple drivers: smax must also be set, '
                                                                                     'and the remaining drivers will get equal s that combine to max s')
    parser.add_argument('--driver_s_max', '-smax', required=False, default=1, help='max selection coef if all drivers acquired, only set if driver count specified as > 1')
    parser.add_argument('--s_beta_param_a', '-a', required=False, default=1,
                        help='value to use for "a" parameter of beta distribution for selection coefficients if s not set')
    parser.add_argument('--s_beta_param_b', '-b', required=False, default=100,
                        help='value to use for "b" parameter of beta distribution for selection coefficients if s not set')
    parser.add_argument('--seed', '-rs', required=False, default=7, help='random seed to use')
    parser.add_argument('--outsuffix', '-o', required=False, default='')
    args = parser.parse_args()

    lesion_rate = None if args.lesion_rate is None else float(args.lesion_rate)
    if lesion_rate is None:
        if args.lesion_counts is None:
            raise ValueError("either lesion rate of lesion counts must be set")
        else:
            lesion_counts = args.lesion_counts.split(',')
            n_drivers = int(lesion_counts[0])
            n_other = int(lesion_counts[1])
            lesions_str = f'{n_drivers}_drivers_{n_other}_other'

    else:
        n_drivers = None
        n_other = None
        lesions_str = f'rate_{lesion_rate}'

    same_chrom = True if args.drivers_same_chrom == 'True' else False
    same_strand = True if args.drivers_same_strand == 'True' else False

    if same_chrom:
        lesions_str = f'{lesions_str}_same_chrom'
        if same_strand:
            lesions_str = f'{lesions_str}_strand'


    if args.driver_s_first is not None:
        if n_drivers is None:
            raise ValueError("driver s can only be specified with specific driver count")
        driver_s_first = float(args.driver_s_first)
        # driver_s_max = None if args.driver_s_max is None else float(args.driver_s_max)
        driver_s_max = float(args.driver_s_max)
        beta_a = None
        beta_b = None
    else:
        ## use beta
        beta_a = float(args.s_beta_param_a)
        beta_b = float(args.s_beta_param_b)
        driver_s_first = None
        driver_s_max = None

    np.random.seed(int(args.seed))

    cell_dir = f'{args.working_dir}/{args.cell_dir}/'
    outdir = f'{cell_dir}/lesions_{lesions_str}'
    if args.outsuffix:
        outdir = f'{outdir}_{args.outsuffix}'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_handler = logging.FileHandler(f'{outdir}/output.log', 'w')
    logger.addHandler(logger_handler)

    logger.info(f'Drawing lesions for cell in {args.cell_dir}')
    logger.info(f'Lesion rate is: {lesion_rate}')
    logger.info(f'Lesion counts are: drivers = {n_drivers}, other = {n_other}')
    logger.info(f'Driver same chromosome: {same_chrom}, same strand: {same_strand}')
    # logger.info(f'Driver lesions being drawn on {args.strand_type} strands')
    if beta_a is None:
        if n_drivers == 1:
            logger.info(f'Driver lesions are selected to have selection coefecients {driver_s_first}')
        else:
            logger.info(f'Driver lesions are selected to have max selected selection coefecients {driver_s_max}, '
                        f'with one driver having selection coeficient {driver_s_first}')
    else:
        logger.info(f'Driver lesions are have selection coeficients follwoing beta distribution with params ({beta_a}, {beta_b})')
    logger.info(f'Random seed {args.seed}')

    # read cell info
    chr_sizes = np.loadtxt(f'{cell_dir}/chr_sizes.txt', dtype=int)
    driver_locations = read_relavent_driver_info(chr_sizes, cell_dir)

    # draw lesions
    assigned_s, lesions = introduce_lesions_wrapper(chr_sizes, driver_locations, lesion_rate, n_drivers, n_other, driver_s_first, driver_s_max, beta_a, beta_b, same_chrom, same_strand) #, args.strand_type)
    pickle.dump(lesions, open(f'{outdir}/driver_lesions.pickle', 'wb'))
    pickle.dump(assigned_s, open(f'{outdir}/driver_lesion_selection_coefficients.pickle', 'wb'))


if __name__ == '__main__':
    main()