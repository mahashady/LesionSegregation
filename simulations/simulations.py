import os
import datetime
import argparse
import logging
import pickle

import yaml

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from simulation_objects import Cell

CLONE_MIN_PROPORTION = 0.1


def update_cell_information(cells_info, new_counts):
    total = sum(new_counts)
    updated_cell_info = []
    for clone in range(len(cells_info)):
        updated_cell_info.append([cells_info[clone][0], cells_info[clone][1],
                                  cells_info[clone][2], cells_info[clone][3],
                                  new_counts[clone]/total, new_counts[clone]])
    return updated_cell_info

def group_clones_by_drivers(clones_by_id, selection_coefs):
    '''
    :param clones_by_id:
    :return: list of grouped clones,
                each group is a list with
                    first item = list of ids
                    second item = driver info
                    third item = summed proportion
    '''
    groups = []
    groups.append([[clones_by_id[0][1]], clones_by_id[0][0].get_driver_info(selection_coefs), clones_by_id[0][-2]])
    for clone_ind in range(1, len(clones_by_id)):
        clone = clones_by_id[clone_ind]
        clone_drivers = clone[0].get_driver_info(selection_coefs)
        for group_ind in range(len(groups)):
            if clone_drivers == groups[group_ind][1]:
                groups[group_ind][0].append(clone[1])
                groups[group_ind][-1] += clone[-2]
            else:
                groups.append([[clone[1]], clone_drivers, clone[-2]])
    return groups

def get_average_fitness(cells_info):
    out = 0
    for cell in cells_info:
        out += (cell[2]*cell[4])
    return out

def get_cell_propagation_probability(cells_info, w_bar):
    cell_prob = []
    for cell in cells_info:
        prob = cell[4]*(cell[2]/w_bar)
        cell_prob.append(prob)
    return cell_prob


def driver_lesions_remaining(cells, selection_coefs):
    for cell in cells:
        if cell.has_driver_lesions(selection_coefs):
            return True
    return False

def lesions_remaining(cells):
    for cell in cells:
        if cell.has_lesions():
            return True
    return False

def log_cell_info(cell, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, smax):
    logging.info(f'Cell {cell.id}')
    logging.info(f'Driver mutations: {cell.get_driver_mutations(selection_coefs)}')
    logging.info(f'Driver lesion status: {cell.has_driver_lesions(selection_coefs)}')
    logging.info(f'Cell fitness: {cell.get_fitness(selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, smax)}')


def find_mrca(clones):
    if len(clones) == 1:
        return 5
    clone_ids = [i[1] for i in clones]
    for gen in range(len(clone_ids[0])):
        first = clone_ids[0][gen]
        for other_clone in clone_ids[1:]:
            # print(first, other_clone[gen])
            if other_clone[gen] != first:
                return min(gen, 5)
    return min(gen, 5)

def mrca_1_check_all_4(clones_df):
    clone_ids_0 = [i[1][0] for i in clones_df.values]
    clone_ids_1_2 = [i[1][1:3] for i in clones_df.values]
    unique_clone_ids_0 = np.unique(clone_ids_0)
    unique_clone_ids_1_2 = np.unique(clone_ids_1_2)
    if len(unique_clone_ids_0) == 1 and len(unique_clone_ids_1_2) == 4:
        return True
    return False

def division_or_cell_death(cells, selection_coefs, driver_strands, epistatic_interactions, require_two_strands, error_free_rate, error_free_tx, repair_rate, shared_mut_through_repair, repair_process, asymmetric_rate_0, asymmetric_rate, smax, verbose, record_detailed=False):
    updated_cells = []
    for cell in cells:
        if repair_rate != 0:
            cell.repair_lesions(repair_rate, shared_mut_through_repair, repair_process)

        decision = cell.division_decision(selection_coefs, driver_strands, epistatic_interactions, require_two_strands, asymmetric_rate_0, asymmetric_rate, error_free_tx, smax)

        if record_detailed and verbose:
            logging.info(f'\nCell {cell.id}: decision = {decision}')

        if decision != 'death':
            if record_detailed and verbose:
                logging.info('Daughter cell information:')
            cell1, cell2 = cell.divide(error_free_rate, propagate_id=record_detailed)
            if decision == 'birth':
                updated_cells.extend([cell1, cell2])
                if record_detailed and verbose:
                    log_cell_info(cell1, selection_coefs, epistatic_interactions, require_two_strands, error_free_tx, smax)
                    log_cell_info(cell2, selection_coefs, epistatic_interactions, require_two_strands, error_free_tx, smax)
            else:
                updated_cells.append(cell1)
                if record_detailed and verbose:
                    log_cell_info(cell1, selection_coefs, epistatic_interactions, require_two_strands, error_free_tx, smax)

    return updated_cells

def simple_division_or_cell_death(cells_ct_tuples, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, asymmetric_rate_0, asymmetric_rate, smax, verbose, record_detailed=False):
    updated_cells = []
    for cell, ct in cells_ct_tuples:
        new_ct = 0
        for cell_i in range(ct):
            decision = cell.division_decision(selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate, error_free_tx, smax)
            if record_detailed and verbose:
                logging.info(f'\nCell {cell.id} ind {cell_i}: decision = {decision}')
            if decision == 'birth':
                new_ct += 2
            elif decision == 'asymmetric':
                new_ct += 1
        if new_ct > 0:
            updated_cells.append((cell, new_ct))
    return updated_cells


def get_smax(selection_coefs):
    w = 1
    for chrom, pos_dict in selection_coefs.items():
        for pos, s in pos_dict.items():
            w *= (1+s)
    return max(1, w-1)


def cell_division_until_drivers_saturated(chr_sizes, selection_coefs, driver_lesion_dict, driver_strands, epistasis_type, require_two_strands, error_free_rate, error_free_tx, repair_rate, shared_mut_through_repair, repair_process, asymmetric_rate_0, asymmetric_rate, n_track_cells, min_gens, verbose):
    generation = 0
    cell = Cell(chr_sizes, '')
    cell.den_exposure(driver_lesion_dict)

    s_max = get_smax(selection_coefs)

    cell.repair_lesions(repair_rate, shared_mut_through_repair, repair_process)
    drivers, lesions, shared_drivers = cell.get_lesion_counts(selection_coefs)
    logging.info(f'Generation {generation}: lesion_count = {lesions}, driver lesion count = {drivers}')

    drivers_repaired_first_cell = False
    duplex_first_cell = False
    # check all drivers repaired?
    if not cell.has_driver_lesions(selection_coefs):
        drivers_repaired_first_cell = True
    if len(cell.get_driver_mutations(selection_coefs)) > 0:
        duplex_first_cell = True

    generation += 1
    ## NOTE: first cell is repaired already, feeding repair rate as 0 here

    cells = division_or_cell_death(cells=[cell], selection_coefs=selection_coefs, driver_strands=driver_strands, epistatic_interactions=epistasis_type,
                                   require_two_strands=require_two_strands, error_free_rate=error_free_rate, error_free_tx=error_free_tx,
                                   repair_rate=0, shared_mut_through_repair=0, repair_process=repair_process,
                                   asymmetric_rate_0=asymmetric_rate_0, asymmetric_rate=asymmetric_rate, smax=s_max,
                                   verbose=verbose, record_detailed=True)


    if verbose:
        logging.info(f'\nGeneration {generation}:')
        for cell in cells:
            log_cell_info(cell, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, s_max)

    if len(cells) == 0:
        lost_first_cell = True
        return cells, generation, shared_drivers, s_max, drivers_repaired_first_cell, duplex_first_cell, lost_first_cell

    lost_first_cell = False
    # detailed lineage tracking
    while generation < n_track_cells:
        generation += 1

        if verbose:
            logging.info(f'\nGeneration {generation}:')
        cells = division_or_cell_death(cells, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_rate,
                                       error_free_tx, repair_rate, shared_mut_through_repair, repair_process,
                                       asymmetric_rate_0, asymmetric_rate, s_max, verbose, record_detailed=True)

        if len(cells) == 0:
            break

    if verbose:
        logging.info('-------------------Continuing till driver saturation-------------------')
    # group cells with lesions and cells without
    if driver_lesions_remaining(cells, selection_coefs):
        cells_w_lesions = []
        cells_no_lesions = []
        for cell in cells:
            if cell.has_driver_lesions(selection_coefs):
                cells_w_lesions.append(cell)
            else:
                cells_no_lesions.append((cell, 1))
    else:
        cells_w_lesions = []
        cells_no_lesions = [(cell, 1) for cell in cells]

    while driver_lesions_remaining(cells_w_lesions, selection_coefs):
        generation += 1
        # propagate cells with lesions
        cells_w_lesions = division_or_cell_death(cells_w_lesions, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_rate,
                                                error_free_tx, repair_rate, shared_mut_through_repair, repair_process,
                                                asymmetric_rate_0, asymmetric_rate, s_max, verbose, record_detailed=True)
        # propagate_cells_without
        cells_no_lesions = simple_division_or_cell_death(cells_no_lesions, selection_coefs, driver_strands, epistasis_type, require_two_strands,
                                                         error_free_tx, asymmetric_rate_0, asymmetric_rate, s_max, verbose, record_detailed=True)
        if len(cells) == 0:
            break
        # update lesions and no lesions list
        new_cells_w_lesions = []
        for cell in cells_w_lesions:
            if cell.has_driver_lesions(selection_coefs):
                new_cells_w_lesions.append(cell)
            else:
                cells_no_lesions.append((cell, 1))
        cells_w_lesions = new_cells_w_lesions

    if verbose:
        logging.info(f'-------------------Continuing for {min_gens} generations-------------------')

    while generation < min_gens:
        generation += 1
        # use simple update
        cells_no_lesions = simple_division_or_cell_death(cells_no_lesions, selection_coefs, driver_strands, epistasis_type, require_two_strands,
                                                     error_free_tx, asymmetric_rate_0, asymmetric_rate, s_max, verbose)
        if verbose:
            logging.info(f'Generation {generation}, #cells = {len(cells)}, drivers remaning? {driver_lesions_remaining(cells, selection_coefs)}, any lesions? {lesions_remaining(cells)}')
        if len(cells) == 0:
            break
    return cells_no_lesions, generation, shared_drivers, s_max, drivers_repaired_first_cell, duplex_first_cell, lost_first_cell


def clonal_expansion(cells, max_generations, carrying_capacity, selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate, error_free_tx, smax, verbose):
    n_0 = 0
    for cell, ct in cells:
        n_0 += ct

    growth_fn = lambda x, y: n_0*(x**y) # x is the growth rate, y is the time step

    cell_expansion = [[cell, cell.id,
                       cell.get_fitness(selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, smax),
                       cell.get_branching_process_probabilities(selection_coefs, driver_strands, epistasis_type, require_two_strands,
                                                                asymmetric_rate_0, asymmetric_rate, error_free_tx, smax),
                       ct/n_0, ct] for cell, ct in cells]

    ## RECORD CLONES
    logging.info('Surviving Cells:')
    surviving_cells_info = []
    for group in cell_expansion:
        logging.info(f'Cell {group[1]}, fitness = {group[2]}, proportion = {group[4]}')
        surviving_cells_info.append([group[1], group[2], group[3], group[4], group[5]])
    surviving_cells_info_df = pd.DataFrame(surviving_cells_info, columns=['clone id', 'clone fitness', 'b,d,c', 'initial clone proportion', 'initial clone size'])

    for gen in range(1, max_generations+1):
        w_bar = get_average_fitness(cell_expansion)
        N_t = growth_fn(w_bar, gen) # feed in rate
        cell_prob = get_cell_propagation_probability(cell_expansion, w_bar)
        try:
            new_generation = np.random.multinomial(N_t, cell_prob)
        except OverflowError:
            break
        cell_expansion = update_cell_information(cell_expansion, new_generation)
        # remove clones that disappear
        cell_expansion = [clone for clone in cell_expansion if clone[-1] > 0]
        if N_t > carrying_capacity:
            break
    if max_generations == 0:
        gen = 0
    surviving_cells_info_df['final clone proportion'] = 0
    surviving_cells_info_df['final clone size'] = 0

    for group in cell_expansion:
        surviving_cells_info_df.loc[surviving_cells_info_df['clone id'] == group[1], 'final clone proportion'] = group[-2]
        surviving_cells_info_df.loc[surviving_cells_info_df['clone id'] == group[1], 'final clone size'] = group[-1]

    return cell_expansion, gen, surviving_cells_info_df

def full_process(chr_sizes, selection_coefs, driver_lesions, driver_strands, require_two_strands, epistasis_type, error_free_rate, error_free_tx, repair_rate, shared_mut_through_repair, repair_process, asymmetric_rate_0, asymmetric_rate, n_track_cells, min_detailed_gens, max_gens, carrying_capacity, outdir, verbose):
    cells, ngenerations, shared_drivers, smax, drivers_repaired_immediately, duplex_first_cell, lost_first_cell = cell_division_until_drivers_saturated(chr_sizes, selection_coefs, driver_lesions, driver_strands, epistasis_type, require_two_strands, error_free_rate, error_free_tx, repair_rate, shared_mut_through_repair, repair_process, asymmetric_rate_0 ,asymmetric_rate, n_track_cells, min_detailed_gens, verbose)
    # Note cells are now in form tuple (cell, count)
    if len(cells) > 0:
        clones, total_gens, surviving_cells = clonal_expansion(cells, max_gens, carrying_capacity, selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate, error_free_tx, smax, verbose)
        # Remove small clones
        #   First group by drivers
        clones_by_driver = group_clones_by_drivers(clones, selection_coefs)
        clones_to_remove = [j for i in clones_by_driver if i[-1] < CLONE_MIN_PROPORTION for j in i[0]]
        #   Then remove small clones defined by undetectable driver
        clones_filtered = [clone for clone in clones if clone[1] not in clones_to_remove]
        # Also try filtering without grouping
        clones_filtered_no_grouping = [clone for clone in clones if clone[-2] >= CLONE_MIN_PROPORTION]

        logging.info('Done with detailed expansion')
        logging.info(f'Final clone information after {total_gens} generation:')
        logging.info(clones)

        # Get mrca
        mrca_filtered_driver = None if len(clones_filtered) == 0 else find_mrca(clones_filtered)
        mrca_filtered_no_grouping = None if len(clones_filtered_no_grouping) == 0 else find_mrca(clones_filtered_no_grouping)


        clone_info_filtered_driver = []
        driver_counts_list_filtered_driver = []
        for clone in clones_filtered:
            driver_count, mut_count = clone[0].get_mutation_counts(selection_coefs)
            clone_info_filtered_driver.append([clone[1], clone[2], clone[3], clone[4], clone[5], total_gens, (driver_count, mut_count)])
            driver_counts_list_filtered_driver.append(int(driver_count))
        clone_info_filtered_driver_df = pd.DataFrame(clone_info_filtered_driver, columns=['clone id', 'clone fitness', 'birth, death, const rates',
                                                                                          'clone proportion', 'clone size', 'n_generations',
                                                                                          '(driver , all mutation) count'])

        clone_info_filtered_no_grouping = []
        driver_counts_list_filtered_no_grouping = []
        for clone in clones_filtered_no_grouping:
            driver_count, mut_count = clone[0].get_mutation_counts(selection_coefs)
            clone_info_filtered_no_grouping.append([clone[1], clone[2], clone[3], clone[4], clone[5], total_gens, (driver_count, mut_count)])
            driver_counts_list_filtered_no_grouping.append(int(driver_count))
        clone_info_filtered_np_grouping_df = pd.DataFrame(clone_info_filtered_no_grouping, columns=['clone id', 'clone fitness', 'birth, death, const rates',
                                                                                          'clone proportion', 'clone size', 'n_generations',
                                                                                          '(driver , all mutation) count'])
        clone_info_df_pre_filter = pd.DataFrame(np.array(clones, dtype=object)[:, 1:], columns=['clone id', 'clone fitness', 'birth, death rates', 'clone proportion', 'clone size'])

        return True, clone_info_filtered_driver_df, clone_info_filtered_np_grouping_df, clone_info_df_pre_filter, surviving_cells, cells, \
               driver_counts_list_filtered_driver, driver_counts_list_filtered_no_grouping, mrca_filtered_driver, mrca_filtered_no_grouping, \
               drivers_repaired_immediately, duplex_first_cell, lost_first_cell
    else:
        logging.info('No cells remaining')
        return False, None, None, None, None, None, None, None, None, None, drivers_repaired_immediately, duplex_first_cell, lost_first_cell

def plot_histogram(data, axis_labels, plotname, **kwargs):
    fig, ax = plt.subplots()
    ax = sns.histplot(data, **kwargs)
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(plotname, bbox_inches='tight')
    plt.close()

def plot_histogram_mrca(data, axis_labels, plotname, **kwargs):
    ax = sns.histplot(data, binrange=(0, 5), **kwargs)
    perc_to_annotate = [0 for i in range(6)]
    mrca_val, counts = np.unique(data, return_counts=True)

    for m, c in zip(mrca_val, counts):
        perc_to_annotate[m] = round(c/sum(counts), 3)
    for ind in range(6):
        patch = ax.patches[ind]
        ax.annotate(perc_to_annotate[ind],
                    (patch.get_x() + patch.get_width()/3, patch.get_height()+0.02))

    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.ylim(0, 1)
    plt.savefig(plotname, bbox_inches='tight')
    plt.close()

def plot_histogram_clone_size(data, axis_labels, plotname):
    fig, ax = plt.subplots()
    ax = sns.histplot(data=data.reset_index(), x='final clone size', log_scale=True) #, binrange=(0, max(data['final clone size'])))
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(plotname, bbox_inches='tight')
    plt.close()


def plot_scatter_clone_prop(data, axis_labels, plotname):
    fig, ax = plt.subplots()
    ax = sns.scatterplot(data=data, x='initial clone proportion', y='final clone proportion', hue='clone fitness')
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(plotname, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir', '-wd', required=True, help='path to base directory')
    parser.add_argument('--cell_dir', '-c', required=True, help='path to file with cell driver info to use from wd')
    parser.add_argument('--lesions_dir', '-ld', required=True, help='path to drawn lesions to use from wd')
    parser.add_argument('--epistasis_type', '-e', required=False, default='add', choices=['add', 'syn', 'ant'], help='epistasis type to use for multiple drivers')
    parser.add_argument('--driver_strand', '-ds', required=False, default='tx', choices=['tx', 'nontx', 'both'], help='whether to treat driver lesions as on transcribed, non transcribed strand, or a combination'
                                                                                                                      'if both specified 1 driver or 3+, each driver will have 0.5 chance of being on each strand, '
                                                                                                                      'for two drivers, each will be assigned to a different strand')
    parser.add_argument('--two_strands_required', '-t', required=False, default='True', help='whether to require two strands to have a mutation to induce a fitness boost')
    parser.add_argument('--error_free_replication', '-er', required=False, default=0.25, help='rate of error free replication')
    parser.add_argument('--error_free_transcription', '-et', required=False, default=1, choices=['0', '1'], help='rate of error free transcription over lesions (erroneous transcription over a driver lesion will lead to fitness change')
    parser.add_argument('--repair_rate', '-rr', required=False, default=0.5, help='rate of lesion repair')
    parser.add_argument('--shared_mutations_through_repair', '-sr', required=False, default=0.05, help='rate of shared mutations to induce from lesions through repair')
    parser.add_argument('--repair_process', '-rp', required=False, default='lesion_mut', choices=['lesion_mut', 'mut_mut'], help='process for getting a mutation through repair')
    parser.add_argument('--asymmetric_division_prob_start', '-ad0', required=False, default=0, help='probability of asymmetric division default')
    parser.add_argument('--asymmetric_division_prob', '-ad', required=False, default=0, help='probability of asymmetric division with driver')
    parser.add_argument('--n_track_segregations', '-n', required=False, default=5, help='number of generations to track segregation for each cell')
    parser.add_argument('--min_detailed_generations', '-min', required=False, default=10, help='min generations to run detailed process for')
    parser.add_argument('--max_generations', '-max', required=False, default=300, help='max generations to run clonal expansion for')
    parser.add_argument('--carrying_capacity', '-cap', required=False, default=10**8, help='tumor carrying capacity')
    parser.add_argument('--min_tumor_size', '-mt', required=False, default=10**6, help='min number of cells to consider a tumor have formed')
    parser.add_argument('--outprefix', '-o', required=False, default='', help='additional prefix to add to out directory')
    parser.add_argument('--seed_min', '-smin', required=False, default=None, help='min random seed to use in range. if min and max not set, seed will be randomized')
    parser.add_argument('--seed_max', '-smax', required=False, default=None, help='max random seed to use in range')
    parser.add_argument('--sim_count', '-sc', required=False, default=1, help='number of simulations to run')
    parser.add_argument('--verbose', '-v', required=False, choices=['True', 'False'], default=False, help='whether to write detailed files for each simulation')
    args = parser.parse_args()

    # Process arguments
    epistasis_dict = {'add': 'NO EPISTASIS', 'syn':'SYNERGESTIC', 'ant':'ANTAGONISTIC'}
    epistasis_type = epistasis_dict[args.epistasis_type]
    require_two_strands = True if args.two_strands_required == 'True' else False
    error_free_replication_rate = float(args.error_free_replication)
    error_free_transcription_rate = float(args.error_free_transcription)
    repair_rate = float(args.repair_rate)
    shared_mutations_through_repair = float(args.shared_mutations_through_repair)
    asymmetric_rate_0 = float(args.asymmetric_division_prob_start)
    asymmetric_rate = float(args.asymmetric_division_prob)
    n_track_segregation = float(args.n_track_segregations)
    min_detailed_gens = int(args.min_detailed_generations)
    max_gens = int(args.max_generations)
    carrying_capacity = int(args.carrying_capacity)
    min_tumor = int(args.min_tumor_size)
    min_seed = None if (args.seed_min is None or args.seed_min == 'None') else int(args.seed_min)
    max_seed = None if (args.seed_max is None or args.seed_max == 'None') else int(args.seed_max)
    seed_range = None if min_seed is None else range(min_seed, max_seed)
    num_simulations = int(args.sim_count)
    verbose = True if args.verbose == 'True' else False

    if repair_rate < shared_mutations_through_repair:
        raise ValueError("repair rate cannot be less than rate of shared mutations through repair")

    # create directory and set up out files
    sim_dir = f'{args.working_dir}/{args.lesions_dir}/simulation_results/'
    if not os.path.isdir(sim_dir):
        os.mkdir(sim_dir)
    sim_dir = f'{sim_dir}/epistasis_{args.epistasis_type}_driver_on_{args.driver_strand}_strands_' \
              f'efr_{error_free_replication_rate}_eft_{error_free_transcription_rate}_2strands_{require_two_strands}_' \
              f'repair_{repair_rate}_{shared_mutations_through_repair}_{args.repair_process}'
    if not os.path.isdir(sim_dir):
        os.mkdir(sim_dir)
    day = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    outprefix = f'_{args.outprefix}' if args.outprefix != '' else ''
    sim_dir = f'{sim_dir}/{day}{outprefix}/'
    if not os.path.isdir(sim_dir):
        os.mkdir(sim_dir)

    chr_sizes = np.loadtxt(f'{args.working_dir}/{args.cell_dir}/chr_sizes.txt', dtype=int)
    selection_coefs = pickle.load(open(f'{args.working_dir}/{args.lesions_dir}/driver_lesion_selection_coefficients.pickle', 'rb'))
    driver_lesions = pickle.load(open(f'{args.working_dir}/{args.lesions_dir}/driver_lesions.pickle', 'rb'))
    driver_strands = {}
    prev = 'tx'
    for chr, mdict in selection_coefs.items():
        if mdict:
            for pos, s in mdict.items():
                if chr not in driver_strands:
                    driver_strands[chr] = {}
                if args.driver_strand == 'both':
                    if prev == 'tx':
                        driver_strands[chr][pos] = 'nontx'
                        prev = 'nontx'
                    else:
                        driver_strands[chr][pos] = 'tx'
                        prev = 'tx'
                else:
                    driver_strands[chr][pos] = args.driver_strand

    # write general information
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_handler = logging.FileHandler(f'{sim_dir}/output.log', 'w')
    logger.addHandler(logger_handler)
    logger.info(f'Simulation with cell from dir {args.cell_dir}')
    logger.info(f'Epistasis type: {epistasis_type}')
    logger.info(f'Den lesion introduction using the lesion introduced in {args.lesions_dir}')
    logger.info(f'Lesion considered to be on {args.driver_strand} strand')
    logger.info(f'Error free replication rate {error_free_replication_rate}')
    logger.info(f'Error free transcription rate {error_free_transcription_rate}')
    logger.info(f'Rate of repair {repair_rate}')
    logger.info(f'Rate of induced shared mutations from lesions through repair {shared_mutations_through_repair}')
    logger.info(f'Process of producing a mutation through repair {args.repair_process}')
    logger.info(f'Asymmetric division rate set to {asymmetric_rate_0} in base cell and {asymmetric_rate} when there is a driver')
    logger.info(f'Requiring both strands to carry a driver mutation to induce a fitness boost? {require_two_strands}')
    logger.info(f'Tracking cell division in detail until driver lesions are all removed, at least for first {n_track_segregation} generations')
    logger.info(f'Running detailed generations until at least {min_detailed_gens} generations')
    logger.info(f'Running clonal expansion for {max_gens} generations with a carrying capacity of {carrying_capacity}')
    logger.info(f'Random seed range {min_seed} - {max_seed}')
    if verbose:
        logger_handler.close()
        logger.handlers.pop()


    sims_compeleted = []
    sims_min = []

    count_drivers_rep_immediately = 0
    count_duplex_first_cell = 0
    count_lost_first_cell = 0

    count_detectable_min = 0

    clone_info_filtered_driver_all = []
    clone_info_filtered_no_grouping_all = []
    clone_info_pre_filter_all = []
    surviving_cells_all = []

    driver_counts_filtered_driver_all = []
    driver_counts_filtered_no_grouping_all = []

    count_all_1_driver = 0
    count_some_1_driver = 0
    count_no_1_driver = 0

    mrca_list_filtered_driver = []
    mrca_list_filtered_no_grouping = []
    mrca_1 = [] # all of these will be based on filtering without grouping
    mrca_1_all_4 = []
    mrca_1_driver = []
    mrca_2_drivers = []
    mrca_1_2_drivers = []

    fitness_list = []
    tumor_sizes_all_clones = []
    tumor_sizes = []

    random_seeds_used = set()
    for i in range(num_simulations):
        if seed_range is None:
            seed = np.random.randint(0, 10*num_simulations)
            while seed in random_seeds_used:
                seed = np.random.randint(0, 10*num_simulations)
        else:
            seed = seed_range[i]

        random_seeds_used.add(seed)
        np.random.seed(seed)

        if verbose: # set seperate output directory for simulation
            sim_dir_i = f'{sim_dir}/{i}/'
            if not os.path.isdir(sim_dir_i):
                os.mkdir(sim_dir_i)
            logger = logging.getLogger()
            logger.setLevel(logging.INFO)
            logger_handler = logging.FileHandler(f'{sim_dir_i}/output.log', 'w')
            logger.addHandler(logger_handler)
            logger.info(f'Simulation with seed {seed}')
        else:
            if seed%1000 == 0:
                logger.info(f'Simulation with random seed {seed}')

        sim_outprefix = sim_dir_i if verbose else sim_dir
        completed, clone_info_filtered_driver, clone_info_filtered_no_grouping, clone_info_pre_filter, surviving_cell_info, starting_cell_groups, \
        driver_counts_filtered_driver, driver_counts_filtered_no_grouping, mrca_filtered_driver, mrca_filtered_no_grouping, \
        drivers_rep_imm, duplex_first_cell, lost_first_cell = full_process(chr_sizes, selection_coefs, driver_lesions, driver_strands, require_two_strands, epistasis_type,
                                                                                         error_free_replication_rate, error_free_transcription_rate,
                                                                                         repair_rate, shared_mutations_through_repair, args.repair_process,
                                                                                         asymmetric_rate_0, asymmetric_rate,
                                                                                         n_track_segregation, min_detailed_gens, max_gens, carrying_capacity,
                                                                                         sim_outprefix, verbose)

        if drivers_rep_imm: count_drivers_rep_immediately += 1
        if duplex_first_cell: count_duplex_first_cell += 1
        if lost_first_cell: count_lost_first_cell += 1
        if completed:
            clone_info_filtered_driver.insert(0, 'simulation', i)
            clone_info_filtered_driver_all.append(clone_info_filtered_driver)

            clone_info_filtered_no_grouping.insert(0, 'simulation', i)
            clone_info_filtered_no_grouping_all.append(clone_info_filtered_no_grouping)

            clone_info_pre_filter.insert(0, 'simulation', i)
            clone_info_pre_filter_all.append(clone_info_pre_filter)

            surviving_cell_info.insert(0, 'simulation', i)
            surviving_cells_all.append(surviving_cell_info)

            total_clone_size_all = sum(surviving_cell_info['final clone size'].values)
            tumor_sizes_all_clones.append(total_clone_size_all)
            total_clone_size = sum(clone_info_filtered_no_grouping['clone size'].values)

            if total_clone_size > 0:
                tumor_sizes.append(total_clone_size)

            if total_clone_size >= min_tumor:
                count_detectable_min += 1
                sims_min.append(i)
                if 0 in driver_counts_filtered_driver or 0 in driver_counts_filtered_no_grouping:
                    print(i, driver_counts_filtered_driver, driver_counts_filtered_no_grouping)
                driver_counts_filtered_driver_all.extend(driver_counts_filtered_driver)
                driver_counts_filtered_no_grouping_all.extend(driver_counts_filtered_no_grouping)

                mrca_list_filtered_driver.append(mrca_filtered_driver)
                mrca_list_filtered_no_grouping.append(mrca_filtered_no_grouping)

                if mrca_filtered_no_grouping == 1:
                    mrca_1.append(i)
                    if mrca_1_check_all_4(clone_info_filtered_no_grouping):
                        mrca_1_all_4.append(i)
                if sum(driver_counts_filtered_no_grouping) == len(driver_counts_filtered_no_grouping):
                    count_all_1_driver += 1
                    if mrca_filtered_no_grouping is not None:
                        mrca_1_driver.append(mrca_filtered_no_grouping)
                elif 1 not in driver_counts_filtered_no_grouping:
                    count_no_1_driver += 1
                    if mrca_filtered_no_grouping is not None:
                        mrca_2_drivers.append(mrca_filtered_no_grouping)
                else:
                    count_some_1_driver += 1
                    if mrca_filtered_no_grouping is not None:
                        mrca_1_2_drivers.append(mrca_filtered_no_grouping)

                fitness = (clone_info_filtered_no_grouping['clone fitness']-1).values
                fitness_list.extend(fitness)

        if verbose:
            logger_handler.close()
            logger.handlers.pop()


    num_completed = len(sims_compeleted)
    sims_compeleted.append(f'total simulations run = {num_simulations}')
    sims_compeleted.append(f'count that reached clonal expansion = {num_completed}')
    sims_compeleted.append(f'count that reached at least {min_tumor} cells = {count_detectable_min}')
    sims_compeleted.append(f'count MRCA 1 = {len(mrca_1)}')
    sims_compeleted.append(f'count MRCA 1 with all 4 cells contributing to tumor = {len(mrca_1_all_4)}')
    sims_compeleted.append(f'count all clones have 1 driver = {count_all_1_driver}')
    sims_compeleted.append(f'count all clones have more than 1 driver = {count_no_1_driver}')
    sims_compeleted.append(f'count some clones have 1 driver and some have more than 1 driver = {count_some_1_driver}')
    sims_compeleted.append(f'count with driver lesions repaired immediately = {count_drivers_rep_immediately}')
    sims_compeleted.append(f'count with duplex created in first cell = {count_duplex_first_cell}')
    sims_compeleted.append(f'count with first cell undergoing d = {count_lost_first_cell}')
    np.savetxt(f'{sim_dir}/completed_simulations.txt', sims_compeleted, fmt='%s')
    np.savetxt(f'{sim_dir}/simulations_mrca_1.txt', mrca_1, fmt='%s')
    np.savetxt(f'{sim_dir}/simulations_mrca_1_all_4_survived.txt', mrca_1_all_4, fmt='%s')
    clone_info_all_filtered_driver_df = pd.concat(clone_info_filtered_driver_all)
    clone_info_all_filtered_driver_df.to_csv(f'{sim_dir}/combined_clone_info_filtered_driver.csv', index=False)

    clone_info_all_filtered_no_grouping_df = pd.concat(clone_info_filtered_no_grouping_all)
    clone_info_all_filtered_no_grouping_df.to_csv(f'{sim_dir}/combined_clone_info_filtered_no_grouping.csv', index=False)

    clone_info_pre_filter_all_df = pd.concat(clone_info_pre_filter_all)
    clone_info_pre_filter_all_df.to_csv(f'{sim_dir}/combined_clone_info_no_filter.csv', index=False)

    surviving_cells_all_df = pd.concat(surviving_cells_all)
    surviving_cells_all_df.to_csv(f'{sim_dir}/surviving_cells_info.csv', index=False)

    plot_histogram(driver_counts_filtered_driver_all, ['num drivers', 'clone proportion'], f'{sim_dir}/driver_counts_per_clone_filtered_driver.pdf',
                   binrange=(0, max(driver_counts_filtered_driver_all)+1), discrete=True, stat='probability')
    plot_histogram(driver_counts_filtered_no_grouping_all, ['num drivers', 'clone proportion'], f'{sim_dir}/driver_counts_per_clone_filtered_no_grouping.pdf',
                   binrange=(0, max(driver_counts_filtered_no_grouping_all)+1), discrete=True, stat='probability')

    plot_histogram_mrca(mrca_list_filtered_driver, ['MRCA generation', 'tumor proportion'],
                        f'{sim_dir}/mrca_generation_counts_{min_tumor}_cells_filtered_driver.pdf',
                        discrete=True, stat='probability')

    if mrca_list_filtered_no_grouping:
        plot_histogram_mrca(mrca_list_filtered_no_grouping, ['MRCA generation', 'tumor proportion'],
                            f'{sim_dir}/mrca_generation_counts_{min_tumor}_cells_filtered_no_grouping.pdf',
                            discrete=True, stat='probability')

    if mrca_2_drivers or mrca_1_2_drivers:
        if mrca_1_driver:
            plot_histogram_mrca(mrca_1_driver, ['MRCA generation', 'tumor proportion'],
                            f'{sim_dir}/mrca_generation_counts_{min_tumor}_cells_1_driver_tumors_filtered_no_grouping.pdf',
                            discrete=True, stat='probability')
        if mrca_2_drivers:
            plot_histogram_mrca(mrca_2_drivers, ['MRCA generation', 'tumor proportion'],
                            f'{sim_dir}/mrca_generation_counts_{min_tumor}_cells_2plus_driver_tumors_filtered_no_grouping.pdf',
                            discrete=True, stat='probability')
        if mrca_1_2_drivers:
            plot_histogram_mrca(mrca_1_2_drivers, ['MRCA generation', 'tumor proportion'],
                            f'{sim_dir}/mrca_generation_counts_{min_tumor}_cells_1_and_more_driver_clones_tumors_filtered_no_grouping.pdf',
                            discrete=True, stat='probability')

    plot_histogram(fitness_list, ['Fitness advantage', 'clone proportion'], f'{sim_dir}/fitness_advantage_per_clone.pdf', bins=100,
                   binrange=(0, 1), stat='probability')

    # write files for plots
    yaml.dump(driver_counts_filtered_driver_all, open(f'{sim_dir}driver_counts_filtered_driver.yml', 'w'))
    yaml.dump(driver_counts_filtered_no_grouping_all, open(f'{sim_dir}driver_counts_filtered_no_grouping.yml', 'w'))

    yaml.dump(mrca_list_filtered_driver, open(f'{sim_dir}mrca_list_{min_tumor}_filtered_driver.yml', 'w'))
    yaml.dump(mrca_list_filtered_no_grouping, open(f'{sim_dir}mrca_list_{min_tumor}_filtered_no_grouping.yml', 'w'))

    plot_histogram(tumor_sizes, [f'tumor size (clones >= {CLONE_MIN_PROPORTION} of cells)', 'proportion'], f'{sim_dir}/tumor_size_sizable_clones.pdf',
                   log_scale=True, stat='probability')

    plot_scatter_clone_prop(surviving_cells_all_df, ['clone proportion before clonal expansion',
                                          'clone proportion after clone expansion'], f'{sim_dir}/clone_proportion_change.pdf')
    surviving_cells_all_df = surviving_cells_all_df.loc[surviving_cells_all_df['final clone size'] > 0]
    plot_histogram_clone_size(surviving_cells_all_df,
                              ['clone size', 'clone count'], f'{sim_dir}/clone_sizes_all.pdf')


if __name__ == '__main__':
    main()
