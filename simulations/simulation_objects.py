import copy

import numpy as np


from utils import  binomial_rv


NEOPLASM_THRESHOLD = 1.02
TRANSFORMATION_THRESHOLD = 1.05


class Strand():
    def __init__(self, chr_size, lesions = None, mutations = None):
        self.chr_size = chr_size
        if mutations is None:
            self.lesion_locations = set()
            self.mutation_locations = set()
            self.lesions = self.has_lesions()
        else:
            self.lesion_locations = lesions if lesions is not None else set()
            self.mutation_locations = mutations #if mutations is not None else set()
            self.lesions = self.has_lesions()

    def __str__(self):
        return f'lesion locations: {self.lesion_locations}\nmutation locations: {self.mutation_locations}\nlesion status: {self.lesions}'

    def __eq__(self, other):
        '''
        Considers strands with same length and same lesion and mutation locations to be equal
        :param other:
        :return:
        '''
        if isinstance(other, Strand):
            if self.lesion_locations == other.lesion_locations and self.mutation_locations == other.mutation_locations and self.chr_size == other.chr_size:
                return True

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def update_lesions(self, lesions):
        self.lesion_locations = lesions

    def get_lesions(self):
        return copy.deepcopy(self.lesion_locations)

    def get_mutations(self):
        return copy.deepcopy(self.mutation_locations)

    def get_lesions(self):
        return copy.deepcopy(self.lesion_locations)

    def has_lesions(self):
        if self.lesion_locations == None:
            return False
        return len(self.lesion_locations) > 0

    def has_driver_lesions(self, selection_coefs):
        for pos in self.lesion_locations:
            if pos in selection_coefs:
                return True
        return False

    def get_driver_mutations(self, selection_coefs):
        d = []
        for pos in self.mutation_locations:
            if pos in selection_coefs:
                d.append(pos)
        return sorted(d)

    def get_driver_lesions_or_mutations(self, selection_coefs):
        d = set()
        for pos in self.lesion_locations:
            if pos in selection_coefs:
                d.add(pos)
        for pos in self.mutation_locations:
            if pos in selection_coefs:
                d.add(pos)
        return sorted(list(d))

    def get_driver_lesions(self, selection_coefs):
        d = []
        for pos in self.lesion_locations:
            if pos in selection_coefs:
                d.append(pos)
        return sorted(d)

    def replicate(self, error_free_rate):
        '''
        Create a new strand object based on replication of self
        Lesions will lead to mutations 75% of the time
        Mutatoins will always lead to mutations
        No lesions on new strand
        :return:
        '''
        new_mutations = set()
        for lesion in self.lesion_locations:
            wt_chance = binomial_rv(1, error_free_rate)
            # print('replicating wt chance = ', wt_chance)
            if wt_chance == 0:
                new_mutations.add(lesion)
        new_mutations.update(self.mutation_locations)
        return Strand(self.chr_size, lesions=set(), mutations=new_mutations)

    def repair_lesions(self, repair_rate, opp_strand, desired_shared_mutations, repair_process):
        new_mutations = set()
        lesions_to_remove = set()
        opp_strand_mutations = set()
        repair_to_shared_mut_rate = desired_shared_mutations/repair_rate

        for lesion in self.lesion_locations:
            repair = binomial_rv(1, repair_rate)
            if repair == 1:
                repair_to_shared_mut = binomial_rv(1, repair_to_shared_mut_rate)


                if lesion in opp_strand.mutation_locations or lesion in opp_strand.lesion_locations:
                    new_mutations.add(lesion)
                    lesions_to_remove.add(lesion)
                elif repair_to_shared_mut == 1:
                        opp_strand_mutations.add(lesion)
                        if repair_process == 'mut_mut':
                            lesions_to_remove.add(lesion)
                            new_mutations.add(lesion)
                else:
                    lesions_to_remove.add(lesion)

        new_mutations.update(self.mutation_locations)
        new_lesions = self.lesion_locations.difference(lesions_to_remove)

        return new_lesions, new_mutations, opp_strand_mutations

    def get_lesion_counts(self, selection_coefs):
        d = 0
        for pos in self.lesion_locations:
            if pos in selection_coefs:
                d += 1
        return d, len(self.lesion_locations)




class Chromosome():
    def __init__(self, chr_size, strand1=None, strand2=None):
        self.chr_size = chr_size
        if strand1 is None:
            self.strand1 = Strand(self.chr_size)
            self.strand2 = Strand(self.chr_size)
        else:
            self.strand1 = strand1
            self.strand2 = strand2

    def __str__(self):
        return f'Strand 1\n{self.strand1.__str__()}\nStrand 2\n{self.strand2.__str__()}'

    def __eq__(self, other):
        '''
                ## Strand order matters!!!
        :param other:
        :return:
        '''
        if isinstance(other, Chromosome):
            if self.strand1 == other.strand1 and self.strand2 == other.strand2:
                return True

        return False

    def __ne__(self, other):
        return not self.__eq__(other)


    def den_exposure(self, lesions):
        '''

        :param lesions:
        :return:
        '''
        strand1_lesions = lesions[0]
        strand2_lesions = lesions[1]
        self.strand1.update_lesions(strand1_lesions)
        self.strand2.update_lesions(strand2_lesions)
        return strand1_lesions, strand2_lesions

    def repair_lesions(self, repair_rate, desired_shared_mutations, repair_process):
        new_strand1_lesions, new_strand1_mutations, add_strand2_mutations = self.strand1.repair_lesions(repair_rate, self.strand2, desired_shared_mutations, repair_process)
        new_strand2_lesions, new_strand2_mutations, add_strand1_mutations = self.strand2.repair_lesions(repair_rate, self.strand1, desired_shared_mutations, repair_process)
        new_strand1_mutations.update(add_strand1_mutations)
        new_strand2_mutations.update(add_strand2_mutations)
        self.strand1 = Strand(self.chr_size, lesions=new_strand1_lesions, mutations=new_strand1_mutations)
        self.strand2 = Strand(self.chr_size, lesions=new_strand2_lesions, mutations=new_strand2_mutations)






    def has_lesions(self):
        return self.strand1.has_lesions() or self.strand2.has_lesions()

    def has_driver_lesions(self, selection_coefs):
        return self.strand1.has_driver_lesions(selection_coefs) or self.strand2.has_driver_lesions(selection_coefs)

    def get_driver_mutations(self, selection_coefs):
        strand1_drivers = self.strand1.get_driver_mutations(selection_coefs)
        strand2_drivers = self.strand2.get_driver_mutations(selection_coefs)
        drivers = set(strand1_drivers)
        drivers.update(set(strand2_drivers))
        return drivers

    def get_driver_lesions_or_mutations(self, selection_coefs):
        strand1_drivers = self.strand1.get_driver_lesions_or_mutations(selection_coefs)
        strand2_drivers = self.strand2.get_driver_lesions_or_mutations(selection_coefs)
        drivers = set(strand1_drivers)
        drivers.update(set(strand2_drivers))
        return drivers

    def get_two_strand_driver_mutations(self, selection_coefs):
        strand1_drivers = self.strand1.get_driver_mutations(selection_coefs)
        strand2_drivers = self.strand2.get_driver_mutations(selection_coefs)
        drivers = set(strand1_drivers).intersection(set(strand2_drivers))
        return drivers

    def get_two_strand_driver_lesions_or_mutations(self, selection_coefs):
        strand1_drivers = self.strand1.get_driver_lesions_or_mutations(selection_coefs)
        strand2_drivers = self.strand2.get_driver_lesions_or_mutations(selection_coefs)
        drivers = set(strand1_drivers).intersection(set(strand2_drivers))
        return drivers

    def get_driver_mutations_mark_strand(self, selection_coefs):
        strand1_drivers = set(self.strand1.get_driver_mutations(selection_coefs))
        strand2_drivers = set(self.strand2.get_driver_mutations(selection_coefs))
        both = strand1_drivers.intersection(strand2_drivers)
        strand1_drivers.difference_update(both)
        strand2_drivers.difference_update(both)
        return strand1_drivers, strand2_drivers, both

    def get_driver_lesions_or_mutations_mark_strand(self, selection_coefs):
        strand1_drivers = set(self.strand1.get_driver_lesions_or_mutations(selection_coefs))
        strand2_drivers = set(self.strand2.get_driver_lesions_or_mutations(selection_coefs))
        both = strand1_drivers.intersection(strand2_drivers)
        strand1_drivers.difference_update(both)
        strand2_drivers.difference_update(both)
        return strand1_drivers, strand2_drivers, both

    def get_mutations(self):
        strand1_mut = self.strand1.get_mutations()
        strand2_mut = self.strand2.get_mutations()
        muts = set(strand1_mut)
        muts.update(strand2_mut)
        return muts

    def get_driver_lesions(self, selection_coefs):
        strand1_drivers = self.strand1.get_driver_lesions(selection_coefs)
        strand2_drivers = self.strand2.get_driver_lesions(selection_coefs)
        shared = np.intersect1d(strand1_drivers, strand2_drivers)
        drivers = set(strand1_drivers)
        drivers.update(set(strand2_drivers))
        return drivers, shared

    def get_lesions(self):
        strand1_mut = self.strand1.get_lesions()
        strand2_mut = self.strand2.get_lesions()
        muts = set(strand1_mut)
        muts.update(strand2_mut)
        return muts

    def get_lesion_counts(self, selection_coefs):
        '''

        :return: drivers, all
        '''
        d1, a1 = self.strand1.get_lesion_counts(selection_coefs)
        d2, a2 = self.strand2.get_lesion_counts(selection_coefs)
        return d1+d2, a1+a2


    def divide(self, error_free_rate):
        strand1_rep = self.strand1.replicate(error_free_rate)
        strand2_rep = self.strand2.replicate(error_free_rate)
        chrom1 = Chromosome(self.chr_size, self.strand1, strand1_rep)
        chrom2 = Chromosome(self.chr_size, self.strand2, strand2_rep)
        return chrom1, chrom2


class ChromosomePair():
    def __init__(self, chr_size, chrom=None, homolog=None):
        # print('making pair ', chr_size)
        self.chr_size = chr_size
        if chrom is None:
            self.chrom = Chromosome(self.chr_size)
            self.homolog = Chromosome(self.chr_size)

        else:
            self.chrom = chrom
            self.homolog = homolog


    def __str__(self):
        return f'Chromosome\n{self.chrom.__str__()}\nHomolog\n{self.homolog.__str__()}'

    def __eq__(self, other):
        '''
            Pair order , and chrom-homolog order matter
        :param other:
        :return:
        '''
        if isinstance(other, ChromosomePair):
            if self.chrom == other.chrom and self.homolog == other.homolog:
                return True
            elif self.chrom == other.homolog and self.homolog == other.chrom:
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)


    def den_exposure(self, lesions):
        chrom_lesions = lesions[0]
        homolog_lesions = lesions[1]
        chrom_lesions = self.chrom.den_exposure(chrom_lesions)
        homolog_lesions = self.homolog.den_exposure(homolog_lesions)
        return chrom_lesions, homolog_lesions

    def repair_lesions(self, repair_rate, desired_shared_mutation_rate, repair_process):
        self.chrom.repair_lesions(repair_rate, desired_shared_mutation_rate, repair_process)
        self.homolog.repair_lesions(repair_rate, desired_shared_mutation_rate, repair_process)


    def has_lesions(self):
        return self.chrom.has_lesions() or self.homolog.has_lesions()

    def has_driver_lesions(self, selection_coefs):
        return self.chrom.has_driver_lesions(selection_coefs) or self.homolog.has_driver_lesions(selection_coefs)


    def get_driver_mutations(self, selection_coefs):
        chrom_drivers = self.chrom.get_driver_mutations(selection_coefs)
        homolog_drivers = self.homolog.get_driver_mutations(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        return sorted(list(drivers))

    def get_two_strand_driver_mutations(self, selection_coefs):
        chrom_drivers = self.chrom.get_two_strand_driver_mutations(selection_coefs)
        homolog_drivers = self.homolog.get_two_strand_driver_mutations(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        return sorted(list(drivers))

    def get_driver_mutations_mark_strand(self, selection_coefs):
        chrom_drivers1, chrom_drivers2, chrom_both = self.chrom.get_driver_mutations_mark_strand(selection_coefs)
        homolog_drivers1, homolog_drivers2, homolog_both = self.homolog.get_driver_mutations_mark_strand(selection_coefs)
        drivers1 = set(chrom_drivers1)
        drivers1.update(set(homolog_drivers1))
        drivers2 = set(chrom_drivers2)
        drivers2.update(set(homolog_drivers2))
        both = set(chrom_both)
        both.update(set(homolog_both))
        return drivers1, drivers2, both


    def get_driver_lesions_or_mutations(self, selection_coefs):
        chrom_drivers = self.chrom.get_driver_lesions_or_mutations(selection_coefs)
        homolog_drivers = self.homolog.get_driver_lesions_or_mutations(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        return sorted(list(drivers))

    def get_two_strand_driver_lesions_or_mutations(self, selection_coefs):
        chrom_drivers = self.chrom.get_two_strand_driver_lesions_or_mutations(selection_coefs)
        homolog_drivers = self.homolog.get_two_strand_driver_lesions_or_mutations(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        return sorted(list(drivers))

    def get_driver_lesions_or_mutations_mark_strand(self, selection_coefs):
        chrom_drivers1, chrom_drivers2, chrom_both = self.chrom.get_driver_lesions_or_mutations_mark_strand(selection_coefs)
        homolog_drivers1, homolog_drivers2, homolog_both = self.homolog.get_driver_lesions_or_mutations_mark_strand(selection_coefs)
        drivers1 = set(chrom_drivers1)
        drivers1.update(set(homolog_drivers1))
        drivers2 = set(chrom_drivers2)
        drivers2.update(set(homolog_drivers2))
        both = set(chrom_both)
        both.update(set(homolog_both))
        return drivers1, drivers2, both

    def get_driver_lesions_fitness(self, selection_coefs):
        chrom_drivers, _ = self.chrom.get_driver_lesions(selection_coefs)
        homolog_drivers, _ = self.homolog.get_driver_lesions(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        return sorted(list(drivers))

    def get_driver_info(self, selection_coefs):
        drivers = self.get_driver_mutations(selection_coefs)
        out = []
        for driver in drivers:
            s = selection_coefs[driver]
            out.append([s])
        return out

    def get_mutations(self):
        chrom_mut = self.chrom.get_mutations()
        homolog_mut = self.homolog.get_mutations()
        muts = set(chrom_mut)
        muts.update(homolog_mut)
        return sorted(list(muts))

    def get_driver_lesions(self, selection_coefs):
        chrom_drivers, chrom_shared = self.chrom.get_driver_lesions(selection_coefs)
        homolog_drivers, homolog_shared = self.homolog.get_driver_lesions(selection_coefs)
        drivers = set(chrom_drivers)
        drivers.update(set(homolog_drivers))
        shared = set(chrom_shared)
        shared.update(set(homolog_shared))
        return sorted(list(drivers)), sorted(list(shared))

    def get_lesions(self):
        chrom_mut = self.chrom.get_lesions()
        homolog_mut = self.homolog.get_lesions()
        muts = set(chrom_mut)
        muts.update(homolog_mut)
        return sorted(list(muts))

    def get_lesion_counts(self, selection_coefs):
        '''

        :return: drivers, all
        '''
        d1, a1 = self.chrom.get_lesion_counts(selection_coefs)
        d2, a2 = self.homolog.get_lesion_counts(selection_coefs)
        return d1+d2, a1+a2


    def divide(self, error_free_rate):
        chrom1, chrom2 = self.chrom.divide(error_free_rate)
        homolog1, homolog2 = self.homolog.divide(error_free_rate)
        coin_flip = binomial_rv(1, 0.5)
        if coin_flip == 1:

            pair1 = ChromosomePair(self.chr_size, chrom1, homolog1)
            pair2 = ChromosomePair(self.chr_size, chrom2, homolog2)
        else:

            pair1 = ChromosomePair(self.chr_size, chrom1, homolog2)
            pair2 = ChromosomePair(self.chr_size, chrom2, homolog1)
        return pair1, pair2

class Cell():
    def __init__(self, chr_size_list, id, chromosomes = None):
        '''
        Diploid, for each size in chr_size_list, create a pair (identical).
        Assume pairs already given if chromosomes list is passed in
        :param chr_size_list: length should be equal to one set of chromsomes (1 length per homologous pair)
        :param chromosomes: list of two lists: [[chromosomes], [homologs]]. Assumes homologous pairs are in the same position in each list
        '''
        self.neoplasm = None
        self.transformed = None
        self.chr_size_list = chr_size_list
        self.id = id
        if chromosomes is None:
            chromosomes_list = [ChromosomePair(i) for i in chr_size_list]
            self.chromosomes = chromosomes_list
        else:
            self.chromosomes = chromosomes
            # self.update_state()

    def __str__(self):
        cell_info = f'{self.id}\n'
        chrom_info = '\n'.join(['Chromosome'+str(i+1)+':\n'+self.chromosomes[i].__str__() for i in range(len(self.chromosomes))])
        return cell_info+chrom_info

    def __eq__(self, other):
        '''
        CHROMOSOME ORDER MATTERS
        :param other:
        :return:
        '''
        if isinstance(other, Cell):
            if self.chromosomes == other.chromosomes:
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def flatten_chromosome_lists(self):
        return [chrom for sublist in self.chromosomes for chrom in sublist]

    def den_exposure(self, lesions):
        out_lesions = {}
        for chromosome_pair_ind in range(len(self.chromosomes)):
            pair_lesions = lesions[chromosome_pair_ind]
            chromosome_pair = self.chromosomes[chromosome_pair_ind]
            pair_lesions = chromosome_pair.den_exposure(pair_lesions)
            out_lesions[chromosome_pair_ind] = pair_lesions
        return out_lesions

    def get_lesion_counts(self, selection_coefs):
        '''

        :return: driver lesions, all lesions
        '''
        lesions = 0
        drivers = 0
        shared_drivers = 0
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            lesions += len(chrom_pair.get_lesions())
            if chrom_pair_ind in selection_coefs:
                chrom_drivers, chrom_shared_drivers = chrom_pair.get_driver_lesions(selection_coefs[chrom_pair_ind])
            else:
                chrom_drivers, chrom_shared_drivers = chrom_pair.get_driver_lesions({})
            drivers += len(chrom_drivers)
            shared_drivers += len(chrom_shared_drivers)
        return drivers, lesions, shared_drivers

    def repair_lesions(self, repair_rate, desired_shared_mutation_rate, repair_process):
        if repair_rate > 0:
            for chrom_pair in self.chromosomes:
                chrom_pair.repair_lesions(repair_rate, desired_shared_mutation_rate, repair_process)

    def get_mutation_counts(self, selection_coefs):
        '''

        :return: driver lesions, all lesions
        '''
        muts = 0
        drivers = 0
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            muts += len(chrom_pair.get_mutations())
            drivers += len(chrom_pair.get_driver_mutations(selection_coefs[chrom_pair_ind]))
        return drivers, muts

    def get_driver_mutations(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            muts = chrom_pair.get_driver_mutations(selection_coefs[chrom_pair_ind])
            if len(muts) > 0:
                out[chrom_pair_ind] = muts
        return out

    def get_two_strand_driver_mutations(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            muts = chrom_pair.get_two_strand_driver_mutations(selection_coefs[chrom_pair_ind])
            if len(muts) > 0:
                out[chrom_pair_ind] = muts
        return out

    def get_driver_mutations_mark_strand(self, selection_coefs):
        out1 = {}
        out2 = {}
        out_both = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            mut_strand1, mut_strand2, mut_both = chrom_pair.get_driver_mutations_mark_strand(selection_coefs[chrom_pair_ind])
            if len(mut_strand1) > 0:
                out1[chrom_pair_ind] = mut_strand1
            if len(mut_strand2) > 0:
                out2[chrom_pair_ind] = mut_strand2
            if len(mut_both) > 0:
                out_both[chrom_pair_ind] = mut_both
        return out1, out2, out_both


    def get_driver_lesions_or_mutations_mark_strand(self, selection_coefs):
        out1 = {}
        out2 = {}
        out_both = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            mut_strand1, mut_strand2, mut_both = chrom_pair.get_driver_lesions_or_mutations_mark_strand(selection_coefs[chrom_pair_ind])
            if len(mut_strand1) > 0:
                out1[chrom_pair_ind] = mut_strand1
            if len(mut_strand2) > 0:
                out2[chrom_pair_ind] = mut_strand2
            if len(mut_both) > 0:
                out_both[chrom_pair_ind] = mut_both
        return out1, out2, out_both

    def get_driver_lesions(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            lesions = chrom_pair.get_driver_lesions_fitness(selection_coefs[chrom_pair_ind])
            if len(lesions) > 0:
                out[chrom_pair_ind] = lesions
        return out

    def get_driver_lesions_or_mutations(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            positions = chrom_pair.get_driver_lesions_or_mutations(selection_coefs[chrom_pair_ind])
            if len(positions) > 0:
                out[chrom_pair_ind] = positions
        return out

    def get_two_strand_driver_lesions_or_mutations(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            muts = chrom_pair.get_two_strand_driver_lesions_or_mutations(selection_coefs[chrom_pair_ind])
            if len(muts) > 0:
                out[chrom_pair_ind] = muts
        return out

    def get_mutations(self):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            muts = chrom_pair.get_mutations()
            if len(muts) > 0:
                out[chrom_pair_ind] = muts
        return out

    def get_driver_info(self, selection_coefs):
        out = {}
        for chrom_pair_ind in range(len(self.chromosomes)):
            chrom_pair = self.chromosomes[chrom_pair_ind]
            info = chrom_pair.get_driver_info(selection_coefs[chrom_pair_ind])
            out[chrom_pair_ind] = info
        return out

    def get_fitness(self, selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, smax):
        '''
        Gets fitness for each chromosome, computes overall fitness of cell accordingly
        :return:
        '''
        def helper(to_iterate, epistasis_type, selection_coefs, n_drivers, skip=0):
            to_return = 1
            antagonostic_found = False
            for chrom, chrom_drivers in to_iterate.items():
                for mut in chrom_drivers:
                    s = selection_coefs[chrom][mut]
                    error_free = binomial_rv(1, skip)
                    if error_free != 1:
                        if epistasis_type == 'NO EPISTASIS':
                            to_return *= (1 + s)
                        elif epistasis_type == 'ANTAGONISTIC':
                            if not antagonostic_found:
                                to_return *= (1 + s)
                                antagonostic_found = True
                        elif epistasis_type == 'SYNERGESTIC':
                            if n_drivers > 1:
                                to_return *= (1 + s)
            return to_return

        def process_drivers_tx_strand(strand1_drivers, strand2_drivers, both, driver_strands):
            # filter down to tx strand lesions or mutations
            updated_strand1 = {}
            updated_strand2 = {}
            for chr, pos_list in strand1_drivers.items():
                for pos in pos_list:
                    if driver_strands[chr][pos] == 'tx':
                        if chr in updated_strand1:
                            updated_strand1[chr].add(pos)
                        else:
                            updated_strand1[chr] = {pos}
            for chr, pos_list in strand2_drivers.items():
                for pos in pos_list:
                    if driver_strands[chr][pos] == 'tx':
                        if chr in updated_strand2:
                            updated_strand2[chr].add(pos)
                        else:
                            updated_strand2[chr] = {pos}
            drivers_out = both
            drivers_out.update(updated_strand1)
            drivers_out.update(updated_strand2)
            return drivers_out

        if require_two_strands:
            if error_free_tx == 1: # lesions don't count
                drivers = self.get_two_strand_driver_mutations(selection_coefs)
            elif error_free_tx == 0: # lesions count
                drivers = self.get_two_strand_driver_lesions_or_mutations(selection_coefs)
            else:
                raise ValueError("error free tx must be 0 or 1")
        else:
            if error_free_tx == 1: # lesions don't count
                drivers1, drivers2, both = self.get_driver_mutations_mark_strand(selection_coefs)
                drivers = process_drivers_tx_strand(drivers1, drivers2, both, driver_strands)
            elif error_free_tx == 0: # lesions count
                drivers1, drivers2, both = self.get_driver_lesions_or_mutations_mark_strand(selection_coefs)
                # print(drivers1, drivers2, both)
                drivers = process_drivers_tx_strand(drivers1, drivers2, both, driver_strands)
                # print(drivers)
        n_drivers = 0
        for chrom, pos in drivers.items():
            n_drivers += len(pos)
        fitness = helper(drivers, epistasis_type, selection_coefs, n_drivers)
        s = (fitness-1)/smax
        fitness = 1+s
        return fitness

    def get_branching_process_probabilities(self, selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate_fit, error_free_tx, smax):

        fitness = self.get_fitness(selection_coefs, driver_strands, epistasis_type, require_two_strands, error_free_tx, smax)

        s = fitness-1
        c = asymmetric_rate_0 if s == 0 else asymmetric_rate_fit
        b = (s - c + 1)/2
        d = 1 - b - c
        # print(b, d, c)
        if b+c > 1:
            raise ValueError("provided asymmetric rate too high for driver fitness")
        return b, d, c


    def division_decision(self, selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate_fit, error_free_tx, smax):
        '''
        w = 1+log(b/d)
        b+d = 1-asymmetric_rate
        :return: current options: 'death', 'division'
        '''
        birth_rate, death_rate, c_rate = self.get_branching_process_probabilities(selection_coefs, driver_strands, epistasis_type, require_two_strands, asymmetric_rate_0, asymmetric_rate_fit, error_free_tx, smax)
        choice = np.random.choice(a=['asymmetric', 'birth', 'death'],
                                  replace=True,
                                  p=[c_rate, birth_rate, death_rate])
        return choice



    def has_driver_lesions(self, selection_coefs):
        for chrom_pair_ind in range(len(self.chromosomes)):
            if self.chromosomes[chrom_pair_ind].has_driver_lesions(selection_coefs[chrom_pair_ind]):
                return True
        return False

    def has_lesions(self):
        for chrom_pair in self.chromosomes:
            if chrom_pair.has_lesions():
                return True
        return False

    def divide(self, error_free_rate, propagate_id=False):
        '''
        :param error_free_rate:
        :param propagate_id:    if False, daughter cells are given same id as parent
                                if True, daughter cells are given suffix A or B from parent id
        :return:
        '''
        daughter_1 = []
        daughter_2 = []
        daughter_1_id = self.id+'A' if propagate_id else self.id
        daughter_2_id = self.id+'B' if propagate_id else self.id
        for chrom_pair in self.chromosomes:
            chrom_pair_1, chrom_pair_2 = chrom_pair.divide(error_free_rate)
            coin_flip = binomial_rv(1, 0.5)
            if coin_flip == 0:
                daughter_1.append(chrom_pair_1), daughter_2.append(chrom_pair_2)
            else:
                daughter_1.append(chrom_pair_2), daughter_2.append(chrom_pair_1)
        return Cell(self.chr_size_list, daughter_1_id, daughter_1), Cell(self.chr_size_list, daughter_2_id, daughter_2)

