from collections import namedtuple
import pandas as pd
import pickle

Position = namedtuple('Position', ['genome', 'chrom', 'pos'])

class PositionedSNP:
    '''
    Represent a positioned SNP: two alleles and several positions in the genomes where this SNP appears
    '''
    def __init__(self, allele_1, allele_2):
        assert allele_1 < allele_2, "Alleles should be given in a canonical order"
        self.alleles = [allele_1, allele_2]
        self.positions = set()

    def add_pos (self, position):
        self.positions.add(position)

    def merge(self, other):
        '''
        Creates a new PositionedSNP with the same allele and positions merged
        :param other: other PositionedSNP
        :return: new PositionedSNP with the same allele and positions merged
        '''
        assert self.alleles == other.alleles # we should not merge SNPs that does not have the same alleles
        new_PositionedSNP = PositionedSNP(*self.alleles)
        new_PositionedSNP.positions = self.positions.union(other.positions)
        return new_PositionedSNP


PositionedSNPIndex = namedtuple('PositionedSNPIndex', ['position', 'allele_1', 'allele_2'])
class PositionedSNPs:
    def __init__(self):
        self.PositionedSNPIndex_to_PositionedSNP = {}  # index the positioned SNPs

    def add_SNPs_from_csv(self, csv_file, genome_1, genome_2):
        '''
        :param csv_file: a csv file with the SNPs computed by get_SNPs_using_mummer rule
        :param genome_1: a string with genome_1 name
        :param genome_2: a string with genome_2 name
        '''
        snps_dataframe = pd.read_csv(csv_file, sep = "\t")

        # populate self.PositionedSNPIndex_to_PositionedSNP
        for i, row in snps_dataframe.iterrows():
            # get the data
            position_1 = Position(genome=genome_1, chrom=row["ref_chrom"], pos=row["ref_pos"])
            position_2 = Position(genome=genome_2, chrom=row["query_chrom"], pos=row["query_pos"])
            allele_1 = min(row["ref_sub"], row["query_sub"])
            allele_2 = max(row["ref_sub"], row["query_sub"])
            assert allele_1 != allele_2, "Allele SNPs have the same base"
            PositionedSNPIndex_1 = PositionedSNPIndex(position=position_1, allele_1=allele_1, allele_2=allele_2)
            PositionedSNPIndex_2 = PositionedSNPIndex(position=position_2, allele_1=allele_1, allele_2=allele_2)

            # get the previous positioned SNPs, if any
            previous_PositionedSNP_1 = self.PositionedSNPIndex_to_PositionedSNP.get(PositionedSNPIndex_1)
            previous_PositionedSNP_2 = self.PositionedSNPIndex_to_PositionedSNP.get(PositionedSNPIndex_2)

            '''
            # verifies if previous_PositionedSNP_1 and previous_PositionedSNP_2 are the same
            previous_PositionedSNPs_are_the_same = \
                previous_PositionedSNP_1 is previous_PositionedSNP_2 or \
                previous_PositionedSNP_1 is None or \
                previous_PositionedSNP_2 is None

            # making sure they correspond to the same object
            if not previous_PositionedSNPs_are_the_same:
                import pdb; pdb.set_trace()
            assert previous_PositionedSNPs_are_the_same, "Bug: previous_PositionedSNPs are not the same"
            '''

            if previous_PositionedSNP_1 is None and previous_PositionedSNP_2 is None:
                # we need to create a current_PositionedSNP
                current_PositionedSNP = PositionedSNP(allele_1, allele_2)

                # and associate it to these positions in the index
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_1] = current_PositionedSNP
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_2] = current_PositionedSNP
            elif previous_PositionedSNP_1 is None:
                # update self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_1]
                current_PositionedSNP = previous_PositionedSNP_2
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_1] = current_PositionedSNP
            elif previous_PositionedSNP_2 is None:
                # update self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_2]
                current_PositionedSNP = previous_PositionedSNP_1
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_2] = current_PositionedSNP
            elif not(previous_PositionedSNP_1 is previous_PositionedSNP_2):
                # if both SNPs do not point to the same PositionedSNP, it means they have to be merged
                current_PositionedSNP = previous_PositionedSNP_1.merge(previous_PositionedSNP_2)

                # we associate these positions in the index to the new merged PositionedSNP
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_1] = current_PositionedSNP
                self.PositionedSNPIndex_to_PositionedSNP[PositionedSNPIndex_2] = current_PositionedSNP

                # and also all the previous positions now point to this PositionedSNP
                for position in current_PositionedSNP.positions:
                    self.PositionedSNPIndex_to_PositionedSNP[position] = current_PositionedSNP
            else:
                assert previous_PositionedSNP_1 is previous_PositionedSNP_2 and previous_PositionedSNP_1 is not None and previous_PositionedSNP_2 is not None
                current_PositionedSNP = previous_PositionedSNP_1


            # add the positions to current_PositionedSNP
            current_PositionedSNP.add_pos(position_1)
            current_PositionedSNP.add_pos(position_2)

    def get_nb_SNPs_in_pangenome(self):
        return len(set(self.PositionedSNPIndex_to_PositionedSNP.values()))

    def get_nb_SNPs_that_can_be_found_with_a_given_genome(self, genome):
        all_positioned_SNPs = set(self.PositionedSNPIndex_to_PositionedSNP.values())
        nb = 0
        for positioned_SNP in all_positioned_SNPs:
            for position in positioned_SNP.positions:
                if position.genome == genome:
                    nb += 1
                    break
        return nb

    @classmethod
    def load (cls, filename):
        '''
        :param filename: contains the json to load
        '''
        with open(filename, "rb") as fin:
            positionedSNPs = pickle.load(fin)
        return positionedSNPs

    def serialize(self, filename):
        with open(filename, "wb") as fout:
            pickle.dump(self, fout)





# TODO: add unit test
'''
positionedSNPs = PositionedSNPs()
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.2.mummer.csv', 'genome.1', 'genome.2')
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.1/genome.1-SEP-genome.3.mummer.csv', 'genome.1', 'genome.3')
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.2/genome.2-SEP-genome.3.mummer.csv', 'genome.2', 'genome.3')
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.1.mummer.csv', 'genome.0', 'genome.1')
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.2.mummer.csv', 'genome.0', 'genome.2')
positionedSNPs.add_SNPs_from_csv('assemblies_sample_out/genome.0/genome.0-SEP-genome.3.mummer.csv', 'genome.0', 'genome.3')
positionedSNPs.serialize("assemblies_sample_out/positionedSNPs")
positionedSNPs_loaded = positionedSNPs.load("assemblies_sample_out/positionedSNPs")
'''