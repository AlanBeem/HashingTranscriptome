from random import SystemRandom


class RegSeqGenerator:
    """this does not yield any results, """
    def __init__(self, master_regulatory_sequence_length, regulatory_letters_list):
        self.master_reg_seq = ""
        self.sequence_elements = regulatory_letters_list.copy()
        self.regulatory_sequences = []
        for rsg_i in range(master_regulatory_sequence_length):
            self.master_reg_seq += str(
                self.sequence_elements[SystemRandom().randrange(0, len(self.sequence_elements), 1)])
        self.walk_index = 0

    def get_random_reg_seq(self, reg_seq_length: int):
        self.regulatory_sequences.append("".join([self.sequence_elements[SystemRandom().randrange(0, len(self.sequence_elements))] for rsg_i in range(reg_seq_length)]))
        return self.regulatory_sequences[-1]

    def get_subseq_of_master_reg_seq(self, reg_seq_length: int, index_of_starting_nucleotide: int):
        """this treats the master regulatory sequence indices as circular and returns a subsequence|length ==
        reg_seq_length"""
        self.regulatory_sequences.append("".join([self.master_reg_seq[rsg_j % len(self.master_reg_seq)] for rsg_j in
                        range(index_of_starting_nucleotide, index_of_starting_nucleotide + reg_seq_length + 1, 1)]))
        return self.regulatory_sequences[-1]

    def get_walked_random_reg_seq(self, reg_seq_length: int, random_walk_start: float, random_walk_end: float,
                                  sticky_index: float):
        """this treats the master regulatory sequence indices as circular and returns a sequence composed of random
        sequence elements and nucleotides from random walks of self.master_reg_seq: str"""
        gwrrs_walking = False
        gwrrs_walk_index = 0
        gwrrs_sequence = ""
        while len(gwrrs_sequence) < reg_seq_length:
            if not gwrrs_walking:
                if SystemRandom().random() < random_walk_start:
                    gwrrs_walking = True
                    if SystemRandom().random() < sticky_index:
                        gwrrs_walk_index = self.walk_index
                    else:
                        gwrrs_walk_index = SystemRandom().randrange(len(self.master_reg_seq))
            else:
                if SystemRandom().random() < random_walk_end:
                    gwrrs_walking = False
                    self.walk_index = gwrrs_walk_index
            if gwrrs_walking:
                gwrrs_sequence += self.master_reg_seq[gwrrs_walk_index % len(self.master_reg_seq)]
                gwrrs_walk_index += 1
            else:
                gwrrs_sequence += self.sequence_elements[SystemRandom().randrange(0, len(self.sequence_elements))]
        self.regulatory_sequences.append(gwrrs_sequence)
        return self.regulatory_sequences[-1]