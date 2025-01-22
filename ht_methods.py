from ht_adts import HashTableAppend
from random import SystemRandom
from Bio import SeqIO


nt = ['A', 'U', 'G', 'C']  # going with U as uracil for this, alternatively, use T for compatibility
# nt = ['A', 'T', 'G', 'C']
bases = [["A", "U"], ["G", "C"]]
# bases = [["A", "T"], ["G", "C"]]
    

def all_seq(N: int, seq_elements: list[any] =nt, starting_sequences: list[str] =[]) -> list[str]:
    """returns a list of all possible sequences of sequence elements (=nt) of length N."""
    if starting_sequences == []:
        seq_list = ['']
    else:
        seq_list = starting_sequences
    while len(seq_list[0]) < N:
        temp_list = []
        for s in seq_list:
            for e in seq_elements:
                temp_list.append(s + e)
        seq_list = temp_list
    return seq_list


def get_mirbase_seeds(n: int=6, filename: str ='data_sources/mature.fa', species_id: str ='hsa'):  # |''
    """use an empty string for species_id to get a dictionary with all {species : list of seeds}"""
    species_dict = HashTableAppend() # {species id : list[str]}
    with open(filename, 'r') as mature:
        for line in mature:
            # >hsa-...
            species = line[1:4] if line.startswith('>') else species
            if not line.startswith('>'):  # 0 123456 7 8... (6nt seed)
                if species == '' or species == species_id:
                    species_dict.add(species, line[1:n + 1])
    return species_dict.get(species_id) if species_id != '' else species_dict


# def get_RNA_activation_sequence_regions():
#     pass


def get_transcripts(filename: str ='data_sources/06-18-2024-gencode.v46.transcripts.fa'):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    transcripts = []
    max_length = 0
    sum_length = 0
    for fasta in fasta_sequences:
        transcripts.append(str(fasta.seq).replace('T', 'U'))  # make them RNA sequences
        if len(str(fasta.seq)) > max_length:
            max_length = len(str(fasta.seq))
        sum_length += len(str(fasta.seq))
        # break # for testing w one transcript sequence
    print("Maximum transcript length: " + str(max_length))
    print("Average transcript length: " + str(sum_length/len(transcripts)))
    print(f"Number of transcripts: {len(transcripts)}")
    return transcripts


# processes GTF file
def _get_three_prime_utr_lengths(filename: str ='data_sources/06-18-2024-gencode.v46.annotation.gtf'):
    """returns a list of integers, either -1: indicating no 3'-UTR, or an integer = length of 3'-UTR"""
    read_line = ""
    current_transcript_bounds = [-1, -1]
    current_stop_codon_bounds = [-1, -1]
    current_utr_bounds = [-1, -1]
    indexed_utr_lengths = []
    gtf_file = open(filename, 'r')
    for i in range(5):
        gtf_file.readline()  # This advances file read to correct row to apply following logic
    while gtf_file.readable():
        read_line = gtf_file.readline().strip()
        split_line = read_line.split("\t")
        if read_line == "":
            break
        else:  #                             ### by comparison to UTRsequences from TargetScan, this works
            if split_line[2] == "transcript":
                indexed_utr_lengths.append(-1)  # Later, this will be interpreted as no 3'-UTR
                current_transcript_bounds[0] = int(split_line[3])
                current_transcript_bounds[1] = int(split_line[4])
            elif split_line[2] == "stop_codon":
                current_stop_codon_bounds[0] = int(split_line[3])
                current_stop_codon_bounds[1] = int(split_line[4])
            elif split_line[2] == "UTR":  # the sequences making up a 3'-UTR are listed separately, and must be combined
                this_utr_has_stop_codon = False
                this_utr_has_transcript_bound = False
                current_utr_bounds[0] = int(split_line[3])
                current_utr_bounds[1] = int(split_line[4])
                if current_utr_bounds[0] == current_stop_codon_bounds[0] or current_utr_bounds[0] == current_stop_codon_bounds[1]:
                    this_utr_has_stop_codon = True
                if current_utr_bounds[1] == current_stop_codon_bounds[0] or current_utr_bounds[1] == current_stop_codon_bounds[1]:
                    this_utr_has_stop_codon = True
                if this_utr_has_stop_codon:
                    if current_utr_bounds[0] == current_transcript_bounds[0] or current_utr_bounds[0] == current_transcript_bounds[1]:
                        this_utr_has_transcript_bound = True
                    if current_utr_bounds[1] == current_transcript_bounds[0] or current_utr_bounds[1] == current_transcript_bounds[1]:
                        this_utr_has_transcript_bound = True
                if this_utr_has_stop_codon and this_utr_has_transcript_bound:
                    indexed_utr_lengths[-1] = abs(current_utr_bounds[0] - current_utr_bounds[1]) + 1

    return indexed_utr_lengths


def get_three_prime_utrs(fasta_file: str='data_sources/06-18-2024-gencode.v46.transcripts.fa',
                         gtf_file: str='data_sources/06-18-2024-gencode.v46.annotation.gtf'):
    """processes .gtf file for 3' lengths, and .fa file for transcript sequences"""
    three_prime_lengths = _get_three_prime_utr_lengths(gtf_file)
    three_primes = []
    max_length = 0
    sum_length = 0
    avg_divisor = 0
    three_prime_index = 0
    for fasta in SeqIO.parse(open(fasta_file),'fasta'):
        three_primes.append('')
        if three_prime_lengths[three_prime_index] > 3:  # ignores 3'-UTRs that are only a stop codon
            string = str(fasta.seq)
            s_len = len(string)
            string = string[s_len - three_prime_lengths[three_prime_index] + 3:]  # so that the string does not include stop codon
            if s_len > max_length:
                max_length = s_len
            sum_length += s_len
            three_primes[-1] = (str(string).replace('T', 'U'))  # make them RNA sequences, assign to otherwise ''
            avg_divisor += 1
        three_prime_index += 1
        # break # for testing w one transcript sequence
    print("Maximum 3'-UTR length: " + str(max_length))
    print("Average 3'-UTR length: " + str(sum_length/avg_divisor))
    print(f"Number of 3'-UTRs: {avg_divisor}")
    return three_primes, three_prime_lengths


def get_random_sequence_lists(random_1_0_file: str="data_sources/download_random-org_06-23-2024", transcripts: list[str] =[''], utr3s: list[str] =['']):
    # floor(average 3'-UTR length: 1214.313) = 1214
    # random bits downloaded as a text file with a single string of ASCII 1's and 0's, from RANDOM.ORG, 06-23-2024 (couldn't
    # find a good website for darkness of space / (quantum?) noise possibly re Bootes Void- but maybe I was remembering this
    
    random_bit_string = open(random_1_0_file, 'r').read()
    
    random_sequences_100 = []
    k = 0
    for _ in range(100):
        the_random_sequence = ""
        sequence_length = 1214
        for j in range(sequence_length):
            the_random_sequence += bases[
                int(random_bit_string[k % len(random_bit_string)])][
                int(random_bit_string[(k + 1) % len(random_bit_string)])]
            k = k + 2
        random_sequences_100.append(the_random_sequence)
    
    random_sequences_1000 = []  # lengths: |each utr3|
    k = 0
    for u in utr3s:
        if u != '':
            the_random_sequence = ""
            for j in range(len(u)):
                the_random_sequence += bases[
                    int(random_bit_string[k % len(random_bit_string)])][
                    int(random_bit_string[(k + 1) % len(random_bit_string)])]
                k = k + 2
            random_sequences_1000.append(the_random_sequence)
    
    random_sequences_10k = []  # lengths: |each transcript|
    k = 0
    for t in transcripts:
        the_random_sequence = ""
        for _ in t:
            the_random_sequence += bases[
                int(random_bit_string[k % len(random_bit_string)])][
                int(random_bit_string[(k + 1) % len(random_bit_string)])]
            k = k + 2
        random_sequences_10k.append(the_random_sequence)
    
    return [random_sequences_100, random_sequences_1000, random_sequences_10k]

# def pickle_all_current_objects  TODO



class RegSeqGenerator:
    """this does not yield any results, but returns them."""
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


