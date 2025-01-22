from ht_adts import HashTableIncrement
from ht_methods import get_transcripts, get_three_prime_utrs, get_random_sequence_lists, get_mirbase_seeds, all_seq, RegSeqGenerator
import csv


n_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
count_folder = 'count_data'
# This writes files if run as main

def get_dicts_from_files(count_filepath: str ='count_data'):
    # , loc_filepath: str ='', related to a multi-day runtime approach, pretty much only good for Jaccard Similarity of targets,
    # it was a useful experiment-
    
    # setup
    labels = ['transcripts', 'utr3s', 'r100', 'r1k', 'r10k']  #, 'seq_gen_1']  RegSeqGenerator to be used in ipynb, may write to this location

    # for counts {seq: int}  # writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)}) 
    a_count_data = [dict() for _ in n_list]
    b_count_data = [dict() for _ in n_list]
    c_count_data = [[dict() for _ in n_list] for _ in ['r100', 'r1k', 'r10k']]
    count_data = [a_count_data, b_count_data, c_count_data[0], c_count_data[1], c_count_data[2]]  # corresponding to labels

    # read count data into dictionaries
    for label in labels:
        for i in range(len(n_list)):
            open_file = open(f"{count_filepath}/{f"random/{label}" if label.startswith('r') else label}/{n_list[i]}nt_subseq_{label}.csv", 'r')
            for line in open_file.readlines():
                if line == '\n':
                    break
                line = line.split(',')
                count_data[labels.index(label)][i].update({line[0]:int(line[1].strip('\n'))})
            open_file.close()

    return labels, count_data


if __name__ == '__main__':
   ## Data sources, also writes down txt files 
    # (a)
    transcripts = get_transcripts()
    with open('count_data/transcripts.txt', 'w') as utr3s_file:
        for t in transcripts:
            utr3s_file.write(t + '\n')
    print("transcripts.txt written")
    # (b)
    three_prime_utrs, three_prime_utrs_lengths = get_three_prime_utrs()  # length == len(transcripts), '' indicates no 3'-UTR
    with open('count_data/utr3s_length_sequence.txt', 'w') as utr3s_file:
        for i in range(len(three_prime_utrs)):
            utr3s_file.write(str(three_prime_utrs_lengths[i]) + ',' + three_prime_utrs[i] + '\n')
    print("utr3s_length_sequence.txt written")
    # (c)
    rand_seqs_list = get_random_sequence_lists(transcripts=transcripts, utr3s=three_prime_utrs) # these are matched to (a) (b)
    with open('count_data/random/r100.txt', 'w') as r100_file:
        for r in rand_seqs_list[0]:
            r100_file.write(str(len(r)) + ',' + r + '\n')
    print("r100.txt written")
    with open('count_data/random/r1k.txt', 'w') as r1k_file:
        for r in rand_seqs_list[1]:
            r1k_file.write(str(len(r)) + ',' + r + '\n')
    print("r1k.txt written")
    with open('count_data/random/r10k.txt', 'w') as r10k_file:
        for r in rand_seqs_list[2]:
            r10k_file.write(str(len(r)) + ',' + r + '\n')
    print("r10k.txt written")
    # (f)
    # hsa_mirbase_seeds = get_mirbase_seeds()
   ##
    print("Data loaded")
    #
   ## Apply Algorithms 1, write to files;
    # (a), (b)
    all_transcripts_n_data_increment = [HashTableIncrement() for _ in n_list]  # (sequence, occurences)
    utr3_n_data_increment = [HashTableIncrement() for _ in n_list]  # (sequence, occurences)
    # in principle, this could be accomplished by a single dictionary, but, this is employed for writing them down separately
    for j in range(len(transcripts)):
        t_len = len(transcripts[j])
        for i in range(t_len):
            for n in n_list:  # 1 through 12; so n - 1 = index of n_data
                if i + n <= t_len:
                    all_transcripts_n_data_increment[n - 1].add(transcripts[j][i: i + n])  # slices of str
            if three_prime_utrs[j] != '' and i >= t_len - three_prime_utrs_lengths[j]:
                for n in n_list:
                    if i + n <= t_len:
                        utr3_n_data_increment[n - 1].add(transcripts[j][i:i + n])
    print("Transcripts and 3'-UTR's processed")
    all_sequences = []  # it's not the most efficient, but this writes down a 0 where applicable
    # # all sequences for N=1: all_seq(1) -> ['A', 'U', 'G', 'C']
    for i in range(len(n_list)):
        # for each in n_data, write a csv
        n_d = all_transcripts_n_data_increment[i]
        n = n_list[i]
        all_sequences.append(all_seq(n_list[i], starting_sequences=all_sequences[-1] if len(all_sequences) > 0 else []))
        with open(f"{count_folder}/transcripts/{n}nt_subseq_transcripts.csv", 'w') as csv_file:  # {data_output_location}/
            writer = csv.DictWriter(csv_file, fieldnames=['Sequence', 'Occurences'])
            for s in all_sequences[-1]:
                writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)})
    print("Transcripts counts written")
    for i in range(len(n_list)):
        # for each in n_data, write a csv
        n_d = utr3_n_data_increment[i]
        seq = all_sequences[i]
        n = n_list[i]
        with open(f"{count_folder}/utr3s/{n}nt_subseq_utr3s.csv", 'w') as csv_file:  # {data_output_location}/
            writer = csv.DictWriter(csv_file, fieldnames=['Sequence', 'Occurences'])
            for s in seq:
                writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)})
    print("3'-UTR's counts written")
    # (c) (1)
    r100_n_data_increment = [HashTableIncrement() for _ in n_list]  # (sequence, occurences)
    id = 0
    for t in rand_seqs_list[0]:  # zip(range(len(rand_seqs_list[0])), 
        t_len = len(t)
        for i in range(t_len):
            for n in n_list:  # eg 1 through 12; so -1 -> index in n_data..
                if i + n <= t_len:
                    r100_n_data_increment[n - 1].add(t[i: i + n])
        id += 1
    print("r100 processed")
    for i in range(len(n_list)):
        # for each in n_data, write a csv
        n_d = r100_n_data_increment[i]
        seq = all_sequences[i]
        n = n_list[i]
        with open(f"{count_folder}/random/r100/{n}nt_subseq_r100.csv", 'w') as csv_file:  # {data_output_location}/
            writer = csv.DictWriter(csv_file, fieldnames=['Sequence', 'Occurences'])
            for s in seq:
                writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)})
    print("r100 counts written")
    # (c) (2)
    # r1k is really UTR lengths
    r1k_n_data_increment = [HashTableIncrement() for _ in n_list]  # (sequence, occurences)
    id = 0
    for t in rand_seqs_list[1]:
        t_len = len(t)
        for i in range(t_len):
            for n in n_list:  # eg 1 through 12; so -1 -> index in n_data..
                if i + n <= t_len:
                    r1k_n_data_increment[n - 1].add(t[i: i + n])
        id += 1
    print("r1k processed")
    for i in range(len(n_list)):
        # for each in n_data, write a csv
        n_d = r1k_n_data_increment[i]
        seq = all_sequences[i]
        n = n_list[i]
        with open(f"{count_folder}/random/r1k/{n}nt_subseq_r1k.csv", 'w') as csv_file:  # {data_output_location}/
            writer = csv.DictWriter(csv_file, fieldnames=['Sequence', 'Occurences'])
            for s in seq:
                writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)})
    print("r1k counts written")
    # (c) (3)
    # r10k is really transcripts average length
    r10k_n_data_increment = [HashTableIncrement() for _ in n_list]  # (sequence, occurences)
    id = 0
    for t in rand_seqs_list[2]:
        t_len = len(t)
        for i in range(t_len):
            for n in n_list:  # eg 1 through 12; so -1 -> index in n_data..
                if i + n <= t_len:
                    r10k_n_data_increment[n - 1].add(t[i: i + n])
        id += 1
    print("r10k processed")
    for i in range(len(n_list)):
        # for each in n_data, write a csv
        n_d = r10k_n_data_increment[i]
        seq = all_sequences[i]
        n = n_list[i]
        with open(f"{count_folder}/random/r10k/{n}nt_subseq_r10k.csv", 'w') as csv_file:  # {data_output_location}/
            writer = csv.DictWriter(csv_file, fieldnames=['Sequence', 'Occurences'])
            for s in seq:
                writer.writerow({'Sequence':s, 'Occurences':n_d.get(s, 0)})
    print("r10k counts written")

    print("Count data files written")
# seq_gen_1 - try to tweak RegSeqGenerator to produce sequences such that Alg.3(seq_gen_1) -> looks like Alg.3(3'-UTR) -> ...
### Looks like this use of the random data might not be random enough, as occurring 12nt sequences seem non-uniform, or non-randomly distributed

