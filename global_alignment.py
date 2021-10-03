import _global_alignment as gal
import numpy as np 

def global_alignment(subject, query, match=10, mis_match=-1, gap_extension=-1, gap_open=-20): 
    result = gal.align(subject.encode('utf-8'),query.encode('utf-8'), match, mis_match, gap_extension, gap_open, 1)
    a, b   = "".join(map(chr, result[3])), "".join(map(chr, result[4]))
    return a,b 

if __name__ == "__main__":
    import Bio
    import collections
    from tqdm import tqdm
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    def read_fastq(fastq_name):
        #Read fastq file
        seq_dict = collections.defaultdict(dict)
        with open(fastq_name.replace("'","").replace("\\","")) as f:
            n = 0 
            for line in tqdm(f, total=224892):
                if line[0] == "@" and n % 4 == 0:
                    key = line[1:].rstrip() 
                    key = key.split(" ")[0] 
                    key = key.replace(":","_")
                    seq_dict[key]["key"] = line[1:].rstrip() 
                elif n%4 == 1:
                    seq_dict[key]["seq"] = line.rstrip() 
                elif n%4 == 2:
                    seq_dict[key]["option"] = line.rstrip() 
                elif n%4 == 3:
                    seq_dict[key]["quality"] = line.rstrip() 
                n += 1
        return seq_dict 

    path = "../../sample_data/Target-ACE_sequence_data/Rep1/AID/polyC/chr19_8211972/"
    ref  = "TAACTTACGGAGTCGCTCTACGCTAAGTCCCCGTGTAAACAGAGCTGAACCTGCAGGCAGGTAAGAGTGTCCCCGGCCTGTGCCCCCCCACCTCCAGACGGCGGTAGCACTCACGTACGACACACTGGCTCCTGCCAGGCAGGACCTAAAGAATCCCATCC" 
    fastq_dict = read_fastq(path + "/R1.fastq")
    
    with open("result.txt", "w") as o:
        for key in tqdm(fastq_dict):
            query = fastq_dict[key]["seq"] 
            a, b  = global_alignment(ref, query, match=10, mis_match=-1, gap_extension=-1, gap_open=-20) 
            #print(a)
            #print(b) 
            #print(">", file=o)
            #print("reference_sequence", a, sep=",", file=o)
            #print(key, b, sep=",", file=o)

    """
    x = "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG"
    y = "CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG"
    #x = "ATGAGTCTCTCTGATAAGGACAAGGCTG"
    #y = "CTGTCTCCTGCCGACAAGACCAACGTCA"
    #x = "AAAA"
    #y = "TTTAAAATTT"
    #x = "ATTTTTGC"
    #y = "ATTTTGC"
    #y = "ACACT"
    #x = "AAT" 
    a, b = global_alignment(x, y, match=1, mis_match=-1, gap_extension=-1, gap_open=-3) 
    print(a)
    print(b)
    print()

    for a in pairwise2.align.globalms(x, y, 1, -1, -3, -1, force_generic=True, penalize_extend_when_opening=True, penalize_end_gaps=True):
        print(format_alignment(*a))
    """
