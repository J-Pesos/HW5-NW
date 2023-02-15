# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    print('Main v6')
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    # Initialize NeedlemanWunsch class.
    NW = NeedlemanWunsch(sub_matrix_file = './substitution_matrices/BLOSUM62.mat', gap_open = -10, gap_extend = -1)
    
    # Calculate all alignment score possibilities against humans.
    gg_score, hs_gg, gg = NW.align(hs_seq, gg_seq)
    mm_score, hs_mm, mm = NW.align(hs_seq, mm_seq)
    br_score, hs_br, br = NW.align(hs_seq, br_seq)
    tt_score, hs_tt, tt = NW.align(hs_seq, tt_seq)

    # Concatenate scores into a list.
    aligned_scores = [(gg_score, 'Gallus gallus'), (mm_score, 'Mus musculus'), (br_score, 'Balaeniceps rex'), (tt_score, 'tursiops truncatus')]
    print(aligned_scores)
    # Sort list of scores.
    sorted_scores = sorted(aligned_scores, key = lambda x: x[1], reverse = True)
    # Print out species in order of the highest alignment score.
    print('List of species in order of most similar to human BRD2 sequence.')
    print([score[1] for score in sorted_scores])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print('List of alignment scores and what sequence was aligned to human BRD2 sequence.')
    print( [ (score[0], score[1]) for score in sorted_scores ] )
    
if __name__ == "__main__":
    main()
