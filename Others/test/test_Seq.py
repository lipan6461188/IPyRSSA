#########
#########   Test Seq.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import Seq

#####################
#  reverse_comp(sequence)
#####################

print Seq.reverse_comp("TAGCTAGCTGGTTAGTTCTATC")
print Seq.reverse_comp("TAGCTAatgcatTAGTTCTATC")
print Seq.reverse_comp("TAGCTAGCT---TAGTTC--TC")
print Seq.reverse_comp("TANNNNNNNNGTTAGTTCTATC")

#####################
#  flat_seq(sequence, lineLen=60)
#####################
import General
seqFn = "test_seq.fasta"
fasta = General.load_fasta(seqFn, rem_tVersion=False)

print Seq.flat_seq(fasta['ENST00000580210.5'])
print Seq.flat_seq(fasta['ENST00000580210.5'], lineLen=10)
print Seq.flat_seq(fasta['ENST00000580210.5'], lineLen=100)
print Seq.flat_seq("ACAGATTGTT")

#####################
#  format_gene_type(gene_type)
#####################
raw_gene_types = ['TEC', 'snRNA', 'processed_transcript', 'protein_coding', 'retained_intron', 'miRNA', 'antisense', 'processed_pseudogene', 'C_region', 'telomerase_RNA', 'tRNA', 'misc_RNA', 'lnc_RNA']

for raw in raw_gene_types:
    print raw+" ==> "+Seq.format_gene_type(raw)

#####################
#  seqClass(object)
#####################

seqFn = "test_seq.fasta"
seqFetcher = Seq.seqClass(seqFn)
seqFetcher.fetch('ENST00000580210.5', 0, 10, '+')
seqFetcher.fetch('ENST00000580210.5', 1, 10, '-')

#mm10:
#chr1    4938576 4940710 -       Gm37079=ENSMUSG00000102653.1    ENSMUST00000194114.1    TEC     4938576-4940710
#chr1    43742915        43743526        +       Gm29157=ENSMUSG00000100635.1    ENSMUST00000188753.1    lincRNA 43742915-43743526

seqFn = "/150T/zhangqf/GenomeAnnotation/genome/mm10.fa"
seqFetcher = Seq.seqClass(seqFn)
print seqFetcher.fetch('chr1', 4938576-1, 4940710, '-')
print seqFetcher.fetch('chr1', 43742915-1, 43743526, '+')

print seqFetcher.has("chr1")
print seqFetcher.seq_len_dict()

#####################
#  lift_genome(lifter, chrID, chrStart, chrEnd, chrStrand, verbose=False)
#####################

from pyliftover import LiftOver
lifter = LiftOver("mm9", "mm10")
print Seq.lift_genome(lifter, "chr12", 27163765, 27163786, "+", verbose=False)
print Seq.lift_genome(lifter, "chrXX", 27163765, 27163786, "+", verbose=False)
print Seq.lift_genome(lifter, "chrXX", 27163765, 27163786, "+", verbose=True)

