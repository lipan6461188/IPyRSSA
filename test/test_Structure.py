#########
#########   Test Cluster.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import Structure

import General
seqFn = "test_seq.fasta"
fasta = General.load_fasta(seqFn, rem_tVersion=False); print fasta.keys()[:10]
ShapeFn = "test_shape.out"
shape = General.load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None); print shape.keys()[:5]


#####################
#  predict_structure(sequence, shape_list=[], bp_constraint=[], mfe=True, clean=True, si=-0.6, sm=1.8, md=None, verbose=False)
#####################

seq_topredict = fasta['ENST00000429728.1']
shape_topredict = shape['ENST00000429728.1']

print Structure.predict_structure(seq_topredict)
print Structure.predict_structure(seq_topredict, shape_topredict, verbose=True)

dot_list = Structure.predict_structure(seq_topredict, shape_topredict, mfe=False)
dot_list = Structure.predict_structure(seq_topredict, shape_topredict, mfe=False, md=100)

#####################
#  bi_fold(seq_1, seq_2, local_pairing=False, mfe=True, clean=True, verbose=False)
#####################

seq_1 = fasta['ENST00000429728.1'][:100]
seq_2 = fasta['ENST00000429728.1'][-100:]

print Structure.bi_fold(seq_1, seq_2)
print Structure.bi_fold(seq_1, seq_2, local_pairing=True, verbose=True)
interactions = Structure.bi_fold(seq_1, seq_2, local_pairing=True, mfe=False)

#####################
#  search_TT_cross_linking(sequence, dot)
#####################

seq_1 = fasta['ENST00000429728.1'][:100]
seq_2 = fasta['ENST00000429728.1'][-100:]

seq_combine = seq_1 + "III" +  seq_2
dot =  Structure.bi_fold(seq_1, seq_2)

print Structure.search_TT_cross_linking(seq_combine, dot)

#####################
#  dyalign(seq_1, seq_2, shape_list_1=[], shape_list_2=[], clean=True, thread_nums=1, verbose=False)
#####################

seq_1 = "AGCTCTCTTGTGATCACCCACATTGGGTGAGGCATAAAATGAGTG"
seq_2 = "GGCTTCTCTCGTGATCACCCACATTGGGTGAGACATAAGATGAGTC"
shape_list_2 = ['0.129', '0.129', '0.0', '0.258', '0.193', '0.0', '0.064', '0.129', '0.129', '0.258', '0.064', '0.193', '0.258', '0.193', '0.0', '0.064', '0.0', '0.0', '0.0', '0.129', '0.451', '1.0', '1.0', '1.0', '0.921', '0.0', '0.0', '0.0', '0.0', '0.064', '0.0', '0.0', '0.643', '0.064', '0.064', '0.257', '0.449', '0.257', '0.257', '0.321', '0.257', '0.0', '0.064', '0.064', '0.0', '0.0']

dot_1, dot_2, aligned_seq_1, aligned_seq_2, aligned_dot_1, aligned_dot_2 = Structure.dyalign(seq_1, seq_2)
dot_1, dot_2, aligned_seq_1, aligned_seq_2, aligned_dot_1, aligned_dot_2 = Structure.dyalign(seq_1, seq_2, shape_list_2=shape_list_2, verbose=True)

import Colors
print aligned_seq_1+"\n"+aligned_dot_1+"\n"+aligned_seq_2+"\n"+aligned_dot_2+"\n"+Colors.color_SHAPE(shape_list_2)

#####################
#  multialign(seq_list, shape_list_list=[], si=-0.6, sm=1.8, thread_nums=1, clean=True, verbose=False)
#####################

seq_1 = "GGCTTCTCTCGTGATCACCCACATTGGGTGAGACATAAGATGAGTC"
seq_2 = "AGCTCTCTTGTGATCACCCACATTGGGTGAGGCATAAAATGAGTG"
seq_3 = "AATATCTTTGTGGTCACCCACATTGGGTGAGACAGAAAATGAATC"
seq_4 = "AAATATGTTTGTGGTCACCCACATTGGGTGAGACAGAAAATGAATC"
shape_1 = ['0.129', '0.129', '0.0', '0.258', '0.193', '0.0', '0.064', '0.129', '0.129', '0.258', '0.064', '0.193', '0.258', '0.193', '0.0', '0.064', '0.0', '0.0', '0.0', '0.129', '0.451', '1.0', '1.0', '1.0', '0.921', '0.0', '0.0', '0.0', '0.0', '0.064', '0.0', '0.0', '0.643', '0.064', '0.064', '0.257', '0.449', '0.257', '0.257', '0.321', '0.257', '0.0', '0.064', '0.064', '0.0', '0.0']

dot_list, aligned_seq_list, aligned_dot_list = Structure.multialign([seq_1,seq_2,seq_3,seq_4], shape_list_list=[shape_1])

for seq,dot in zip(aligned_seq_list,aligned_dot_list):
    print seq
    print dot

#####################
#  dot2ct(dot)
#####################

print Structure.dot2ct(".((((..(((((.(((((((...)))))))..)))))...)))).")

#####################
#  dot2bpmap(dot)
#####################

print Structure.dot2bpmap(".((((..(((((.(((((((...)))))))..)))))...)))).")

#####################
#  dot_from_ctFile(ctFile, number=1)
#####################

ctFn = "test_structure.ct"
dot = Structure.dot_from_ctFile(ctFn, number=3)

#####################
#  find_stem_loop(ss, max_loop_len=4, max_stem_gap=3, min_stem_len=5)
#####################

seq_topredict = fasta['ENST00000429728.1']
shape_topredict = shape['ENST00000429728.1']
dot = Structure.predict_structure(seq_topredict)

stemloops = Structure.find_stem_loop(dot, max_loop_len=4, max_stem_gap=3, min_stem_len=5)

#####################
#  find_bulge_interiorLoop(dot, stem)
#####################

print Structure.find_bulge_interiorLoop(dot, stemloops[0])
print Structure.find_bulge_interiorLoop(dot, stemloops[1]), dot[stemloops[1][0]-1:stemloops[1][3]]
print Structure.find_bulge_interiorLoop(dot, stemloops[2])

#####################
#  calcSHAPEStructureScore(dot, shape_list, stem, params={}, report=False)
#####################

sub_shape_0 = shape_topredict[ stemloops[0][0]-1:stemloops[0][3] ]
sub_dot_0 = dot[ stemloops[0][0]-1:stemloops[0][3] ]
sub_shape_1 = shape_topredict[ stemloops[1][0]-1:stemloops[1][3] ]
sub_dot_1 = dot[ stemloops[1][0]-1:stemloops[1][3] ]

print sub_dot_0+"\n"+Colors.color_SHAPE(sub_shape_0)
print sub_dot_1+"\n"+Colors.color_SHAPE(sub_shape_1)

print Structure.calcSHAPEStructureScore(dot, shape_topredict, stemloops[0])
print Structure.calcSHAPEStructureScore(dot, shape_topredict, stemloops[1])

#####################
#  sliding_score_stemloop(sequence, shape_list=None, max_loop_len=8, max_stem_gap=3, min_stem_len=5, start=0, end=None, wSize=200, wStep=100)
#####################

seq_topredict = fasta['ENST00000429728.1']
shape_topredict = shape['ENST00000429728.1']

stemloops = Structure.sliding_score_stemloop(seq_topredict, max_stem_gap=2)
stemloops = Structure.sliding_score_stemloop(seq_topredict, shape_topredict, max_stem_gap=2)

#####################
#  multi_alignment(seq_list, clean=True, verbose=False)
#####################

import General
multi_seqs = General.load_fasta("multi_alignment.fa")

aligned_seq_list = Structure.multi_alignment(multi_seqs.values(), clean=True, verbose=False)

#####################
#  kalign_alignment(seq_list, clean=True, verbose=False)
#####################

import General
multi_seqs = General.load_fasta("multi_alignment.fa")

aligned_seq_list = Structure.kalign_alignment(multi_seqs.values(), clean=True, verbose=False)

#####################
#  align_find(aligned_seq, sub_seq)
#####################

print Structure.align_find(aligned_seq_list[0], "CAGCTCCCGGGCACGCCAGCCCCCGGGTTGGGCTGGGTACGGCTGGGGCT")

#####################
#  locate_homoseq(seq_dict, main_id, align_func, main_region=None, sub_seq=None)
#####################

import General
multi_seqs = General.load_fasta("multi_alignment.fa")

homoseq = Structure.locate_homoseq(multi_seqs, main_id="XM_005250976.3", align_func=Structure.kalign_alignment, main_region=None, sub_seq="CATCGAGGAGGAGATCCGCGTGGTGCGCCTGCAGTTGGAGGCCACCGAGCGCCAGCGTGGC")

#####################
#  dot_to_alignDot(dot, aligned_seq)
#  shape_to_alignSHAPE(shape_list, aligned_seq)
#####################

dot = '......(((.(((.(((((((...)))))))..))).)))......'
aligned_seq = '-GGCTTC-TCTCGTGATCACCCACATTGGGTGAGACATAAGAT-G-AGTC-'
shape_list = ['0.129', '0.129', '0.0', '0.258', '0.193', '0.0', '0.064', '0.129', '0.129', '0.258', '0.064', '0.193', '0.258', '0.193', '0.0', '0.064', '0.0', '0.0', '0.0', '0.129', '0.451', '1.0', '1.0', '1.0', '0.921', '0.0', '0.0', '0.0', '0.0', '0.064', '0.0', '0.0', '0.643', '0.064', '0.064', '0.257', '0.449', '0.257', '0.257', '0.321', '0.257', '0.0', '0.064', '0.064', '0.0', '0.0']

print Structure.dot_to_alignDot(dot, aligned_seq)
print Structure.shape_to_alignSHAPE(shape_list, aligned_seq)

#####################
#  annotate_covariation(ref_aligned_seq, input_aligned_seq, ref_aligned_dot, anno_loop=False)
#####################

seq_1 = "GGCTTCTCTCGTGATCACCCACATTGGGTGAGACATAAGATGAGTC"
seq_2 = "AGCTCTCTTGTGATCACCCACATTGGGTGAGGCATAAAATGAGTG"
seq_3 = "AATATCTTTGTGGTCACCCACATTGGGTGAGACAGAAAATGAATC"
seq_4 = "AAATATGTTTGTGGTCACCCACATTGGGTGAGACAGAAAATGAATC"
shape_1 = ['0.129', '0.129', '0.0', '0.258', '0.193', '0.0', '0.064', '0.129', '0.129', '0.258', '0.064', '0.193', '0.258', '0.193', '0.0', '0.064', '0.0', '0.0', '0.0', '0.129', '0.451', '1.0', '1.0', '1.0', '0.921', '0.0', '0.0', '0.0', '0.0', '0.064', '0.0', '0.0', '0.643', '0.064', '0.064', '0.257', '0.449', '0.257', '0.257', '0.321', '0.257', '0.0', '0.064', '0.064', '0.0', '0.0']

dot_list, aligned_seq_list, aligned_dot_list = Structure.multialign([seq_1,seq_2,seq_3,seq_4], shape_list_list=[shape_1])

for seq,dot in zip(aligned_seq_list,aligned_dot_list):
    annot_seq = Structure.annotate_covariation(aligned_seq_list[1], seq, aligned_dot_list[1], anno_loop=True)
    print annot_seq



