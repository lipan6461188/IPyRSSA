#########
#########   Test Visual.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import Visual

import General
ShapeFn = "test_shape.out"
shape = General.load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None)
dotFn = "test_structure.dot"
dot = General.load_dot(dotFn, rem_tVersion=False)[]
seqFn = "test_seq.fasta"
fasta = General.load_fasta(seqFn, rem_tVersion=False)

sequence = fasta['ENST00000558492.1']
dotBracket = dot['ENST00000558492.1'][1]
shape_list = shape['ENST00000558492.1']

#####################
#  Plot_RNAStructure_Shape(sequence, dot, shape_list, mode='label', correctT=True, highlight_region=[], title="", wait=True, VARNAProg=VARNAProg)
#####################

print Visual.Plot_RNAStructure_Shape(sequence, dotBracket, shape_list, mode='fill')
print Visual.Plot_RNAStructure_Shape(sequence, dotBracket, shape_list, mode='label')
print Visual.Plot_RNAStructure_Shape(sequence, dotBracket, shape_list, mode='heatmap')

#####################
#  Plot_RNAStructure_Base(sequence, dot, mode='fill', correctT=True, highlight_region=[], title="", wait=True, VARNAProg=VARNAProg)
#####################

print Visual.Plot_RNAStructure_Base(sequence, dotBracket, mode='fill')
print Visual.Plot_RNAStructure_Base(sequence, dotBracket, mode='label')

#####################
# Plot_RNAStructure_highlight(sequence, dot, hg_base_list=[], mode='fill', correctT=True, highlight_region=[], title="", wait=True, VARNAProg=VARNAProg):
#####################

print Visual.Plot_RNAStructure_highlight(sequence, dotBracket, hg_base_list=range(20,90), mode='fill')
print Visual.Plot_RNAStructure_highlight(sequence, dotBracket, hg_base_list=range(20,50)+range(90,120), mode='fill')


