#########
#########   Test Seq.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import Colors

#####################
#  f(text, fc='red', bc='default', ft='normal')
#####################

print Colors.f("Hello ", fc='red') + Colors.f("World", fc='cyan')
print Colors.f("Hello ", fc='white', bc='red') + Colors.f("World", fc='cyan')
print Colors.f("Hello ", fc='white', bc='red', ft='blink') + Colors.f("World", fc='cyan')

#####################
#  color_SHAPE(shape_list, cutoff=[0.3, 0.5, 0.7])
#  color_Seq_SHAPE(sequence, shape_list, cutoff=[0.3, 0.5, 0.7])
#####################

import General
seqFn = "test_seq.fasta"
fasta = General.load_fasta(seqFn, rem_tVersion=False)
ShapeFn = "test_shape.out"
shape = General.load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None)

print Colors.color_SHAPE(shape['ENST00000429728.1'][40:190])
print Colors.color_Seq_SHAPE(fasta['ENST00000429728.1'][40:190], shape['ENST00000429728.1'][40:190])

