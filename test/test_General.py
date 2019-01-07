#########
#########   Test General.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import General

#####################
#  load_fasta(seqFn, rem_tVersion=False)
#####################

seqFn = "test_seq.fasta"
fasta = General.load_fasta(seqFn, rem_tVersion=False); print fasta.keys()[:10]
fasta = General.load_fasta(seqFn, rem_tVersion=True); print fasta.keys()[:10]

#####################
#  write_fasta(Fasta, seqFn)
#####################

seqFn = "test_seq_out.fasta"
General.write_fasta(fasta, seqFn)

#####################
#  load_dot(dotFn, rem_tVersion=False)
#####################

dotFn = "test_structure.dot"
dot = General.load_dot(dotFn, rem_tVersion=False); print dot.keys()
dot = General.load_dot(dotFn, rem_tVersion=True); print dot.keys()

#####################
#  load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None)
#####################

ShapeFn = "test_shape.out"
shape = General.load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None); print shape.keys()[:5]
shape = General.load_shape(ShapeFn, rem_tVersion=True, min_RPKM=None); print shape.keys()[:5]

#####################
#  init_pd_rect(rowNum, colNum, rowNames=[], colNames=[], init_value=None)
#  init_list_rect(rowNum, colNum, init_value=0)
#####################

frame = General.init_pd_rect(rowNum=10, colNum=10, rowNames=[], colNames=[], init_value=None)
print frame

frame = General.init_list_rect(rowNum=10, colNum=10, init_value=0)
print frame

#####################
#  find_all_match(pattern, string)
#####################

General.find_all_match("ENSG", "Gene ENSG is from ENSG002822")

#####################
#  bi_search(item, sorted_list)
#####################

sorted_list = sorted([random.randint(0,20) for i in range(20)])
print sorted_list
General.bi_search(8, sorted_list)

#####################
#  calc_gini(list_of_values)
#  calc_shape_gini(shape_list, min_num=10)
#####################

ShapeFn = "test_shape.out"
shape = General.load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None); print shape.keys()[:5]
for tid in shape:
    gini = General.calc_shape_gini(shape[tid], min_num=20)
    if gini != -1:
        print tid, gini

#####################
#  require_exec(exec_command, warning="", exception=True)
#####################

General.require_exec("muscle", warning="", exception=True)
General.require_exec("muscle2", warning="Failed", exception=False)


