#-*- coding:utf-8 -*-

import sys, os, shutil, re, logging, subprocess, string, io, argparse, bisect, concurrent
from os.path import exists, join, getsize, isfile, isdir, abspath, basename
import numpy as np
import pandas as pd
from tqdm.auto import trange, tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

def load_fasta(seqFn, rem_tVersion=False, load_annotation=False, full_line_as_id=False):
    """
    seqFn               -- Fasta file or input handle (with readline implementation)
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    load_annotation     -- Load sequence annotation
    full_line_as_id     -- Use the full head line (starts with >) as sequence ID. Can not be specified simutanouly with load_annotation

    Return:
        {tid1: seq1, ...} if load_annotation==False
        {tid1: seq1, ...},{tid1: annot1, ...} if load_annotation==True
    """
    if load_annotation and full_line_as_id:
        raise RuntimeError("Error: load_annotation and full_line_as_id can not be specified simutanouly")
    if rem_tVersion and full_line_as_id:
        raise RuntimeError("Error: rem_tVersion and full_line_as_id can not be specified simutanouly")

    fasta = {}
    annotation = {}
    cur_tid = ''
    cur_seq = ''
    
    if isinstance(seqFn, str):
        IN = open(seqFn)
    elif hasattr(seqFn, 'readline'):
        IN = seqFn
    else:
        raise RuntimeError(f"Expected seqFn: {type(seqFn)}")
    for line in IN:
        if line[0] == '>':
            if cur_seq != '':
                fasta[cur_tid] = re.sub(r"\s", "", cur_seq)
                cur_seq = ''
            data = line[1:-1].split(None, 1)
            cur_tid = line[1:-1] if full_line_as_id else data[0]
            annotation[cur_tid] = data[1] if len(data)==2 else ""
            if rem_tVersion and '.' in cur_tid: 
                cur_tid = ".".join(cur_tid.split(".")[:-1])
        elif cur_tid != '':
            cur_seq += line.rstrip()
    
    if isinstance(seqFn, str):
        IN.close()
    
    if cur_seq != '':
        fasta[cur_tid] = re.sub(r"\s", "", cur_seq)
    
    if load_annotation:
        return fasta, annotation
    else:
        return fasta

def load_stockholm(stoFn):
    """
    Read stockholm file
    
    Return:
        [ (id2seq_dict, "...(((...)))...", "AGCTGACG..AGCTG"), ... ]
    """
    from Bio import AlignIO
    from Bio.Alphabet import generic_rna
    
    alignment_list = []
    for record in AlignIO.parse(stoFn, "stockholm"):
        alignObjs = list(iter(record))
        id2seq = { alignObj.id:str(alignObj.seq) for alignObj in alignObjs }
        refStr = record.column_annotations.get("secondary_structure", "")
        refAnnot = record.column_annotations.get("reference_annotation", "")
        alignment_list.append((id2seq, refStr, refAnnot))
    
    return alignment_list

def write_fasta(fasta, seqFn, annotation={}, line_num=60):
    """
    fasta               -- A dictionary { tid1: seq1, tid2: seq2, ... }
    seqFn               -- Fasta file
    annotation          -- {id1:annot1, ...}
    """
    if line_num <= 0:
        line_num = 1_000_000_000_000
    
    OUT = open(seqFn, 'w')
    for tid in fasta:
        if tid in annotation:
            print(f">{tid} {annotation[tid]}", file=OUT)
        else:
            print(f">{tid}", file=OUT)
        full_seq = fasta[tid]
        start = 0
        while start<len(full_seq):
            print(full_seq[start:start+line_num], file=OUT)
            start+=line_num
    
    OUT.close()

def load_dot(dotFn, rem_tVersion=False, load_annotation=False):
    """
    dotFn               -- Dot file
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    load_annotation     -- Load sequence annotation
    
    Return:
        {tid1: [seq1, dot1], ...} if load_annotation==False
        {tid1: [seq1, dot1], ...},{tid1: annot1, ...} if load_annotation==True
    """
    Dot = {}
    annotation = {}
    cur_tid = ""
    
    for line in open(dotFn):
        if line[0] == '>':
            data = line[1:].split()
            cur_tid, annot_str = data[0], " ".join(data[1:])
            annotation[cur_tid] = annot_str
            if rem_tVersion and '.' in cur_tid: 
                cur_tid = ".".join(cur_tid.split(".")[:-1])
            Dot[cur_tid] = []
        else:
            content = line.strip()
            if content:
                Dot[cur_tid].append( content.split()[0] )
    
    ## check
    for tid in Dot:
        if len(Dot[tid]) != 2:
            sys.stderr.writelines("Format Error: "+cur_tid+"\n")
            raise NameError("Format Error: "+cur_tid)
    
    if load_annotation:
        return Dot, annotation
    else:
        return Dot

def write_dot(dot, dotFn):
    """
    dot                 -- A dictionary { tid1: (seq, dot), tid2: (seq, dot), ... }
    dotFn               -- A dot file
    """
    
    OUT = open(dotFn, 'w')
    for trans_id in dot:
        OUT.writelines('>%s\n%s\n%s\n' % (trans_id, dot[trans_id][0], dot[trans_id][1]))
    
    OUT.close()

def load_shape(ShapeFn, min_ratio=0, min_valid_count=0, rem_tVersion=False, min_RPKM=None):
    """
    ShapeFn             -- Standard icSHAPE file
    min_ratio           -- Mimimun ratio of valid shape
    min_valid_count     -- Mimimun count of valid shape
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    min_RPKM            -- Minimum RPKM
    """
    SHAPE = {}
    
    for line in open(ShapeFn):
        data = line.strip().split()
        transID, transLen, transRPKM = data[0], int(data[1]), data[2]
        if rem_tVersion and '.' in transID: 
            transID = ".".join(transID.split(".")[:-1])
        
        if min_RPKM and float(transRPKM) < minRPKM: continue
        
        total = len(data[3:])
        valid_num = total - data[3:].count('NULL')
        if valid_num>=min_valid_count and valid_num/total>=min_ratio:
            SHAPE[ transID ] = data[3:]
    
    return SHAPE

def write_shape(Shape, ShapeFn, RPKM={}):
    """
    Shape               -- {tid: [score1, score2, ...]}
    ShapeFn             -- Path of file to save
    RPKM                -- RPKM dictionary
    """
    OUT = open(ShapeFn, 'w')
    for tid in Shape:
        print(f"{tid}\t{len(Shape[tid])}\t{RPKM.get(tid, '*')}\t", file=OUT, end="")
        print("\t".join(['NULL' if d=='NULL' else str(d) for d in Shape[tid]]), file=OUT)
    OUT.close()

def load_SHAPEMap(shapeFn, relocate=False, loadAll=False):
    """
    Read SHAPE Map file produced by shapemapper2
    
    relocate            -- Locate the column index or use default output from shapemapper
    loadAll             -- Load all informative columns, or only load sequence and shape
    
    Return { key => list/sequence }
       keys:
        seq                 -- Sequence
        mod_list            -- A list of number of mutation in treated sample
        mod_cov_list        -- A list of coverage of mutation in treated sample
        dmso_list           -- A list of number of mutation in untreated sample
        dmso_cov_list       -- A list of coverage of mutation in untreated sample
        hq_pro_list         -- A list of raw shape reactivity
        hq_std_list         -- A list of stderr of raw shape reactivity
        shape_pro_list      -- A list of normalized shape reactivity
        shape_std_list      -- A list of stderr of normalized shape reactivity
    """
    
    def fmtN(raw):
        return "NULL" if raw=="nan" else float(raw)
    
    IN = open(shapeFn)
    header = IN.readline()
    
    if relocate:
        headers = header.strip().split()
        seq_col_idx    = headers.index("Sequence")
        tr_mod_col_idx = headers.index("Modified_mutations")
        tr_cov_col_idx = headers.index("Modified_effective_depth")
        co_mod_col_idx = headers.index("Untreated_mutations")
        co_cov_col_idx = headers.index("Untreated_effective_depth")
        hq_pro_col_idx = headers.index("HQ_profile")
        hq_std_col_idx = headers.index("HQ_stderr")
        try: no_pro_col_idx = headers.index("Norm_profile")
        except ValueError: no_pro_col_idx = -1
        try: no_std_col_idx = headers.index("Norm_stderr")
        except ValueError: no_pro_col_idx = -1
    else:
        seq_col_idx = 1
        tr_mod_col_idx = 2
        tr_cov_col_idx = 4
        co_mod_col_idx = 9
        co_cov_col_idx = 11
        hq_pro_col_idx = 25
        hq_std_col_idx = 26
        no_pro_col_idx = 27
        no_std_col_idx = 28
    
    sequence = ""
    mod_list = []
    mod_cov_list = []
    dmso_list = []
    dmso_cov_list = []
    hq_pro_list = []
    hq_std_list = []
    shape_pro_list = []
    shape_std_list = []
    
    for line in IN:
        data = line.strip().split()
        
        sequence += data[seq_col_idx]
        if no_pro_col_idx != -1:
            shape_pro_list.append( fmtN(data[no_pro_col_idx]) )
        
        if loadAll:
            mod_list.append( int(data[tr_mod_col_idx]) )
            mod_cov_list.append( int(data[tr_cov_col_idx]) )
            dmso_list.append( int(data[co_mod_col_idx]) )
            dmso_cov_list.append( int(data[co_cov_col_idx]) )
            hq_pro_list.append( fmtN(data[hq_pro_col_idx]) )
            hq_std_list.append( fmtN(data[hq_std_col_idx]) )
            if no_pro_col_idx != -1:
                shape_std_list.append( fmtN(data[no_std_col_idx]) )
    
    shapemap = { "seq":sequence, "mod_list":mod_list, "mod_cov_list":mod_cov_list, "dmso_list":dmso_list, "dmso_cov_list":dmso_cov_list, 
        "hq_pro_list":hq_pro_list, "hq_std_list":hq_std_list, "shape_pro_list":shape_pro_list, "shape_std_list":shape_std_list }
    
    return shapemap

def load_ct(ctFn, load_all=False):
    """
    Read ct file
    
    ctFn                -- ct file name
    load_all            -- load all ct from ct file, or load the first one
    
    Return:
        [seq,dotList,length] if load_all==False
        {1:[seq,dotList,length], ...} if load_all==True
    """
    Ct = {}
    
    ID = 1
    ctList = []
    seq = ""
    last_id = 0
    
    seqLen = 0
    headline = ""
    
    for line in open(ctFn):
        line = line.strip()
        if line[0]=='#':
            continue
        data = line.strip().split()
        if not data[0].isdigit():
            raise RuntimeError("cf file format Error: the first item should be a digit")
        elif seqLen==0:
            seqLen = int(data[0])
            headline = line.strip()
        elif int(data[0])!=last_id+1:
            raise RuntimeError("ct file format error...")
        else:
            left_id = int(data[0])
            right_id = int(data[4])
            seq += data[1]
            if right_id != 0 and left_id<right_id:
                ctList.append((left_id, right_id))
            last_id += 1
            if left_id == seqLen:
                #print(data, last_id+1)
                Ct[ID] = [seq, ctList, seqLen, headline]
                assert seqLen==len(seq)
                last_id = 0
                seq = ""
                ctList = []
                ID += 1
                seqLen = 0
                if not load_all:
                    return Ct[1]
    
    if seq:
        Ct[ID] = [seq, ctList, seqLen]
        if seqLen != left_id:
            raise RuntimeError("ct file format error...")
    
    return Ct

def write_ct(Fasta, Dot, ctFn):
    """
    Fasta           -- A dictionary { tid1: seq1, tid2: seq2, ... }
    Dot             -- A dictionary { tid1: dotbracket1, tid2: dotbracket2, ... }
    ctFn            -- .ct file
    
    Save dot-bracket structure to .ct file
    """
    import Structure
    OUT = open(ctFn, 'w')
    for tid in set(Fasta)&set(Dot):
        seq = Fasta[tid]
        dot = Dot[tid]
        assert len(seq) == len(dot)
        bpmap = Structure.dot2bpmap(dot)
        OUT.writelines( str(len(seq))+"\t"+tid+"\n" )
        for i in range(len(seq)):
            paired_base = bpmap.get(i+1, 0)
            OUT.writelines("%5d %c %7d %4d %4d %4d\n" % ( i+1,seq[i],i,i+2,paired_base,1 ))
    OUT.close()

def init_pd_rect(rowNum, colNum, rowNames=[], colNames=[], init_value=None):
    """
    rowNum             -- Number of rows
    colNum             -- Number of columns
    rowNames           -- Name of rows
    colNames           -- Name of columns
    init_value         -- Initialization value
    
    Initialize a pandas rect
    """
    import pandas as pd
    import numpy as np
    
    if colNames:
        assert(len(colNames)==colNum)
    else:
        colNames = np.arange(colNum)
    
    if rowNames:
        assert(len(rowNames)==rowNum)
    else:
        rowNames = np.arange(rowNum)
    
    df = pd.DataFrame(np.zeros((rowNum, colNum)), index=rowNames, columns=colNames)
    if init_value == None:
        return df
    else:
        df.iloc[:,:] = init_value
        return df

def init_list_rect(rowNum, colNum, init_value=0):
    """
    rowNum             -- Number of rows
    colNum             -- Number of columns
    init_value         -- Initialization value
    
    Initialize a pandas rect
    """
    import copy
    rect = []
    
    for i in range(rowNum):
        row = []
        for j in range(colNum):
            row.append(copy.deepcopy(init_value))
        rect.append(row)
    
    return rect

def find_all_match(pattern, string):
    """
    pattern            -- Regex expression
    string             -- String
    
    Find all position range and substring
    """
    import re
    
    matches = []
    for item in re.finditer(pattern, string):
        s, e = item.start(), item.end()
        matches.append( ((s+1, e), string[s:e]) )
    
    return matches

def bi_search(item, sorted_list, retern_index=False):
    """
    pattern             -- Regex expression
    sorted_list         -- A increasingly sorted list
    retern_index        -- If retern_index==True, the index will be returned
    
    Return:
        if retern_index == False
            True if item in sorted_list
            False if item not in sorted_list
        else
            sorted_list.index(item)
    """
    start = 0   
    end = len(sorted_list) - 1
    
    while start <= end:
        middle = (start + end) // 2
        if sorted_list[middle] < item:
            start = middle + 1
        
        elif sorted_list[middle] > item:
            end = middle - 1
        
        else:
            return middle if retern_index else True
    
    return -1 if retern_index else False

def calc_gini(list_of_values):
    """
    list_of_values          -- A list of float values
    
    Return -1 if failed
    """
    length = len(list_of_values)
    total = sum(list_of_values)
    if total == 0: 
        return -1
    
    Sorted_Array = sorted(list_of_values)
    accum, giniB = 0, 0
    
    for i in Sorted_Array:
        accum += i
        giniB += accum - i / 2.0
    
    fair_area = accum * length / 2.0
    return (fair_area - giniB) / fair_area

def calc_shape_gini(shape_list, min_num=10):
    """
    shape_list          -- A list of SHAPE scores
    min_num             -- Miminum number of scores
    
    Return -1 if failed
    """
    float_shape = [ float(shape) for shape in shape_list if shape != 'NULL' ]
    
    if len(float_shape) > min_num:
        return calc_gini(float_shape)
    
    return -1

def require_exec(exec_command, warning="", exception=True):
    """
    exec_command            -- Shell command
    warning                 -- Print warning if command not found
    exception               -- Raise an exception if command not found
                               if exception is False, no warning showed
    
    Test if command in the PATH
    
    Return full path  if command found
    """
    import distutils.spawn
    exec_path = distutils.spawn.find_executable(exec_command)
    
    if not warning:
        warning = "Error: %s not found in PATH" % (exec_command, )
    
    if not exec_path and exception:
        sys.stderr.writelines(warning+"\n")
        raise NameError(warning)
    
    return exec_path

def calc_shape_structure_positive_rate(dot, shape_list, cutoff):
    """
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    cutoff              -- A cutoff to discriminate single-stranded bases and double-stranded bases
    
    Calculate the positive rate between shape scores and secondary structure
    
    Return [true postive rate, false positive rate]
    """
    
    Pos_Num = 0
    True_Pos = 0
    False_Pos = 0
    Neg_Num = 0
    for idx, code in enumerate(list(dot)):
        if shape_list[idx] != 'NULL':
            if code != ".":
                Pos_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    True_Pos += 1
                else:
                    pass
            else:
                Neg_Num += 1
                if float(shape_list[idx]) <= cutoff:
                    False_Pos += 1
                else:
                    pass
    
    return 1.0*True_Pos/Pos_Num, 1.0*False_Pos/Neg_Num

def calc_shape_structure_ROC(dot, shape_list, start=0.0, step=0.01, stop=1.0):
    """
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    step                -- Cutoff step
    
    Calculate the ROC points structure and shape scores
    
    Return [point1, point2, point3,...]
    """
    
    assert(len(dot)==len(shape_list))
    
    ROC = []
    cutoff = start-step
    while cutoff < stop + step:
        TPR, FPR = calc_shape_structure_positive_rate(dot, shape_list, cutoff)
        ROC.append( (FPR, TPR) )
        cutoff += step
    
    return ROC

def roc_curve(dot, shape_list, min_len=15):
    """
    Calculate the FPR, TPR with sklearn
    
    Return
        fpr, tpr, thresholds
        if valid length < min_len:
            return None, None, None
    """
    import sklearn, sklearn.metrics
    import numpy as np
    
    y_true = np.array([alpha for alpha in dot])=="."
    y_score = np.array([ np.nan if d=='NULL' else float(d) for d in shape_list], dtype=float)
    
    y_true_ = y_true[~np.isnan(y_score)]
    y_score_ = y_score[~np.isnan(y_score)]
    
    if len(y_true_)>=min_len:
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(y_true_, y_score_)
    else:
        fpr, tpr, thresholds = None, None, None
    return fpr, tpr, thresholds

def calc_AUC(ROC):
    """
    ROC                 -- ROC point list
    
    Return AUC
    """
    import sklearn
    import sklearn.metrics
    
    x = [it[0] for it in ROC]
    y = [it[1] for it in ROC]
    return sklearn.metrics.auc(x, y, reorder=False)

def calc_AUC_v2(dot, shape_list):
    """
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    
    Calculate the AUC between structure and shape
    
    Return [point1, point2, point3,...]
    """
    from sklearn.metrics import roc_curve, auc
    import numpy as np
    
    assert len(dot) == len(shape_list)
    assert len(dot) > 20
    
    dot_array = np.array(list(dot))
    shape_array = np.array(shape_list, dtype=str)
    
    dot_array = dot_array[shape_array!='NULL']
    shape_array = shape_array[shape_array!='NULL']
    shape_array = shape_array.astype(float)
    
    unpaired = (dot_array=='.')
    
    FPR, TPR, _ = roc_curve(unpaired, shape_array)
    AUC = auc(FPR, TPR)
    
    return AUC

def seq_entropy(sequence, seq_type='prot', allow_X=False):
    """
    Give a sequence, calculate the entropy (0-4)
    sequence: nuc or prot sequence
    seq_type: nuc or prot
    """
    import numpy as np
    assert seq_type in ('nuc', 'prot')
    
    nuc_aatypes  = ['A','T','C','G']
    prot_aatypes = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    if allow_X:
        prot_aatypes.append('X')
        nuc_aatypes.append('N')
    
    new_seq = ""
    for r in sequence:
        if seq_type == 'nuc':
            r = r.replace('U', 'T')
            if r not in nuc_aatypes:
                if allow_X:
                    r = 'N'
                else:
                    raise RuntimeError(f"Expected nuc type: {r}")
        else:
            if r not in prot_aatypes:
                if allow_X:
                    r = 'X'
                else:
                    raise RuntimeError(f"Expected prot res type: {r}")
        new_seq += r
    
    fb = nuc_aatypes if seq_type == 'nuc' else prot_aatypes
    m = {}
    for b1 in fb:
        for b2 in fb:
            m[b1+b2] = 0
    for i in range(len(new_seq)-2):
        m[new_seq[i:i+2]] += 1
    
    m = list(m.values())
    ms = sum(m)
    prob = [ii/ms for ii in m]
    entropy = sum([ -p*np.log2(p) for p in prob if p!=0 ])
    return entropy

def Run_catchKI(command, folder_list):
    """
    Run shell command, if Ctrl+C, then remove folders in folder_list
    """
    import signal
    sig = os.system(command)
    if sig == signal.SIGINT:
        for folder in folder_list:
            shutil.rmtree(folder)
        raise RuntimeError("Ctrl + C KeyboardInterrupt.")
    return sig

class ColorFormatter(logging.Formatter):
    
    def __init__(self, fmt):
        
        from IPyRSSA import Colors
        self.formats = {
            logging.DEBUG:    Colors.f(fmt, 'lightgray'),
            logging.INFO:     Colors.f(fmt, 'lightgray'),
            logging.WARNING:  Colors.f(fmt, 'yellow'),
            logging.ERROR:    Colors.f(fmt, 'red', ft='normal'),
            logging.CRITICAL: Colors.f(fmt, 'red', ft='bold')
        }
    
    def format(self, record):
        log_fmt   = self.formats.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def get_logger(name=None, file=None, use_color=False):
    """
    Get a logger object
    
    Parameters
    ----------------
    name: The node name
    file: Output the log to this file
    use_color: Print log with color
    
    Example
    ----------------
    logger = get_logger(name='lipan', file='/tmp/1.log', use_color=True)
    
    logger.info("this is a information")
    # [2023-06-20 14:21:04,836::MYTEST::INFO] this is a information
    
    logger.debug('this is a debug text')
    # [2023-06-20 14:24:24,398::MYTEST::DEBUG] this is a debug text
    
    logger.warning("this is a warning text")
    # [2023-06-20 14:21:04,837::MYTEST::WARNING] this is a warning text
    
    logger.error("this is a error text")
    # [2023-06-20 14:21:05,080::MYTEST::ERROR] this is a error text
    
    logger.critical("this is a critical error text")
    # [2023-06-20 14:39:21,941::root::CRITICAL] this is a critical error text
    
    # 等价于
    logger.log(logging.INFO, 'this is a error text')
    logger.log(logging.DEBUG, 'this is a error text')
    logger.log(logging.ERROR, 'this is a error text')
    logger.log(logging.WARNING, 'this is a error text')
    logger.log(logging.CRITICAL, 'this is a error text')
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    if use_color:
        formatter = ColorFormatter('[%(asctime)s::%(name)s::%(levelname)s] %(message)s')
    else:
        formatter = logging.Formatter('[%(asctime)s::%(name)s::%(levelname)s] %(message)s')
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    if file is not None:
        file_handler = logging.FileHandler(file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

class persistent_locals(object):
    """
    A class to record the local variables from a function. It is a good debuger.
    
    Ref: https://stackoverflow.com/questions/9186395/
    
    Example
    ------------
    @persistent_locals
    def func():
        local1 = 1
        local2 = 2

    func()
    print(func.locals)
    """
    def __init__(self, func):
        self._locals = {}
        self.func = func

    def __call__(self, *args, **kwargs):
        def tracer(frame, event, arg):
            if event=='return':
                self._locals = frame.f_locals.copy()

        # tracer is activated on next call, return or exception
        sys.setprofile(tracer)
        try:
            # trace the function call
            res = self.func(*args, **kwargs)
        finally:
            # disable tracer and replace with old one
            sys.setprofile(None)
        return res

    def clear_locals(self):
        self._locals = {}

    @property
    def locals(self):
        return self._locals

def set_seed_all(seed:int=0, random_seed:bool=True, numpy_seed:bool=True, torch_seed:bool=True):
    """
    Set seed for random, numpy and torch
    """
    if random_seed:
        import random
        random.seed(seed)
    if numpy_seed:
        import numpy as np
        np.random.seed(seed)
    if torch_seed:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

def recursive_list(root_path, recursive=True, include_dir=False, startswith=None, endswith=None, regex=None):
    """
    Recurlive list all files in PATH
    """
    def pass_filter(file):
        if startswith is not None and not file.startswith(startswith):
            return False
        if endswith is not None and not file.endswith(endswith):
            return False
        if regex is not None and re.match(regex, file) is None:
            return False
        return True

    matches = []
    for file in os.listdir(root_path):
        full_file = join(root_path, file)
        if isfile(full_file):
            if pass_filter(file):
                matches.append(full_file)
        elif isdir(full_file) and recursive:
            if include_dir and pass_filter(file):
                matches.append(full_file)
            matches += recursive_list(full_file, recursive, startswith=startswith, endswith=endswith, regex=regex)
    return matches

def get_monomer_tmscore(pdb_pd, pdb_gt, method='TMscore', seq_depend_realn=True):
    """
    Calculate monomer TMscore
    
    Parameters
    -----------------
    pdb_pd: predicted PDB file
    pdb_gt: ground truth PDB file
    method: TMscore or USalign
    seq_depend_realn: realign two structures with sequence
    
    Return
    ----------------
    Tuple of (tmscore, seqid)
    """

    assert method in ('TMscore', 'USalign')
    
    if method == 'TMscore' and shutil.which('TMscore') is None:
        print(f"TMscore not in PATH. Please install at first")
        print(f"wget https://zhanggroup.org/TM-score/TMscore.cpp -O TMscore.cpp && g++ -O3 -o TMscore TMscore.cpp")
    
    if method == 'USalign' and shutil.which('USalign') is None:
        print(f"USalign not in PATH. Please install at first")
        print(f"wget https://zhanggroup.org/US-align/bin/module/USalign.cpp -O USalign.cpp && g++ -O3 -o USalign USalign.cpp")
    
    cmd = f"{method} {pdb_pd} {pdb_gt} -outfmt 1"
    if method == 'TMscore':
        if seq_depend_realn:
            cmd += ' -seq'
        else:
            print(f"Warning: you are using TMscore without sequence realign!")
    elif method == 'USalign':
        if seq_depend_realn:
            cmd += ' -TMscore 5'
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        print(output)
        return None, None
    
    lines = output.split('\n')
    min_i = min([ i for i in range(len(lines)) if lines[i].startswith('>') ])
    head1, seq1, head2, seq2 = lines[min_i:min_i+4]
    assert len(seq1) == len(seq2)
    
    tmscore = None
    seqid = None
    l = None
    for item in head2.split():
        if item.startswith('TM-score='):
            tmscore = float(item[9:])
        elif item.startswith('seqID='):
            seqid = float(item[6:])
    
    if tmscore is None or seqid is None:
        print(output)
    
    return (tmscore, seqid)

def get_monomer_gdtscore(pdb_pd, pdb_gt, seq_depend_realn=True):
    """
    Calculate monomer GDT-TS and GDT-HA
    
    Parameters
    -----------------
    pdb_pd: predicted PDB file
    pdb_gt: ground truth PDB file
    seq_depend_realn: realign two structures with sequence
    
    Return
    ----------------
    Tuple of (GDT-TS-score, GDT-HA-score)
    """
    
    cmd = f"TMscore {pdb_pd} {pdb_gt} -outfmt 0"
    if seq_depend_realn:
        cmd += ' -seq'
    status, output = subprocess.getstatusoutput(cmd)
    if status != 0:
        print(output)
        return None, None
    lines = output.split('\n')
    gdt_scores = []
    for line in lines:
        if line.startswith(('GDT-TS-score=', 'GDT-HA-score=')):
            gdt = float(line.split()[1])
            gdt_scores.append(gdt)
    
    return gdt_scores

def get_monomer_lddt(pdb_pd, pdb_gt):
    """
    Calculate monomer lDDT
    
    Parameters
    -----------------
    pdb_pd: predicted PDB file
    pdb_gt: ground truth PDB file
    
    Return
    ----------------
    np.ndarray: lDDT for each residues
    """
    import torch, af2_features
    
    def permute_final_dims(tensor, inds: List[int]):
        zero_index = -1 * len(inds)
        first_inds = list(range(len(tensor.shape[:zero_index])))
        return tensor.permute(first_inds + [zero_index + i for i in inds])
    
    def lddt(
        all_atom_pred_pos,
        all_atom_positions,
        all_atom_mask,
        cutoff: float = 15.0,
        eps: float = 1e-10,
        per_residue: bool = True):
        
        import torch
        
        all_atom_pred_pos = all_atom_pred_pos if torch.is_tensor(all_atom_pred_pos) else torch.as_tensor(all_atom_pred_pos)
        all_atom_positions = all_atom_positions if torch.is_tensor(all_atom_positions) else torch.as_tensor(all_atom_positions)
        all_atom_mask = all_atom_mask if torch.is_tensor(all_atom_mask) else torch.as_tensor(all_atom_mask)
        
        n = all_atom_mask.shape[-2]
        dmat_true = torch.sqrt(eps + torch.sum((all_atom_positions[..., None, :] - all_atom_positions[..., None, :, :]) ** 2, dim=-1))
        
        dmat_pred = torch.sqrt(eps + torch.sum((all_atom_pred_pos[..., None, :] - all_atom_pred_pos[..., None, :, :]) ** 2, dim=-1))
        dists_to_score = ((dmat_true < cutoff) * all_atom_mask * permute_final_dims(all_atom_mask, (1, 0)) * (1.0 - torch.eye(n, device=all_atom_mask.device)))
        
        dist_l1 = torch.abs(dmat_true - dmat_pred)
        
        score = (
            (dist_l1 < 0.5).type(dist_l1.dtype)
            + (dist_l1 < 1.0).type(dist_l1.dtype)
            + (dist_l1 < 2.0).type(dist_l1.dtype)
            + (dist_l1 < 4.0).type(dist_l1.dtype)
        )
        score = score * 0.25
        
        dims = (-1,) if per_residue else (-2, -1)
        norm = 1.0 / (eps + torch.sum(dists_to_score, dim=dims))
        score = norm * (eps + torch.sum(dists_to_score * score, dim=dims))
        
        return score
    
    def lddt_ca(
        all_atom_pred_pos,
        all_atom_positions,
        all_atom_mask,
        cutoff: float = 15.0,
        eps: float = 1e-10,
        per_residue: bool = True):
        """
        compute lddt score for Ca only.
        """
        import torch
        
        all_atom_pred_pos = all_atom_pred_pos if torch.is_tensor(all_atom_pred_pos) else torch.as_tensor(all_atom_pred_pos)
        all_atom_positions = all_atom_positions if torch.is_tensor(all_atom_positions) else torch.as_tensor(all_atom_positions)
        all_atom_mask = all_atom_mask if torch.is_tensor(all_atom_mask) else torch.as_tensor(all_atom_mask)
        
        from alphafold.common import residue_constants
        ca_pos = residue_constants.atom_order["CA"]
        all_atom_pred_pos = all_atom_pred_pos[..., ca_pos, :]  # [*, N_res, 3]
        all_atom_positions = all_atom_positions[..., ca_pos, :]  # [*, N_res, 3]
        all_atom_mask = all_atom_mask[..., ca_pos : (ca_pos + 1)]  # keep dim: [*, N_res, 1]
        
        return lddt(
            all_atom_pred_pos,
            all_atom_positions,
            all_atom_mask,
            cutoff=cutoff,
            eps=eps,
            per_residue=per_residue,
        )
    
    pd_prot = af2_features.read_pdb(pdb_pd)
    gt_prot = af2_features.read_pdb(pdb_gt, full_padding=True)
    L1 = pd_prot.atom_positions.shape[0]
    L2 = gt_prot.atom_positions.shape[0]
    s1 = pd_prot.seq(True)
    s2 = gt_prot.seq(True)
    
    assert L1 == L2, f"Expect same length, but got {L1} and {L2}"
    assert s1 == s2, f"Expect same sequence, but got {s1} and {s2}"
    
    lddt = lddt_ca(pd_prot.atom_positions, gt_prot.atom_positions, gt_prot.atom_mask).cpu().numpy()
    return lddt

def parallel_get_monomer_metrics(pd_gt_pairs, method='TMscore', seq_depend_realn=True, metric_type='TM-score', num_process=25, disable_tqdm=False):
    """
    Calculate monomer TMscore with parrallel
    
    Parameters
    -----------------
    pd_gt_pairs: [(pd_pdb_1, gt_pdb_1), (pd_pdb_2, gt_pdb_2), ...]
    method: TMscore or USalign
    seq_depend_realn: realign two structures with sequence
    metric_type: TM-score, GDT, lDDT
    num_process: number of process
    disable_tqdm: disable tqdm progress bar
    
    Return
    ----------------
    List of results from get_monomer_tmscore funcs
    """
    assert metric_type in ('TM-score', 'GDT', 'lDDT')
    results = [None]*len(pd_gt_pairs)
    with ThreadPoolExecutor(max_workers=num_process) as executor:
        tasks = {}
        for i, (pd_pdb, gt_pdb) in enumerate(pd_gt_pairs):
            assert exists(pd_pdb), pd_pdb
            assert exists(gt_pdb), gt_pdb
            if metric_type == 'TM-score':
                kwds = { 
                        'pdb_pd': pd_pdb,
                        'pdb_gt': gt_pdb,
                        'method': method, 
                        'seq_depend_realn': seq_depend_realn
                }
                handle = executor.submit(get_monomer_tmscore, **kwds)
            elif metric_type == 'GDT':
                # get_monomer_gdtscore(pdb_pd, pdb_gt, seq_depend_realn=True)
                kwds = { 
                        'pdb_pd': pd_pdb,
                        'pdb_gt': gt_pdb,
                        'seq_depend_realn': seq_depend_realn
                }
                handle = executor.submit(get_monomer_gdtscore, **kwds)
            elif metric_type == 'lDDT':
                # get_monomer_lddt(pdb_pd, pdb_gt)
                kwds = { 
                        'pdb_pd': pd_pdb,
                        'pdb_gt': gt_pdb
                }
                handle = executor.submit(get_monomer_lddt, **kwds)
            tasks[handle] = i
        for handle in tqdm(concurrent.futures.as_completed(tasks), total=len(tasks), disable=disable_tqdm, dynamic_ncols=True, leave=False):
            i = tasks[handle]
            results[i]= handle.result()
    return results

##############################################
### MSA file related
##############################################

def write_msa_txt(msa, file_or_stream, q_seq=None, sort=False, annotations=None):
    """
    Write MSA to msa .txt file
    
    Parmeters
    --------------
    msa: list of msa sequences
    file_or_stream: file or stream to save results
    q_seq: query sequence. If not given, the first sequence in msa is the query sequence
    sort: Sort input MSA according to identity
    annotations: Annotations of each sequences (excluding query seqs). len(annotations)==len(msa)
    
    Return
    --------------
    None
    """
    if q_seq is not None:
        if len(msa) > 0:
            assert len(q_seq) == len(msa[0])
    else:
        q_seq, msa = msa[0], msa[1:]
    assert '-' not in q_seq, f"Expect no gap in q_seq, but got {q_seq}"
    
    if annotations is not None:
        assert isinstance(annotations, (list, tuple))
        assert len(annotations) == len(msa)
    
    id_arr = np.array([ np.mean([ r1==r2 for r1,r2 in zip(q_seq, seq) ]) for seq in msa ], dtype=np.float64)
    if sort:
        id_order = np.argsort(id_arr)[::-1]
        id_arr   = id_arr[id_order]
        msa      = [ msa[i] for i in id_order ]
        if annotations is not None:
            annotations = [ annotations[i] for i in id_order ]
    
    OUT = file_or_stream if hasattr(file_or_stream, 'write') else open(file_or_stream, 'w')
    print(q_seq, end='\t\n', file=OUT)
    for i in range(len(msa)):
        if annotations is None or annotations[i] is None or annotations[i] == "":
            print(msa[i], round(id_arr[i],3), sep='\t', file=OUT)
        else:
            print(msa[i], round(id_arr[i],3), annotations[i], sep='\t', file=OUT)
    if not hasattr(file_or_stream, 'write'):
        OUT.close()

def load_msa_txt(file_or_stream, load_id=False, load_annot=False, sort=False):
    """
    Read msa txt file
    
    Parmeters
    --------------
    file_or_stream: file or stream to read (with read method)
    load_id: read identity and return
    
    Return
    --------------
    msa: list of msa sequences, the first sequence in msa is the query sequence
    id_arr: Identity of msa sequences
    annotations: Annotations of msa sequences
    """
    msa = []
    id_arr = []
    annotations = []
    
    if hasattr(file_or_stream, 'read'):
        lines = file_or_stream.read().strip().split('\n')
    else:
        with open(file_or_stream) as IN:
            lines = IN.read().strip().split('\n')
        # lines = open(file_or_stream).read().strip().split('\n')
    
    for idx,line in enumerate(lines):
        data = line.strip().split()
        if idx == 0:
            assert len(data) == 1, f"Expect 1 element for the 1st line, but got {data} in {file_or_stream}"
            q_seq = data[0]
        else:
            if len(data) >= 2:
                id_arr.append( float(data[1]) )
            else:
                assert len(q_seq) == len(data[0])
                id_ = round(np.mean([ r1==r2 for r1,r2 in zip(q_seq, data[0]) ]), 3)
                id_arr.append(id_)
            msa.append( data[0] )
            if len(data) >= 3:
                annot = " ".join(data[2:])
                annotations.append( annot )
            else:
                annotations.append(None)
    
    id_arr = np.array(id_arr, dtype=np.float64)
    if sort:
        id_order = np.argsort(id_arr)[::-1]
        msa      = [ msa[i] for i in id_order ]
        id_arr   = id_arr[id_order]
        annotations = [ annotations[i] for i in id_order ]
    msa = [q_seq] + msa
    
    outputs = [ msa ]
    if load_id:
        outputs.append( id_arr )
    if load_annot:
        outputs.append( annotations )
    if len(outputs) == 1:
        return outputs[0]
    return outputs

def load_a3m(file, rem_lowercase=True, load_annotation=False):
    """
    file                -- A3M file or input handle (with readline implementation)
    rem_lowercase       -- Remove lowercase alphabet
    
    Return:
        - name2seq: dict
        - name2id:  dict
        if load_annotation is True:
            - name2annotation: dict 
    """
    name2seq, name2annotation = load_fasta(file, load_annotation=True)
    table = str.maketrans('', '', string.ascii_lowercase)
    
    query_name = next(iter(name2seq))
    query_seq  = name2seq[query_name]
    assert not any([ str.islower(r) or r == '-' for r in query_seq ]), f"Expect all residues are upper-case and no gap in query_seq, but got {query_seq}"
    
    name2id = {}
    for name, seq in name2seq.items():
        align_seq = seq.translate(table)
        assert len(align_seq) == len(query_seq), f"Expect same length, but got {len(align_seq)} and {len(query_seq)}"
        name2id[name] = round(np.mean([ r1==r2 for r1,r2 in zip(align_seq, query_seq) ]), 3)
        if rem_lowercase:
            name2seq[name] = align_seq
    
    if load_annotation:
        return name2seq, name2id, name2annotation
    else:
        return name2seq, name2id

def get_msa_perRes_Neff(msa, seq_id_cutoff=0.8, device='cpu', disable_tqdm=True):
    """
    Calculate Neff of MSA
    Refer to https://academic.oup.com/bioinformatics/article/36/4/1091/5556814 for more details
    Code: https://github.com/neftlon/af2-plddt-contextualized/blob/5aa87397d4deebbc2c7b95f1881873182f9e4f6f/af22c/ref_scores.py#L102
    
    Parameters
    ---------------
    msa: List of MSA sequences
    seq_id_cutoff: pairwise sequence idensity cutoff
    device: cpu or cuda
    disable_tqdm: Disable tqdm
    
    Return
    ---------------
    Neff: np.ndarray. N_eff for each residues
    """
    import torch
    import af2_features
    gap_tok_id = af2_features.rc.restypes_with_x_and_gap.index('-')
    
    num_seq = len(msa)
    len_seq = len(msa[0])
    
    # [N_seq, L_seq]
    msa_enc = torch.from_numpy(np.array([ af2_features.seq2aatype(seq) for seq in msa ], dtype=np.int32)).to(device)
    pair_seq_id = torch.zeros([num_seq, num_seq], dtype=torch.float32, device=device)
    
    for i in trange(num_seq, disable=disable_tqdm):
        pair_seq_id[i] = torch.mean((msa_enc[i][None, :] == msa_enc).float(), axis=1)
    
    inv_n_eff_weights = 1 / (pair_seq_id > seq_id_cutoff).sum(1)
    Neff = torch.sum(inv_n_eff_weights[:, None] * ( msa_enc != gap_tok_id ).float(), axis=0)
    Neff = Neff.cpu().numpy()
    
    return Neff

def print_msa(msa, len_per_row=100, annotations=None, colors=None, OUT=sys.stdout, print_id_matrix:bool=False):
    """
    Print MSA as text
    
    Parameters
    ---------------
    msa: list of protein sequences
    len_per_row: number of residues per row
    annotations: sequence titles. If provided, the length must equal to msa
    colors: Colors of sequences. If provided, the length must equal to msa
    OUT: Output stream
    
    Return
    ---------------
    None
    """
    from IPyRSSA import Colors
    assert isinstance(msa, (tuple, list))
    if annotations is None:
        annotations = [f'seq_{i}' for i in range(len(msa))]
    else:
        assert len(msa) == len(annotations)
        assert all([ isinstance(it, str) for it in annotations ])
    
    if colors is None:
        colors = ['default'] * len(msa)
    else:
        assert len(msa) == len(colors)
    
    n_seq = len(msa)
    l_seq = len(msa[0])
    # len_aln = len(msa[0])
    max_annot_len = max([len(it) for it in annotations])
    
    info = ""
    if print_id_matrix:
        import seqdb_op
        id_matrix  = pd.DataFrame(np.zeros([n_seq, n_seq], dtype=np.float32), index=annotations, columns=annotations)
        cov_matrix = pd.DataFrame(np.zeros([n_seq, n_seq], dtype=np.float32), index=annotations, columns=annotations)
        for i in range(0, n_seq):
            for j in range(0, n_seq):
                identity, coverage = seqdb_op.get_seq_identity_and_coverage(msa[i], msa[j], remove_flanking_gap=True, remove_continous_gap=10, precision=3)
                id_matrix.iloc[i, j]  = identity
                cov_matrix.iloc[i, j] = coverage
        info += "Identity\n" + str(id_matrix) + '\n\n'
        info += "Coverage\n" + str(cov_matrix) + '\n\n'
    
    aln_seq_starts = [0] * len(msa)
    for start in range(0, l_seq, len_per_row):
        end = min(start + len_per_row, l_seq)
        for idx,seq in enumerate(msa):
            seq_start   = aln_seq_starts[idx]
            seq_tok_num = end - start - seq[start:end].count('-')
            next_seq_start = seq_start + seq_tok_num
            aln_seq_starts[idx] = next_seq_start
            seq_end = max(next_seq_start - 1, seq_start)
            d = 8
            info += annotations[idx].ljust(max_annot_len + 1) + str(seq_start).rjust(5) + ' ' + \
                    Colors.f(seq[start:end], colors[idx]) + ' ' + \
                    str(seq_end).ljust(5) + '\n'
        info += '\n'
    print(info, file=OUT, flush=True)

def read_msa_block(msa_block_file, load_id=False, load_annot=False, sort=False):
    """
    MSA block file:
        >msa1
        ABCDEFGHIJ
        ABC-EF-H-J 0.8 seq1
        -BCD-F-HIJ 0.8 seq2
        >msa2
        ABCDEFGHIJ
        ABC-EF-H-J 0.8 seq1
        -BCD-F-HIJ 0.8 seq2
    
    Return
    ------------
    it: a iterative object
        next(it) return same object as load_msa_txt
    """
    with open(msa_block_file) as IN:
        stop_read = False
        while not stop_read:
            header = IN.readline()
            assert header[0] == '>', f"Expect header startswith >, but got {header}"
            msa_txt = ""
            while True:
                loc = IN.tell()
                line = IN.readline()
                if line == "":
                    stop_read = True
                    break
                elif line[0] == '>':
                    IN.seek(loc)
                    break
                else:
                    msa_txt += line
            yield load_msa_txt(io.StringIO(msa_txt), load_id=load_id, load_annot=load_annot, sort=sort)

##############################################
### Statistics
##############################################

def get_statistics(values, be_sorted=False, beautiful_print=True, precision=3):
    """
    Get statistics values of values
    
    Return
    -----------
    qs: list
    statistics: list
    """
    import copy
    if len(values) == 0:
        return [], []
    values = copy.deepcopy(values)
    N = len(values)
    if isinstance(values, (list, tuple)):
        values = np.array(values)
    
    if not be_sorted:
        values.sort()
    
    #qs = [0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
    met = ['Min', 'Q1', 'Q5', 'Q10', 'Q25', 'Median', 'Q75', 'Q90', 'Q95', 'Q99', 'Max', 'Mean', '#']
    statistics = []
    for m in met:
        if m == 'Min':
            statistics.append( round(values[0], precision) )
        elif m == 'Max':
            statistics.append( round(values[-1], precision) )
        elif m == 'Median':
            statistics.append( round(values[ N//2 ], precision) )
        elif m == 'Mean':
            statistics.append( round(values.mean(), precision) )
        elif m == '#':
            statistics.append( len(values) )
        elif m.startswith('Q'):
            q = int(m[1:]) / 100
            i = min(int( len(values)*q ), len(values)-1)
            statistics.append( round(values[i], precision) )
    
    if beautiful_print:
        statistics_str = [ str(d) for d in statistics ]
        max_len = max([ len(d) for d in statistics_str ]) + 2
        max_len = max(max_len, 8)
        header = ""
        content = ""
        for m,v in zip(met, statistics_str):
            header  += m.center(max_len, " ")
            content += v.center(max_len, " ")
        print(header)
        print(content)
    
    return met, statistics

##############################################
### Tensor information
##############################################

def print_tensor_info(tensor, print_source:bool = False, print_content:bool = False):
    """
    Print tensor information
    """
    import torch
    import numpy as np
    from inspect import getframeinfo, stack
    
    caller = getframeinfo(stack()[1][0])
    filename = caller.filename
    if print_source:
        print(f'File "{os.path.abspath(caller.filename)}", line {caller.lineno}')
    
    if isinstance(tensor, np.ndarray):
        print(f"{filename} ({caller.lineno}): type(tensor)={type(tensor)}, tensor.dtype={tensor.dtype}, tensor.shape={tensor.shape}")
    elif isinstance(tensor, torch.Tensor):
        print(f"{filename} ({caller.lineno}): type(tensor)={type(tensor)}, tensor.dtype={tensor.dtype}, tensor.shape={tensor.shape}, tensor.device={tensor.device}")
    else:
        print(f"type(tensor) = {type(tensor)}")
    
    if print_content:
        print(tensor)

def print_dict(dict_var: dict, color:bool = True, indent:int = 0, max_elem:int = None, max_layer:int = None, layer:int = 0):
    """
    Print dictionary content
    
    Parameters
    --------------
    dict_var: dict
    color: bool. Print with color
    indent: indent
    max_elem: max element to print
    max_layer: max recursive layers
    layer: current layer
    """
    import torch
    import numpy as np
    def mtype(var):
        return f"{type(var).__module__}.{type(var).__qualname__}"
    if color:
        from IPyRSSA.Colors import f
        fk = lambda x: f(x, 'blue')
    else:
        fk = lambda x: x
    if not isinstance(dict_var, dict):
        raise RuntimeError(f"Expect input dict but got {type(dict_var)}")
    if max_elem is None:
        max_elem = 1_000_0000_000
    max_elem = max_elem + 1 if max_elem % 2 == 1 else max_elem
    assert max_elem >= 1
    for i,(k,v) in enumerate(dict_var.items()):
        indent_prefix = " " * indent
        prefix = indent_prefix + fk(f"{k}")
        if i < max_elem // 2 or i >= len(dict_var) - max_elem // 2:
            if isinstance(v, (list, tuple)):
                print(f"{prefix}: type={mtype(v)}, len(v)={len(v)}")
            elif isinstance(v, dict):
                print(f"{prefix}: type={mtype(v)}, len(v)={len(v)}")
                if max_layer is None or max_layer > layer:
                    print_dict(v, color=color, indent=indent+4, max_elem=max_elem, max_layer=max_layer, layer=layer+1)
                elif max_layer is not None:
                    print(f"{indent_prefix}    ....")
            elif isinstance(v, (int, float)):
                print(f"{prefix}: type={mtype(v)}, v={v}")
            elif isinstance(v, (str, bytes)):
                print(f"{prefix}: type={mtype(v)}, v[:20]={v[:20]}")
            elif isinstance(v, np.ndarray):
                print(f"{prefix}: type={mtype(v)}, v.dtype={v.dtype}, v.shape={v.shape}")
            elif isinstance(v, torch.Tensor):
                print(f"{prefix}: type={mtype(v)}, v.dtype={v.dtype}, v.shape={v.shape}, v.device={v.device}")
            elif isinstance(v, object):
                print(f"{prefix}: type={mtype(v)}")
            else:
                print(f"{prefix}: type={mtype(v)}")
        elif i == max_elem // 2:
            print(f"{indent_prefix}....")
            # break

##############################################
### TXT index
##############################################

def build_txt_index(big_file, index_file, sep='\t', skip_head=0, n_parts=1000, disable_tqdm=False):
    """
    Build index file for text-format file. The first column must be numbers and sorted
    """
    linecount = int(subprocess.getoutput(f'wc -l {big_file}').split()[0])
    number_list  = []
    number_bytes = []
    start_bytes = 0
    last_index  = -999999
    with open(big_file) as IN:
        for i,line in enumerate(tqdm(IN, total=linecount, dynamic_ncols=True, disable=disable_tqdm)):
            if i < skip_head:
                start_bytes += len(line)
                continue
            index = int(line[:line.index(sep)])
            if last_index > index:
                raise RuntimeError(f"Error: last_index={last_index}, index={index} at line {i}")
            number_list.append(index)
            number_bytes.append(start_bytes)
            start_bytes += len(line)
            last_index = index
    assert len(number_list) > n_parts
    index_list = []
    part_size = len(number_list) // n_parts
    for i in range(0, n_parts):
        j = part_size * i
        while j > 0 and number_list[j-1] == number_list[j]:
            j -= 1
        index_list.append(( number_list[j], number_bytes[j] ))
    with open(index_file, 'w') as OUT:
        for index,start_bytes in index_list:
            print(index, start_bytes, sep='\t', file=OUT)

def read_txt_with_index(big_file, index, items, sep='\t', max_lines=1000000, disable_tqdm=False):
    """
    Read big file with index
    """
    data = []
    with open(big_file) as IN:
        for item in tqdm(items, disable=disable_tqdm):
            i = bisect.bisect_left(index.index_array, item)
            if i >= len(index.index_array) or (index.index_array[i] > item and i > 0):
                i -= 1
            _ = IN.seek(index.start_bytes_array[i])
            for line in IN:
                cur_item = int(line[:line.index(sep)])
                if cur_item > item:
                    break
                elif cur_item == item:
                    data.append(line)
                    if len(data) >= max_lines:
                        break
    return data

def read_index(index_file):
    """
    Read index file build from build_txt_index
    """
    index_array, start_bytes_array = [], []
    with open(index_file) as IN:
        for line in IN:
            index, byte = [ int(x) for x in line.strip().split() ]
            index_array.append(index)
            start_bytes_array.append(byte)
    index_array = np.array(index_array)
    start_bytes_array = np.array(start_bytes_array)
    index = argparse.Namespace(**{ 'index_array': index_array, 'start_bytes_array': start_bytes_array })
    return index

##############################################
### DDP
##############################################

def print_deepspeed_ddp_env(only_once=False):
    attrs = ['MASTER_ADDR', 'MASTER_PORT', 'CROSS_RANK', 'CROSS_SIZE', 'LOCAL_SIZE', 'LOCAL_RANK', 'WORLD_SIZE', 'RANK']
    rank = int(os.environ['RANK'])
    if not only_once or rank == 0:
        # '\x1b[0;33;49m222\x1b[0m'
        print("\x1b[0;33;49mDeepSpeed DDP Environment Variables\x1b[0m")
        for attr in attrs:
            print(f"{attr}: {os.environ.get(attr, None)}")

def print_torch_ddp_env(only_once=False):
    attrs = ['MASTER_ADDR', 'MASTER_PORT', 'WORLD_SIZE', 'RANK', 'LOCAL_WORLD_SIZE', 'LOCAL_RANK', 'GROUP_WORLD_SIZE', 'GROUP_RANK', 'ROLE_WORLD_SIZE', 'ROLE_RANK']
    rank = int(os.environ['RANK'])
    if not only_once or rank == 0:
        print("\x1b[0;33;49mTorch DDP Environment Variables\x1b[0m")
        for attr in attrs:
            print(f"{attr}: {os.environ.get(attr, None)}")


