#-*- coding:utf-8 -*-

import sys, commands

def load_fasta(seqFn, rem_tVersion=False):
    """
    seqFn               -- Fasta file
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    """
    Fasta = {}
    cur_tid = ''
    
    for line in open(seqFn):
        if line[0] == '>':
            cur_tid = line[1:].split()[0]
            if rem_tVersion and '.' in cur_tid: 
                cur_tid = ".".join(cur_tid.split(".")[:-1])
            Fasta[ cur_tid ] = ''
        else:
            Fasta[ cur_tid ] += line.strip()
    
    return Fasta

def write_fasta(Fasta, seqFn):
    """
    Fasta               -- A dictionary { tid1: seq1, tid2: seq2, ... }
    seqFn               -- Fasta file
    """
    from Seq import flat_seq
    
    OUT = open(seqFn, 'w')
    for trans_id in Fasta:
        print >>OUT, '>%s\n%s' % (trans_id, flat_seq(Fasta[trans_id]))
    
    OUT.close()

def load_dot(dotFn, rem_tVersion=False):
    """
    dotFn               -- Dot file
    rem_tVersion        -- Remove version information. ENST000000022311.2 => ENST000000022311
    """
    Dot = {}
    cur_tid = ""
    
    for line in open(dotFn):
        if line[0] == '>':
            cur_tid = line[1:].split()[0]
            if rem_tVersion and '.' in cur_tid: 
                cur_tid = ".".join(cur_tid.split(".")[:-1])
            Dot[cur_tid] = []
        else:
            content = line.strip()
            if content:
                Dot[cur_tid].append( content )
    
    ## check
    for tid in Dot:
        if len(Dot[tid]) != 2:
            print >>sys.stderr, "Format Error: "+cur_tid
            raise NameError("Format Error: "+cur_tid)
    
    return Dot

def load_shape(ShapeFn, rem_tVersion=False, min_RPKM=None):
    """
    ShapeFn             -- Standard icSHAPE file
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
        
        SHAPE[ transID ] = data[3:]
    
    return SHAPE

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
    rect = []
    
    for i in range(rowNum):
        row = []
        for j in range(colNum):
            row.append(init_value)
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

def bi_search(item, sorted_list):
    """
    pattern            -- Regex expression
    sorted_list        -- A increasingly sorted list
    """
    start = 0   
    end = len(sorted_list) - 1
    
    while start <= end:
        middle = (start + end) / 2
        if sorted_list[middle] < item:
            start = middle + 1
        
        elif sorted_list[middle] > item:
            end = middle - 1
        
        else:   
            return True
    
    return False

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
        print >>sys.stderr, warning
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

def calc_shape_structure_ROC(dot, shape_list, step=0.01):
    """
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    step                -- Cutoff step
    
    Calculate the ROC points structure and shape scores
    
    Return [point1, point2, point3,...]
    """
    
    assert(len(dot)==len(shape_list))
    
    ROC = []
    cutoff = -step
    while cutoff < 1.0 + step:
        TPR, FPR = calc_shape_structure_positive_rate(dot, shape_list, cutoff)
        ROC.append( (FPR, TPR) )
        cutoff += step
    
    return ROC

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


