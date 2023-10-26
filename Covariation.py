#-*- coding:utf-8 -*-
"""

This module is designed for RNA base pairing covariation calling

########### Example

import Covariation
query_seq = "GATTTAAGTGAATAGCTTGGCTATCTCACTTCCCCTCGTTCTCTTGCAGAACTTTGATTTTAACGAACTTAAATAAAAGCCCTGTTGTTTAGCGTATCGTTGCACTTGTCTGGTGGGATTGTGGCATTAATTTGCCTGCTCATCTAGGCAGTGGACATATGCTCAACACTGGGTATAATTCTAATTGAATACTATTTTTCAGTTAGAGCGTCGTGTCTCTTGTACGTCTCGGTCACAATACACGGTTTCGTCCGGTGCGTGGCAATTCGGGGCACATCATGTCTTTCGTGGCTGGTGTGACCGCGCAAGGTGCGCGCGGTACGTATCGAGCAGCGCTCAACTCTGAAAAACATCAAGACCATGTGTCTCTAACTGTGCCACTCTGTGGTTCAGGAAACCTGGTTGAAAAACTTTCACCATGGTTCATGGATGGCGAAAATGCCTATGAAGTGGTGAAGGC"
query_dot = '.....(((((((((((...))))).))))))......(((((.....))))).(((.......)))............((((.(.((((.(((((((.(((.((((.((((((((((.....((((......))))..)))))).)))))))))))))))))).))))).)))).....................(((((((...((((((..((.(((..(((((((((((((((..(((.(((......)))))).)))))....)))).(((((((.(((......))))))))))(((((((.......))))))))))))).)))))))))))....))))))).......(((((.((.((.....(((.(((((...))))).))))).)).)))))......((((((((((..((((((...((((....)))))))))))))))))))).'
seqdbFn = "examples/Rfam_sequence.fa"
workdir_root = '/Share2/home/zhangqf7/tmp/rfam'
covary_bps = Covariation.call_covariation(query_seq, query_dot, "MERS_5UTR", seqdbFn, workdir=None,
                    nohmm=True, cmsearchE=1, cpu=20, use_LSF=True, 
                    LSF_parameters={}, progress=True, clean=False)

"""

from . import Structure, General, Colors

import os, sys

def dot2sto(dot, modelname, outfile, refSeq=None, GS_DE=None, mode='w'):
    """
    dot             -- { seqname:[aligned_seq, aligned_dot], ... }
    modelname       -- CM name
    outfile         -- Write the model to file
    refSeq          -- Reference sequence, for #=GC RF line
    GS_DE           -- The annotation for sequence, input a dictionary
                       =GS NC000283.2 DE Bat coronavirus .....
    mode            -- Cover or append
    """
    OUT = open(outfile, mode)
    print("# STOCKHOLM 1.0\n", file=OUT)
    print(f"#=GF ID   {modelname}\n", file=OUT)
    
    if GS_DE:
        overlapped_keys = list(GS_DE.keys() & dot.keys())
        longest_key = sorted(overlapped_keys,key=lambda x: len(x), reverse=True)[0]
        key_len = len(longest_key)+2
        for key in overlapped_keys:
            print(f'#=GS %-{key_len}sDE %s' % (key, GS_DE[key]), file=OUT)
        print('', file=OUT)
    
    common_dot = ""
    maxLen = max(max([len(key) for key in dot]), 15)
    for seqName in dot:
        align_seq, align_dot = dot[seqName]
        assert len(align_seq)==len(align_dot), "Sequence length should have same length with align_dot length"
        if not common_dot:
            common_dot = align_dot
        else:
            assert align_dot==common_dot, "The aligned dot should be same"
        print( seqName+" "*(maxLen+5-len(seqName)), align_seq, sep="", file=OUT )
    print("#=GC SS_cons"+" "*(maxLen+5-12), common_dot, sep="", file=OUT )
    if refSeq:
        print("#=GC RF"+" "*(maxLen+5-7), refSeq, sep="", file=OUT )
    print("\n//", file=OUT)
    OUT.close()

def cmbuild(inStoFn, outCMFn, verbose=False, showCMD=True):
    """
    Build a cm file from stockholm file
    inStoFn                     -- Input stockholm file
    outCMFn                     -- Output CM file
    verbose                     -- Show command and log information
    showCMD                     -- Print the command

    Require: cmbuild
    """
    import General
    import shutil
    
    cmbuild_exe = General.require_exec("cmbuild", exception=True)
    cmd = f"{cmbuild_exe} -F "
    if not verbose:
        cmd += '-o /dev/null '
    cmd += f"{outCMFn} {inStoFn}"
    
    if showCMD:
        import Colors
        print( Colors.f(cmd, fc='yellow') )
    os.system(cmd)

def cmcalibrate(CMFn, cpu=0, verbose=True, showCMD=True, use_LSF=False, LSF_parameters={}):
    """
    Calibrate the CM model
    CMFn                -- CM file
    cpu                 -- How many CPU to use
    verbose             -- Show command and log information
    showCMD             -- Print the command
    use_LSF             -- Submit to LSF if True
    LSF_parameters      -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmcalibrate', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    
    Return:
        Return job object if use_LSF==True
        Return None if use_LSF==False

    Require: cmcalibrate
    """
    import General
    import shutil
    
    cmcalibrate_exe = General.require_exec("cmcalibrate", exception=True)
    cmd = f"{cmcalibrate_exe} "
    if cpu>0:
        cmd += f"--cpu {cpu} "
    
    cmd += CMFn
    if not verbose:
        cmd += " > /dev/null"
    if showCMD:
        import Colors
        print( Colors.f(cmd, fc='yellow') )
    
    if use_LSF:
        import Cluster
        job = Cluster.new_job(command=cmd, 
            queue=LSF_parameters.get('queue', 'Z-ZQF'), 
            cpu=LSF_parameters.get('cpu', 20), 
            job_name=LSF_parameters.get('job_name', 'cmcalibrate'), 
            logFn=LSF_parameters.get('logFn', '/dev/null'),
            errFn=LSF_parameters.get('errFn', '/dev/null'))
        job.get_submit_command()
        job.submit()
        return job
    else:
        os.system(cmd)

def cmsearch(CMFile, seqdbFn, outTXT, outSto, 
    cpu=0, toponly=False, nohmm=False, 
    nohmmonly=False, outputE=20, acceptE=1, 
    cut_ga=False, rfam=False, glocal=False,
    verbose=True, showCMD=True, use_LSF=False, LSF_parameters={}):
    """
    Search CM model from sequence database
    CMFile              -- CM file
    seqdbFn             -- File name of sequence database
    outTXT              -- Output txt file
    outSto              -- Output Stockholm file
    cpu                 -- How many threads to use
    toponly             -- Only search the top(forward) strand
    nohmm               -- skip all HMM filter stages, use only CM (slow)
    nohmmonly           -- never run HMM-only mode, not even for models with 0 basepairs
    outputE             -- report sequences <= this E-value threshold in output  [10.0]  (x>0)
    acceptE             -- consider sequences <= this E-value threshold as significant  [0.01]
    cut_ga              -- use CM's GA gathering cutoffs as reporting thresholds
    rfam                -- set heuristic filters at Rfam-level (fast)
    glocal              -- configure CM for glocal alignment [default: local]
    verbose             -- Show command and log information
    showCMD             -- Print the command
    use_LSF             -- Submit to LSF if True
    LSF_parameters      -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmsearch', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    
    Require: cmsearch
    """
    import General
    import shutil
    
    cmsearch_exe = General.require_exec("cmsearch", exception=True)
    cmd = f"{cmsearch_exe} --notextw "
    if cpu>0:
        cmd += f"--cpu {cpu} "
    if toponly:
        cmd += "--toponly "
    if nohmm:
        cmd += "--nohmm "
    if nohmmonly:
        cmd += "--nohmmonly "
    if cut_ga:
        cmd += "--cut_ga "
    if rfam:
        cmd += "--rfam "
    if glocal:
        cmd += "-g "
    cmd += f"-E {outputE} --incE {acceptE} -o {outTXT} -A {outSto} {CMFile} {seqdbFn}"
    
    if not verbose:
        cmd += " > /dev/null"
    if showCMD:
        import Colors
        print( Colors.f(cmd, fc='yellow') )
    
    if use_LSF:
        import Cluster
        job = Cluster.new_job(command=cmd, 
            queue=LSF_parameters.get('queue', 'Z-ZQF'), 
            cpu=LSF_parameters.get('cpu', 20), 
            job_name=LSF_parameters.get('job_name', 'cmsearch'), 
            logFn=LSF_parameters.get('logFn', '/dev/null'),
            errFn=LSF_parameters.get('errFn', '/dev/null'))
        #print(job.get_submit_command())
        job.submit()
        return job
    else:
        os.system(cmd)

def cmalign(CMFile, seqdbFn, outFn, 
    cpu=0, glocal=False, outformat='Stockholm', mxsize=1028.0, 
    verbose=True, showCMD=True, use_LSF=False, LSF_parameters={}):
    """
    Search CM model from sequence database
    CMFile              -- CM file
    seqdbFn             -- File name of sequence database
    outFn               -- Output file
    cpu                 -- How many threads to use
    glocal              -- Configure CM for global alignment [default: local]
    outformat           -- Output alignment in format <s>  [Stockholm]
                           Stockholm, Pfam, AFA (aligned FASTA), A2M, Clustal, PHYLIP
    mxsize              -- Set maximum allowable DP matrix size to <x> Mb  [1028.0]  (x>0.)
    verbose             -- Show command and log information
    showCMD             -- Print the command
    use_LSF             -- Submit to LSF if True
    LSF_parameters      -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmalign', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    
    Require: cmalign
    """
    import General
    import shutil
    
    cmalign_exe = General.require_exec("cmalign", exception=True)
    cmd = f"{cmalign_exe} "
    if cpu>0:
        cmd += f"--cpu {cpu} "
    if glocal:
        cmd += "-g "
    cmd += f"--outformat {outformat} --mxsize {mxsize} -o {outFn} {CMFile} {seqdbFn}"
    
    if not verbose:
        cmd += " > /dev/null"
    if showCMD:
        import Colors
        print( Colors.f(cmd, fc='yellow') )
    
    if use_LSF:
        import Cluster
        job = Cluster.new_job(command=cmd, 
            queue=LSF_parameters.get('queue', 'Z-ZQF'), 
            cpu=LSF_parameters.get('cpu', 20), 
            job_name=LSF_parameters.get('job_name', 'cmsearch'), 
            logFn=LSF_parameters.get('logFn', '/dev/null'),
            errFn=LSF_parameters.get('errFn', '/dev/null'))
        #print(job.get_submit_command())
        job.submit()
        return job
    else:
        os.system(cmd)

def R_scape(StoFn, outDir, outname=None, maxIdentity=0.985, minIndentity=0.500, 
    F=0.5, gapthresh=0.5, two_set_test=True, fold=False, acceptE=0.05, nseqmin=5,
    verbose=False, showCMD=True):
    """
    StoFn                           -- Stockholm file
    outDir                          -- Output the file to this directory
    outname                         -- File prefix
    maxIdentity                     -- require seqs to have < <x> id  [1.0]  (0<x<=1.0)
    minIndentity                    -- require seqs to have >= <x> id  (0<=x<1.0)
    F                               -- filter out seqs <x*seq_cons residues  (0<x<=1.0)
    gapthresh                       -- keep columns with < <x> fraction of gaps  [0.75]  (0<=x<=1)
    two_set_test                    -- two-set test: basepairs / all other pairs. Requires a given structure
    fold                            -- obtain the structure with maximum covariation
    acceptE                         -- Eval: max expected number of covNBPs allowed  [0.05]  (x>=0)
    nseqmin 						-- minimum number of sequences in the alignment  (n>0)
    verbose                         -- Show command and log information
    showCMD                         -- Print the command
    
    Require: R-scape
    """
    import General
    import shutil
    
    R_scape_exe = General.require_exec("R-scape", exception=True)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    cmd = f"{R_scape_exe} --outmsa --r2rall --outtree --roc --voutput --outnull --consensus "
    cmd += f"--outdir {outDir} -I {maxIdentity} -i {minIndentity} -F {F} --gapthresh {gapthresh} -E {acceptE} --nseqmin {nseqmin} "
    if outname:
        cmd += f"--outname {outname} "
    if two_set_test:
        cmd += "-s "
    if fold:
        cmd += "--fold "
    cmd += f" {StoFn} "
    
    if not verbose:
        cmd += "> /dev/null"
    if showCMD:
        import Colors
        print(Colors.f(cmd, fc='yellow'))
    
    os.system(cmd)

def read_RScape_result(Rscape_cov_fn, RScape_inpit_sto=None, querySeq=None):
    """
    Rscape_cov_fn               -- The XXXX.cov file of R-Scape output
    RScape_inpit_sto            -- The input stockholm file for R-scape
    querySeq                    -- The qeury sequence
    
    if RScape_inpit_sto and querySeq is provided, the coordinate system will be converted into querySeq
    Return:
        [ [left,right,evalue],... ]
    """
    import numpy as np
    cov_pairs = []
    for line in open(Rscape_cov_fn):
        if line[0]=='*':
            data = line.strip().split()
            left = int(data[1])
            right = int(data[2])
            evalue = float(data[4])
            cov_pairs.append((left, right, evalue))
    cov_pairs = sorted(cov_pairs, key=lambda x: x[0])
    
    if RScape_inpit_sto and querySeq:
        import General,Structure
        id2seq, refStr, refSeq = General.load_stockholm(RScape_inpit_sto)[0]
        for id1 in id2seq:
            if 'U' in id2seq[id1]:
                querySeq = querySeq.replace('T','U')
            elif 'T' in id2seq[id1]:
                querySeq = querySeq.replace('U','T')
            aligned_seq1, aligned_seq2 = Structure.multi_alignment( [id2seq[id1].replace('-',''), querySeq] )
            if np.mean([ d1==d2 for d1,d2 in zip(aligned_seq1, aligned_seq2) ]) > 0.9:
                break
        pos2pos = {}
        j = 0
        for i in range(len(id2seq[id1])):
            if id2seq[id1][i] != '-':
                while j<len(aligned_seq1) and aligned_seq1[j] == '-':
                    j += 1
                if j>=len(aligned_seq1):
                    break
                pos2pos[i] = j
                j += 1
        cov_pairs = [ [pos2pos.get(d[0]-1,-1)+1, pos2pos.get(d[1]-1,-1)+1, d[2]] for d in cov_pairs ]
    
    return cov_pairs

def collapse_sequences(id2seq, refSeq, max_identity=0.95, min_match_identity=0.5, max_indel_ratio=0.5):
    """
    Remove sequences with identity larger than max_identity and less than min_identity
    id2seq                      -- {id1:aligned_seq1, id2:aligned_seq2, ...}
    refSeq                      -- Reference aligned sequence
    max_identity                -- Maximum identity (0-1)
    min_match_identity          -- Minimum identity for match region (match or mismatch) (0-1)
    max_indel_ratio             -- Maximum indel ratio (0-1)
                                   A-[-] or [-]-A
    
    Return: 
        {id1:aligned_seq1, id2:aligned_seq2, ...}
    """
    def calc_seq_identity(seq1, seq2):
        assert len(seq1)==len(seq2), f"Error: len(seq1)={len(seq1)}; len(seq2)={len(seq2)}"
        total = 0 # Either of side have nucleotides
        match_num = 0 # Both of side have nucleotides
        exact_match_num = 0 # Both of side have nucleotides and same nucleotide
        for b1,b2 in zip(seq1, seq2):
            if b1!='-' and b2!='-':
                if b1!='-' and b2!='-':
                    if b1==b2:
                        exact_match_num += 1
                    match_num += 1
                total += 1
        #match_num = sum([ 1 for base1,base2 in zip(seq1,seq2) if base1==base2 and base1 in 'ATCG' ])
        #total_num = sum([ 1 for base1,base2 in zip(seq1,seq2) if base1!='-' or base2!='-' ])
        #return match_num/(total_num+0.001)
        return exact_match_num, match_num, total
    
    refSeq = refSeq.replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
    collapsed_id2seq = {}
    for key in id2seq:
        tmp_seq = id2seq[key].replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
        exact_match_num, match_num, total = calc_seq_identity(tmp_seq,refSeq)
        if exact_match_num/total<=max_identity and \
            exact_match_num/(match_num+0.001)>=min_match_identity and \
            (total-match_num)/total<=max_indel_ratio:
            collapsed_id2seq[key] = id2seq[key]
    return collapsed_id2seq

def remove_bpbreak_sequences(id2seq, refStr, maxBpBreakRatio=0.2, maxBpBreakCount=99, maxBpDeletRatio=0.2, maxBpDeletion=99, verbose=False):
    """
    Remove sequences with base pairing violent WC/Wobble rule
    id2seq                      -- {id1:aligned_seq1, id2:aligned_seq2, ...}
    refStr                      -- Reference dot-bracket structure
    maxBpBreakRatio             -- Maximum ratio of violated base pairing (0-1)
                                   A-U => A-C or A-U => A-[-]
    maxBpBreakCount             -- Maximum count of violated base pairing (Integer)
                                   A-U => A-C or A-U => A-[-]
    maxBpDeletRatio             -- Maximum ratio of deleted base pairing (0-1)
                                   A-U => [-]-[-]
    maxBpDeletion               -- Maximum count of deleted base pairing (Integer)
                                   A-U => [-]-[-]
    
    Return: 
        {id1:aligned_seq1, id2:aligned_seq2, ...}           
    """
    import Structure
    def count_bpbreak(sequence, ctList):
        bpbreakCount = 0
        bpDeletCount = 0
        for b1,b2 in ctList:
            if sequence[b1-1]+sequence[b2-1] == '--':
                bpDeletCount += 1
            elif sequence[b1-1]+sequence[b2-1] not in ('AT','TA','CG','GC','GT','TG'):
                bpbreakCount += 1
        return bpbreakCount, bpDeletCount
    ctList = Structure.dot2ct(refStr)
    clean_id2seq = {}
    for key in id2seq:
        tmp_seq = id2seq[key].replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
        bpbreakCount, bpDeletCount = count_bpbreak(tmp_seq,ctList)
        if bpbreakCount<=maxBpBreakCount and \
            bpbreakCount/(len(ctList)+0.001)<=maxBpBreakRatio and \
            bpDeletCount/(len(ctList)+0.001)<=maxBpDeletRatio and \
            bpDeletCount<=maxBpDeletion:
            clean_id2seq[key] = id2seq[key]
        elif verbose:
            print(tmp_seq)
            print(refStr, bpbreakCount, bpDeletCount)
    return clean_id2seq

def remove_gap_columns(id2seq, refSeq=None, refStr=None, minGapRatio=1.0):
    """
    Remove alignment columns with too many gaps
    
    minGapRatio                 -- Remove columns with gap ratio greater than minGapRatio [0-1]
    
    Return {
        'id2seq': {id1:aligned_seq1, id2:aligned_seq2, ...},
        'refSeq': refSeq, # if refSeq != None
        'refStr': refStr  # if refStr != None
    }
    """
    keys = list(id2seq.keys())
    alignment_list = []
    for key in keys:
        alignment_list.append( id2seq[key].replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T') )
    align_columns = collect_columns(alignment_list)
    Len = len(id2seq[keys[0]])
    if refStr is None:
        refBpmap = {}
    else:
        refBpmap = Structure.dot2bpmap(refStr)
    
    columns_to_remove = []
    for i in range(Len):
        if i in columns_to_remove:
            continue
        #print( i, align_columns[i].count('-'), len(align_columns[i]) )
        if align_columns[i].count('-')/len(align_columns[i]) >= minGapRatio:
            columns_to_remove.append( i )
            if i+1 in refBpmap:
                columns_to_remove.append( refBpmap[i+1]-1 )
    
    columns_to_remove.sort()
    #print(columns_to_remove)
    clean_id2seq = {}
    for key in id2seq:
        clean_id2seq[key] = ""
        cur_seq = id2seq[key]
        for i in range(Len):
            if not General.bi_search(i, columns_to_remove):
                clean_id2seq[key] += cur_seq[i]
    
    return_value = { 'id2seq': clean_id2seq }
    
    if refSeq:
        _refSeq = ""
        for i in range(Len):
            if not General.bi_search(i, columns_to_remove):
                _refSeq += refSeq[i]
        return_value['refSeq'] = _refSeq
    
    if refStr:
        _refStr = ""
        for i in range(Len):
            if not General.bi_search(i, columns_to_remove):
                _refStr += refStr[i]
        _Len = len(_refStr)
        #_refStr = Structure.ct2dot(Structure.dot2ct(_refStr), _Len)
        return_value['refStr'] = _refStr
    
    return return_value

def read_sto_DE(stoFn, remove_slash=False):
    """
    Read =GS DE annotation from stoFn
    
    Return:
        { id1: description1, ... }
    """
    seq_desc = {}
    for line in open(stoFn):
        if line.startswith('#=GS'):
            data = line.strip().split()
            tag = data[0]
            flag = data[2]
            if tag=='#=GS' and flag=='DE':
                key = data[1]
                if remove_slash:
                    key = key.split('/')[0]
                seq_desc[key] = seq_desc.get(key, "") + " ".join(data[3:])
    return seq_desc

def get_alignedPos2cleanPos_dict(aligned_seq):
    """
    Return a dictionary of aligned position to clean position
    aligned_seq                 -- Aligned sequence
    
    Return:
        { aligned_pos: clean_pos, ... }
    """
    alignedPos2cleanPos = {}
    i, j = 0, 0
    while i<len(aligned_seq):
        if aligned_seq[i] in ('~',':','.','-'):
            alignedPos2cleanPos[i+1] = None
        else:
            alignedPos2cleanPos[i+1] = j+1
            j += 1
        i += 1
    return alignedPos2cleanPos

def call_covariation(query_seq, query_dot, model_name, seqdbFn, workdir=None,
    nohmm=False, cmsearchE=1, cpu=20, use_LSF=True, 
    LSF_parameters={}, progress=True, clean=False):
    """
    Call covariation has two steps: 
        1. Search homology sequence from sequence databse;
        2. Call covaring base pairs
    
    query_seq               -- Query sequence
    query_dot               -- Query dot-bracket structure
    model_name              -- CM model name and R-scape output file name
    seqdbFn                 -- File name of sequence database
    workdir                 -- Working directory, if not provide a random directory will generate.
    nohmm                   -- cmsearch parameter
    cmsearchE               -- cmsearch E-value
    cpu                     -- Threads for cmcalibrate and cmsearch
    use_LSF                 -- Submit to LSF if True
    LSF_parameters          -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmsearch', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    progress                -- Print the progress
    clean                   -- Clean the directory after running
    
    Return a list of covaring base pairs: 
        [ (left, right), ... ]
    """
    import os, General, Colors, shutil
    
    if workdir is None:
        randID = random.randint(1000000,9000000)
        dirname = f"call_covariation_{randID}"
        workdir = os.path.join(os.environ['HOME'], dirname)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    
    new_seqdbFn = os.path.join(workdir, "seqDB.fa")
    sto_file = os.path.join(workdir, "input.sto")
    cm_file = os.path.join(workdir, "model.cm")
    output_sto = os.path.join(workdir, "output.sto")
    output_txt = os.path.join(workdir, "output.txt")
    R_scape_dir = os.path.join(workdir, "R-scape")
    
    ### Step 0. Prepare combined Fasta file
    if progress:
        print(Colors.f(f"The work directory is: {workdir}", fc='green'))
        print(Colors.f("Step 0. Prepare combined Fasta file", fc='green'))
    General.write_fasta({'input': query_seq}, new_seqdbFn)
    os.system(f'cat {seqdbFn} >> {new_seqdbFn}')
    
    ### Step 1. Build stockholm file
    if progress:
        print(Colors.f("Step 1. Build stockholm file", fc='green'))
    dot2sto({'input': [query_seq, query_dot]}, model_name, sto_file, mode='w')
    
    ### Step 2. Build cm file from sto file
    if progress:
        print(Colors.f("Step 2. Build cm file from sto file", fc='green'))
    cmbuild(sto_file, cm_file, verbose=False, showCMD=progress)
    
    ### Step 3. Calibrate file
    if progress:
        print(Colors.f("Step 3. Calibrate file", fc='green'))
    job = cmcalibrate(cm_file, cpu=cpu, verbose=False, showCMD=progress, use_LSF=use_LSF, LSF_parameters={})
    job.wait()
    
    ### Step 4. Search with CM
    if progress:
        print(Colors.f("Step 4. Search with CM", fc='green'))
    job = cmsearch(cm_file, new_seqdbFn, output_txt, output_sto, 
        cpu=cpu, toponly=True, nohmm=nohmm, 
        nohmmonly=True, outputE=20, acceptE=cmsearchE, 
        cut_ga=False, rfam=False, glocal=False,
        verbose=False, showCMD=progress, use_LSF=use_LSF, LSF_parameters={})
    job.wait()
    
    #### Step 5. R-scape Runing
    if progress:
        print(Colors.f("Step 5. R-scape Runing", fc='green'))
    R_scape(output_sto, R_scape_dir, outname=model_name, maxIdentity=0.985, minIndentity=0.500, 
        F=0.5, gapthresh=0.5, two_set_test=True, fold=False, acceptE=0.05, nseqmin=5, verbose=False, showCMD=progress)
    
    #### Step 6. Read R-scape result
    rscape_file = os.path.join(R_scape_dir,f'{model_name}.cov')
    covary_bps = []
    if os.path.exists(rscape_file):
        rscape_list = read_RScape_result(rscape_file)
        id2seq, refStr, refAnnot = General.load_stockholm(output_sto)[0]
        input_id = [ key for key in id2seq if key.startswith('input') ][0]
        posDict = get_alignedPos2cleanPos_dict(id2seq[input_id])
        for bp in rscape_list:
            covary_bps.append( (posDict[bp[0]], posDict[bp[1]]) )
            l, r = covary_bps[-1]
            bp_base = query_seq[l-1]+query_seq[r-1]
            if bp_base.replace('U','T') not in ('AT','TA','GC','CG','GT','TG'):
                print(f"Warning: {bp_base} is not a RNA base pair")
    
    if clean:
        shutil.rmtree(workdir)
    
    return covary_bps

def calc_MI(left_align_bases, right_align_bases, only_canonical=True, gap_mode='remove'):
    """
    Calculate the mutual information
    left_align_bases            -- ['A', 'C', 'A', 'T', ...]
    right_align_bases           -- ['T', 'G', 'T', 'G', ...]
    only_canonical              -- Only consider AU,UA,GC,CG,UG,GU pairs
    gap_mode                    -- One of remove, penalty, ignore. [remove] mode will remove all base pairs
                                   contains at least one gap; [penalty] mode will give a negative score for base 
                                   pairs with gap; [ignore] mode will treat the gap as a kind of base
    
    Return mutual information
    If return -1, it means too little bases
    """
    import collections
    import numpy as np
    
    left_align_bases = list(left_align_bases)
    right_align_bases = list(right_align_bases)
    assert gap_mode in ('remove', 'penalty', 'ignore'), "gap_mode should be one of remove, penalty, ignore" 
    assert len(left_align_bases) == len(right_align_bases), "Length should be same"
    i = 0
    while i<len(left_align_bases):
        if gap_mode=="remove" and (left_align_bases[i]=='-' or right_align_bases[i]=='-'):
            del left_align_bases[i]
            del right_align_bases[i]
        else:
            left_align_bases[i] = left_align_bases[i].upper().replace('U','T')
            right_align_bases[i] = right_align_bases[i].upper().replace('U','T')
            i += 1
    
    if len(left_align_bases)<5:
        return -1
    
    f1 = collections.Counter(left_align_bases)
    f2 = collections.Counter(right_align_bases)
    
    sum_1 = sum(f1.values())
    sum_2 = sum(f2.values())
    for key in f1:
        f1[key] /= sum_1
    for key in f2:
        f2[key] /= sum_2
    
    fcomb = {}
    for b1,b2 in zip(left_align_bases, right_align_bases):
        fcomb[b1+b2] = fcomb.get(b1+b2, 0) + 1
    
    sum_comb = sum(fcomb.values())
    for key in fcomb:
        fcomb[key] /= sum_comb
    
    MI = 0
    gap_cols = 0
    for key in fcomb:
        if only_canonical and key not in ('AT','TA','GC','CG','GT','TG'):
            continue
        if gap_mode=='penalty' and '-' in key:
            gap_cols += fcomb[key]
            MI -= fcomb[key]
        else:
            b1 = key[0]
            b2 = key[1]
            MI += fcomb[key] * np.log2( (fcomb[key]) / (f1[b1] * f2[b2]) )
    
    return MI

def calc_RNAalignfold(left_align_bases, right_align_bases):
    """
    Calculate the RNAalignfold covariation score
    left_align_bases            -- ['A', 'C', 'A', 'T', ...]
    right_align_bases           -- ['T', 'G', 'T', 'G', ...]
    
    Return the RNAalignfold covariation score
    If return -1, it means too little bases
    """
    import collections
    import numpy as np
    
    BPs = ('AT', 'TA', 'GC', 'CG', 'GT', 'TG')
    
    left_align_bases = list(left_align_bases)
    right_align_bases = list(right_align_bases)
    len1, len2 = len(left_align_bases), len(right_align_bases)
    assert len(left_align_bases) == len(right_align_bases), f"{len1}, {len2} Length should be same"
    Len = len(left_align_bases)
    
    i = 0
    while i<len(left_align_bases):
        left_align_bases[i] = left_align_bases[i].upper().replace('U','T')
        right_align_bases[i] = right_align_bases[i].upper().replace('U','T')
        i += 1
    
    tCount = 0
    tScore = 0
    penalty = 0
    for i in range(Len):
        if left_align_bases[i]+right_align_bases[i] not in BPs:
            penalty += 1
        for j in range(1, Len):
            tCount += 1
            if left_align_bases[i]+right_align_bases[i] not in BPs or left_align_bases[j]+right_align_bases[j] not in BPs:
                continue
            cscore = 2
            if left_align_bases[i]==left_align_bases[j]:
                cscore -= 1
            if right_align_bases[i]==right_align_bases[j]:
                cscore -= 1
            tScore += cscore
    score = tScore / tCount
    #print(penalty/Len)
    score -= penalty/Len
    return score

def calc_RNAalignfold_stack(columns, i, j, min_alignment=5):
    """
    columns             -- [['A','T','C',...], ['A','T','C',...], ...]. Each list contains bases in alignment column
    i,j                 -- The columns to calculate the score
    min_alignment       -- Minimun alignment required
    
    Return RNAalignfold_stack covariation score
    If the count of alignment less than min_alignment, then return -1
    """
    Len = len(columns)
    if len(columns[0])<min_alignment:
        return -1
    score = calc_RNAalignfold(columns[i-1], columns[j-1])
    if i==1 or j==Len:
        score += calc_RNAalignfold(columns[i-1+1], columns[j-1-1])
        score /= 2
    else:
        score = 2*score + calc_RNAalignfold(columns[i-1+1], columns[j-1-1]) + calc_RNAalignfold(columns[i-1-1], columns[j-1+1])
        score /= 4
    return score

def collect_columns(alignment_list):
    """
    alignment_list                  -- [alignSeq1, alignSeq2, alignSeq3...]
    
    Return [ ['A','T','C',...], ['A','T','C',...],... ]
    """
    AlignLen = len(alignment_list[0])
    columns = []
    for i in range(AlignLen):
        columns.append([])
    for alignedSeq in alignment_list:
        for i in range(AlignLen):
            columns[i].append(alignedSeq[i])
    # Check
    for i in range(len(columns)):
        columns[i] = tuple(columns[i])
        assert len(columns[i])==len(columns[0])
    return columns

def calc_covBP_from_sto(stoFn, querySeq=None, query_id_pattern=None, full_seq=None, allpair=False, min_seqs=8, min_score=0.4):
    """
    Calculate the RNAalignfold_stack score for all base pairs in stockholm
    Please note that the refseq in stockholm file should be consistent with querySeq you provide
    
    stoFn                   -- Stockholm file
    querySeq                -- Query sequence, the reference sequence (#=GC RF) will aligned with querySeq
                               If no querySeq provided, the ID match query_id_pattern will used as query sequence
    query_id_pattern        -- The query ID mattern, if querySeq is provided, will be ignored
    full_seq                -- If full_seq is provide and querySeq=None, the sequence index by query_id_pattern will 
                               be search in full_seq, the left and right index will be re-indexed
    allpair                 -- Calculate all base pairs, not only pairs in structure
    min_seqs                -- Minimum sequences in the alignments
    min_score               -- Minimum score to record
    
    Return [ [left, right, score],.... ]
    """
    import Structure, re
    id2seq_dict,refStr,refSeq = General.load_stockholm(stoFn)[0]
    if len(id2seq_dict)<min_seqs:
        return []
    columns = collect_columns([value for key,value in id2seq_dict.items()])
    alignLen = len(refSeq)
    if allpair:
        ct = []
        for i in range(alignLen):
            for j in range(i+1,alignLen):
                ct.append([i, j])
    else:
        ct = Structure.dot2ct(refStr)
    covary_bps = []
    for b1,b2 in ct:
        rnaalignscore2 = calc_RNAalignfold_stack(columns, b1, b2)
        if rnaalignscore2>min_score:
            covary_bps.append( [b1, b2, round(rnaalignscore2,3)] )
    refSeq = refSeq.replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
    if querySeq:
        refSeq_realigned, query_aligned = Structure.multi_alignment([refSeq.replace('-',''), querySeq.upper().replace('U','T')])
        ### 构建一个字典：从refSeq的坐标指向querySeq的坐标
        refseq2query = {}
        pos_ref,pos_refrealign = 0,0
        query_pos = 0
        while pos_ref<len(refSeq):
            while pos_ref<len(refSeq) and refSeq[pos_ref]=='-':
                pos_ref += 1
            if pos_ref==len(refSeq): break
            while refSeq_realigned[pos_refrealign]=='-':
                if query_aligned[pos_refrealign]!='-':
                    query_pos += 1
                pos_refrealign += 1
            #if pos_refrealign==len(query_aligned): break
            refseq2query[pos_ref+1] = query_pos+1
            pos_ref += 1
            pos_refrealign += 1
            query_pos += 1
        i = 0
        while i < len(covary_bps):
            if covary_bps[i][0] in refseq2query and covary_bps[i][1] in refseq2query:
                covary_bps[i][0] = refseq2query[covary_bps[i][0]]
                covary_bps[i][1] = refseq2query[covary_bps[i][1]]
                i += 1
            else:
                del covary_bps[i]
    elif query_id_pattern:
        tid_list = [ name for name in id2seq_dict.keys() if re.match(query_id_pattern, name) ]
        assert len(tid_list)==1, f"Error: {len(tid_list)} tid match pattern"
        query_tid = tid_list[0]
        clean_query_seq = id2seq_dict[query_tid].replace('~','').replace(':','').replace('.','').replace('-','').upper().replace('U','T')
        if full_seq:
            full_seq = full_seq.replace('U','T')
            start_pos = full_seq.find(clean_query_seq)
            if start_pos==-1:
                raise RuntimeError('queryseq not in full_seq')
            if full_seq.find(clean_query_seq,start_pos+1)!=-1:
                raise RuntimeError('queryseq are find in multiple position')
        else:
            start_pos = 0
        aligned_query_seq = id2seq_dict[query_tid].upper().replace('U','T')
        alignedPos2cleanPos = get_alignedPos2cleanPos_dict(id2seq_dict[query_tid])
        i = 0
        while i < len(covary_bps):
            l,r = covary_bps[i][:2]
            if alignedPos2cleanPos[l] and alignedPos2cleanPos[r]:
                nl = alignedPos2cleanPos[l] + start_pos
                nr = alignedPos2cleanPos[r] + start_pos
                covary_bps[i][0] = nl
                covary_bps[i][1] = nr
                if full_seq:
                    assert aligned_query_seq[l-1]==full_seq[nl-1], "Exception: seq not same with full_seq"
                    assert aligned_query_seq[r-1]==full_seq[nr-1], "Exception: seq not same with full_seq"
                i += 1
            else:
                del covary_bps[i]
    else:
        raise RuntimeError('querySeq or query_id_pattern must be specified')
    return covary_bps

def Build_Rfam_seed(input_sto, input_seq, seqDBFastaFn, workdir, toponly=True, iteration=3, blastFilter=None, use_LSF=True, cores=5, verbose=True):
    """
    This function build a rfam seed according to https://docs.rfam.org/en/latest/building-families.html
    Repeat this process "iteration" times:
    Seed -[cmbuild]-> CM -[cmcalibrate]-> CM -[cmsearch]-> Stockholm -[Filter]-> Fasta -[cmalign]-> Stockholm -[cmbuild]-> CM
    
    input_sto               -- Seed stockholm file, contains single sequence and structure, output file of Covariation.dot2sto
    input_seq               -- The sequence in input_sto
    seqDBFastaFn            -- The target sequences, contains all homologous sequences
    workdir                 -- Work directory, all immediate files will be contained
    toponly                 -- Only search the positive strand
    iteration               -- How many times to repeat the cmsearch processes
    blastFilter             -- If set, will use blast result to filter the cmsearch target. The cmsearch coordinate will
                               be compared with blast result. The purpose is to prevent the cmsearch give too much weight 
                               on structure and ignore the sequence
                               { 'full_seq': "ACTG...",  # The full length sequence contain input_seq
                                'blastdb':"/path/to/blastdb", 
                                'flanking': 20, # Max distance between the cmsearch result and blast result
                                'blast_identity': 0.4, # Blast cutoff
                                'min_match': 5, # Minimun target
                                'ignore_seq': ['input'] # Which results to ignore
                                }
    use_LSF                 -- bsub or run in local
    cores                   -- How many cores to use
    verbose                 -- If print the information in the screen
    
    Return:
        Return a stockholm file with final alignment
    """
    import Alignment, Colors
    if blastFilter:
        input_start = blastFilter['full_seq'].find(input_seq)
        assert input_start != -1, f"input_seq not in full_seq, consider maintain consistent letter case and U/T"
        input_end = input_start + len(input_seq)
        start = max(input_start-blastFilter.get('flanking', 20), 0)
        end = min(input_end+blastFilter.get('flanking', 20), len(blastFilter['full_seq']))
        hits = Alignment.blast_sequence_V2(blastFilter['full_seq'][start:end], blastFilter['blastdb'], verbose=verbose, task='blastn', maxhit=5000, perc_identity=blastFilter['blast_identity'] )
        target_valid_range = {}
        for hit in hits:
            if hit.hit_acc not in target_valid_range:
                target_valid_range[hit.hit_acc] = (hit.hit_from, hit.hit_to)
        if len(target_valid_range)<blastFilter.get('min_match', 5):
            if verbose: print(Colors.f("Failed: Too little blast target", 'red'))
            return False
        else:
            if verbose: print(f"Blast: {len(target_valid_range)} results found")
    
    for it_num in range(iteration):
        if verbose: print(f"===========>>> Iteration {it_num}  <<<===========")
        
        ### Define files 
        if it_num==0:
            input_sto = input_sto
        else:
            input_sto = os.path.join(workdir, f'cmalign_output_{it_num-1}.sto')
        input_cm = os.path.join(workdir, f'cmsearch_input_{it_num}.cm')
        cmsearch_output_txt = os.path.join(workdir, f'cmsearch_output_{it_num}.txt')
        cmsearch_output_sto = os.path.join(workdir, f'cmsearch_output_{it_num}.sto')
        cmalign_input_fa = os.path.join(workdir, f'cmalign_input_{it_num}.fa')
        cmalign_output_sto = os.path.join(workdir, f'cmalign_output_{it_num}.sto')
        
        ### Build CM model
        cmbuild(input_sto, input_cm)
        h = cmcalibrate(input_cm, use_LSF=use_LSF, LSF_parameters={'cpu': cores})
        if use_LSF: h.wait()
        
        ### cmsearch
        tmp_fa,tmp_annot = General.load_fasta(seqDBFastaFn, load_annotation=True)
        tmp_fa['input'] = input_seq
        General.write_fasta(tmp_fa, os.path.join(workdir, 'tmp_ref.fa'), tmp_annot)
        del tmp_fa
        h = cmsearch(input_cm, os.path.join(workdir, 'tmp_ref.fa'), \
            cmsearch_output_txt, cmsearch_output_sto, cpu=cores, toponly=toponly, nohmm=False, nohmmonly=False, \
            outputE=20, acceptE=1, cut_ga=False, rfam=False, glocal=False, verbose=True, showCMD=True, use_LSF=use_LSF, 
            LSF_parameters={'cpu': cores})
        if use_LSF: h.wait()
        os.remove(os.path.join(workdir, 'tmp_ref.fa'))
        
        ### Read cmsearch result
        id2seq_dict, refStr, refSeq = General.load_stockholm(cmsearch_output_sto)[0]
        refID = [ key for key in id2seq_dict if key.startswith('input') ][0]
        refSeq = id2seq_dict[refID]
        
        ### Filter results not in blast hits
        if blastFilter:
            for key in list(id2seq_dict.keys()):
                if '/' in key:
                    true_id, start_end = key.rsplit('/', 1)
                    start, end = start_end.split('-')
                    start = int(start)
                    end = int(end)
                    if true_id in target_valid_range:
                        valid_start, valid_end = target_valid_range[true_id]
                        if (end<valid_start or start>valid_end):
                            del id2seq_dict[key]
                    elif true_id in blastFilter.get('ignore_seq',['input']):
                        if verbose: print("Ignored:",true_id)
                    else:
                        if verbose: print("Not in blast hits:",Colors.f(true_id, fc='red'))
                        del id2seq_dict[key]
                else:
                    if verbose: print("Not slash in key:",Colors.f(true_id, fc='green'))
                    del id2seq_dict[key]
        
        ### Remove highly similar sequences
        uniq_seq = collapse_sequences(id2seq_dict, refSeq, max_identity=0.98, min_match_identity=0.2, max_indel_ratio=0.5)
        ### Remove sequences with too many base pairs broken
        clean_id2seq = remove_bpbreak_sequences(uniq_seq, refStr, maxBpBreakRatio=0.2, maxBpDeletRatio=0.5, maxBpBreakCount=99, maxBpDeletion=99)
        ### Remove sequences with nucleotides other than AUCG
        for key in list(clean_id2seq.keys()):
            v = clean_id2seq[key].replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
            if len(v) != v.count('A')+v.count('T')+v.count('C')+v.count('G')+v.count('-'):
                del clean_id2seq[key]
        pure_fasta = { k:v.replace('~','').replace(':','').replace('.','').replace('-','').upper().replace('U','T') for k,v in clean_id2seq.items() }
        ### Remove highly similar sequences
        collapsed_pure_fasta = Alignment.cd_hit_est(pure_fasta, identity=0.98, verbose=True)
        clean_id2seq = { k:clean_id2seq[k] for k in collapsed_pure_fasta }
        if len(clean_id2seq)==0:
            #### No homologous sequence to call covariation in the sequence db
            if verbose: print(Colors.f("Failed: No homologous sequences left", 'red'))
            return False
        ### Remove alignment columns with too many gaps
        returnvalues = remove_gap_columns(clean_id2seq, refSeq=refSeq, refStr=refStr, minGapRatio=1.0)
        
        ### cmalign
        GS_DE = read_sto_DE(cmsearch_output_sto) # read description
        small_seqdb = { k:v.replace('-','') for k,v in returnvalues['id2seq'].items() }
        small_seqdb['input'] = input_seq
        General.write_fasta(small_seqdb, os.path.join(workdir, 'cmalign_input.fa'), GS_DE)
        h = cmalign(input_cm, os.path.join(workdir, 'cmalign_input.fa'), cmalign_output_sto, cpu=cores, glocal=False, 
            outformat='Stockholm', mxsize=1028.0, verbose=True, showCMD=True, use_LSF=use_LSF, LSF_parameters={'cpu': cores})
        if use_LSF: h.wait()
    
    return cmalign_output_sto


    