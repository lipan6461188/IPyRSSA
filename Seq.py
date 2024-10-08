#-*- coding:utf-8 -*-

import os, sys, time, re, random, pickle, copy, gzip, io, yaml, logging, configparser, math, shutil, pathlib, tempfile, hashlib, argparse, json, inspect, urllib, collections, subprocess, requests, platform, multiprocessing, importlib, string, code, warnings, concurrent, gc, functools, types, traceback, base64, bz2, ctypes
from pyliftover import LiftOver
from tqdm.auto import tqdm, trange
from os.path import exists, join, getsize, isfile, isdir, abspath, basename, realpath, dirname
from . import General

# def reverse_comp(sequence):
#     """
#     sequence                -- Raw input sequence
    
#     Return the reverse complimentary sequence
#     """
#     RC_Map = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N',
#                'a':'t', 'c':'g', 't':'a', 'g':'c', '-':'-'}

#     return ''.join(map(lambda x: RC_Map[x], list(sequence[::-1])))

def reverse_comp(sequence):
    """
    sequence                -- Raw input sequence
    
    Return the reverse complimentary sequence
    """
    RC_Map = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N',
               'a':'t', 'c':'g', 't':'a', 'g':'c', '-':'-'}
    
    rc_seq = []
    for i in range(len(sequence)-1, -1, -1):
        rc_seq.append( RC_Map[ sequence[i] ] )
    rc_seq = ''.join(rc_seq)
    
    return rc_seq

def flat_seq(sequence, lineLen=60):
    """
    sequence                -- Raw input sequence
    lineLen                 -- Number of nucleotides for each line
    
    Flat raw long sequence to multiple lines
    """
    flatted_seq = ''
    sequence = sequence.strip()
    
    idx = 0
    while idx < len(sequence):
        flatted_seq += sequence[idx:idx+lineLen]+'\n'
        idx += lineLen
    
    return flatted_seq[:-1]

def format_gene_type(gene_type):
    """
    gene_type               -- Raw gene type
    """
    import re
    
    valid_gene_type = ('pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA', 'mRNA', 'tRNA')
    lncRNA_class = ('3prime_overlapping_ncrna','antisense','lincRNA','non_coding','sense_intronic','sense_overlapping','processed_transcript','lnc_RNA')
    
    if gene_type in valid_gene_type:
        return gene_type;
    
    if re.match('.*pseudogene',gene_type):
        return 'pseudogene';
    
    if gene_type == 'protein_coding':
        return 'mRNA';
    
    if gene_type in lncRNA_class:
        return 'lncRNA';
    
    return 'other'

class seqClass(object):
    def __init__(self, seqFn):
        import pysam
        
        self.genome = pysam.Fastafile(seqFn)
        sys.stdout.writelines("seqClass: input 0-based coordinate -- [start, end)\n")
        self.size_dict = self._seq_len_dict()
        self.key_list = sorted(list(self.size_dict.keys()))

    def fetch(self, chrID, chrStart=None, chrEnd=None, chrStrand="+"):
        """
        chrID           -- Chromosome ID
        chrStart        -- Chromosome start position
        chrEnd          -- Chromosome end position
        chrStrand       -- + or -
        """
        if chrStart is None:
            chrStart = 0
        if chrEnd is None:
            chrEnd = self.size_dict[chrID]

        if chrStrand == '+':
            return self.genome.fetch(chrID, chrStart, chrEnd)
        
        if chrStrand == '-':
            return reverse_comp(self.genome.fetch(chrID, chrStart, chrEnd))
    
    def get(self, chrID):
        """
        Get the full sequence
        """
        return self.genome.fetch(chrID, 0, self.size_dict[chrID])

    def has(self, chrID):
        return chrID in self.genome.references
    
    def _seq_len_dict(self):
        size_dict = {}
        for chr_id, chr_len in zip(self.genome.references, self.genome.lengths):
            size_dict[chr_id] = chr_len
        return size_dict

    def seq_len_dict(self):
        return self.size_dict

class liftOverClass(LiftOver):
    def __init__(self, from_db, to_db):
        """
        from_db         -- 'hg19','hg38','mm9','mm10'
        to_db           -- 'hg19','hg38','mm9','mm10'
        """
        from pyliftover import LiftOver
        
        LiftOver.__init__(self, from_db=from_db, to_db=to_db)
    
    def convert_coor(self, Chr, Pos, Strand='+'):
        """
        Chr             -- Chromosome id
        Pos             -- Genome position
        
        object.convert_coor('chr1', 109303388, '+')
        """
        return self.convert_coordinate(chromosome=Chr, position=Pos, strand=Strand)

def lift_genome(lifter, chrID, chrStart, chrEnd, chrStrand, verbose=False):
    """
    lifter          -- An object of pyliftover.LiftOver
    chrID           -- Chromosome ID
    chrStart        -- Chromosome start position
    chrEnd          -- Chromosome end position
    chrStrand       -- + or -
    verbose         -- Show warning information
    
    Convert genome position in different genome version
    
    Return (chrID, chrStart, chrEnd, chrStrand) if success
    Return (-1, -1, -1, -1) if failed
    """
    start_list = lifter.convert_coordinate(chrID, chrStart, chrStrand)
    end_list = lifter.convert_coordinate(chrID, chrEnd, chrStrand)
    if start_list is None or end_list is None:
        if verbose:
            sys.stderr.writelines("Warning: %s:%s-%s(%s) cannot convert -- chromosome not found\n" % (chrID, chrStart, chrEnd, chrStrand))
        return (-1, -1, -1, -1)
    
    if len(start_list) == 0 or len(end_list) == 0:
        if verbose:
            sys.stderr.writelines("Warning: %s:%s-%s(%s) cannot convert -- not appropriate map position\n" % (chrID, chrStart, chrEnd, chrStrand))
        return (-1, -1, -1, -1)
    
    s_chrID, s_chrPos, s_chrStrand, _ = start_list[0]
    e_chrID, e_chrPos, e_chrStrand, _ = end_list[0]
    
    if s_chrID != e_chrID or s_chrStrand != e_chrStrand:
        if verbose:
            sys.stderr.writelines("Warning: %s:%s-%s(%s) cannot convert -- different start/end chromosome or strand\n" % (chrID, chrStart, chrEnd, chrStrand))
        return (-1, -1, -1, -1)
    
    if s_chrID != chrID or s_chrStrand != chrStrand:
        if verbose:
            sys.stderr.writelines("Warning: %s:%s-%s(%s) cannot convert -- chromosome or strand changed\n" % (chrID, chrStart, chrEnd, chrStrand))
        return (-1, -1, -1, -1)
    
    return (s_chrID, s_chrPos, e_chrPos, s_chrStrand)

def search_subseq_from_genome(Seqer, chrID, start, end, strand, pattern, caller='first', check=False):
    """
    Parameters:
        Seqer: Object of Seq.seqClass()
        chrID: Chromosome ID
        start, end: The genome start pos and end pos. [start, end)
        strand: + or -
        pattern: Sequence or regex
        caller: first, all or callable
                first will return the first result, 
                all will return all result,
                callable function will be called
                    caller(sequence, pattern) must return [ [start, end], [start, end], ... ] or [] if no result
        check: bool. Check the result
    Return:
        if caller is callable or all: [ [genome_start, genome_end], ... ]
        if caller is first: [genome_start, genome_end]
    """
    import re
    
    subseq = Seqer.fetch(chrID, start, end, strand)
    if callable(caller):
        hits = caller(subseq, pattern)
    elif isinstance(caller, str):
        hits = list(re.finditer( pattern,  subseq ))
        if len(hits)>0:
            hits = [ hit.span() for hit in hits ]
    else:
        raise RuntimeError("caller must be first, all or callable")
    
    if strand == '+':
        hits = [ (start+hit[0], start+hit[1]) for hit in hits ]
    else:
        hits = [ (end-hit[0]-(hit[1]-hit[0]), end-hit[0]) for hit in hits ]
    
    if check:
        if len(hits)>0:
            print(f"Check {chrID}:{start}-{end}({strand}) {pattern}")
        for hit in hits:
            new_subseq = Seqer.fetch(chrID, hit[0], hit[1], strand)
            pass_ = not (re.match(pattern, new_subseq) is None)
            print(f"{new_subseq} pass {pass_}")
    
    if caller == 'first':
        hits = hits[0]
    
    return hits


###################################
### Fasta Reader
###################################

class FastaReader:
    def __init__(self, fasta_file, index_file=None, disable_tqdm=False):
        
        self.fasta_file = fasta_file
        self.index_file = index_file
        self.index_file = self.build_index(disable_tqdm=disable_tqdm)
        index           = self.read_fasta_index(disable_tqdm=disable_tqdm)
        self.names      = index[0]
        self.ranges     = index[1]
    
    def build_index(self, disable_tqdm=False):
        """
        Build index file for Fasta
        
        Parameters
        -----------------
        disable_tqdm: Disable tqdm progress bar
        
        Return
        ----------------
        index_file
        """
        
        if self.index_file is None:
            self.index_file = self.fasta_file + '.index'
        
        lc = int(subprocess.getoutput(f'wc -l {self.fasta_file}').split()[0])
        if os.path.exists(self.index_file):
            return self.index_file
        
        with open(self.index_file, 'w') as OUT, open(self.fasta_file) as IN:
            byte_start = byte = 0
            for line in tqdm(IN, dynamic_ncols=True, total=lc, disable=disable_tqdm, desc='Build Fasta Index'):
                if line.startswith('>'):
                    if byte_start < byte:
                        print(header, byte_start, byte, sep='\t', file=OUT, flush=True)
                    header = line[1:].split()[0]
                    byte_start = byte
                byte += len(line)
            if byte_start < byte:
                print(header, byte_start, byte, sep='\t', file=OUT, flush=True)
        
        cmd = f"sort --parallel={os.cpu_count()//2} -S 10G -k1,1 {self.index_file} > {self.index_file}.sorted"
        status, output = subprocess.getstatusoutput(cmd)
        if status != 0:
            print(f"Error in sorting: {output}")
            return None
        else:
            _ = os.system(f"mv {self.index_file}.sorted {self.index_file}")
        
        return self.index_file
    
    def read_fasta_index(self, disable_tqdm=False):
        lc = int(subprocess.getoutput(f'wc -l {self.index_file}').split()[0])
        names = []
        ranges = []
        with open(self.index_file) as IN:
            for line in tqdm(IN, total=lc, dynamic_ncols=True, desc='Read Fasta Index'):
                name, start, end = line.strip().split()
                if len(names) > 0:
                    assert name > names[-1], f"Expect header be sorted. but got {names[-1]} and {name}"
                names.append(name)
                ranges.append((int(start), int(end)))
        return (names, ranges)
    
    def get(self, name=None, i=None, return_all=False):
        """
        Get sequence
        
        Parameters
        -----------------------
        name: str
        return_all: bool. return (header, annot, seq) or only seq
        """
        if name is not None:
            # (names, ranges) = (self.names, self.ranges)
            i = General.bi_search(name, self.names, retern_index=True)
            assert i >= 0, f"Expect {name} in names, but not found"
            start, end = self.ranges[i]
        elif i is not None:
            start, end = self.ranges[i]
        else:
            raise RuntimeError(f"One of name and i must be specified")
        
        with open(self.fasta_file) as IN:
            _ = IN.seek(start)
            content = IN.read(end - start)
        
        lines  = content.strip().split('\n')
        header = lines[0][1:].split()[0]
        annot  = lines[0][1+len(header)+1:]
        seq    = "".join(lines[1:])
        
        if return_all:
            return (header, annot, seq)
        else:
            return seq
    
    def __len__(self):
        return len(self.names)
    
    def __getitem__(self, i):
        if isinstance(i, int):
            return self.get(i=i, return_all=True)
        elif isinstance(i, str):
            return self.get(name=i)
        else:
            raise RuntimeError(f"Expect i be int or str, but got {type(i)}")
    
    def __repr__(self):
        return f"FastaReader object with {len(self.names)} sequences"

###################################
### Clustering sequences
###################################

def cluster_sequences(seq_dict, method='cd-hit', seq_type='prot', id_cutoff=0.80, cov_cutoff=0.5, verbose=False):
    """
    Parameters
    ------------
    seq_dict: {seq_name: 'AXJIJDIWIJDIWD', ...}
    method: cd-hit or mmseqs
    seq_type: prot or dna
    id_cutoff: float
        Similarity id_cutoff
    verbose: bool
        Print the command
    
    Return
    -----------
    cluster: [ [seq_id, seq_length, similarity], ... ]
    """
    assert seq_type in ('prot', )
    assert method in ('cd-hit', 'mmseqs')
    import tempfile, shutil, subprocess
    
    workdir = tempfile.mkdtemp(prefix="clusterSeq_")
    
    if method == 'cd-hit':
        cdhit_bin = shutil.which('cd-hit')
        assert cdhit_bin is not None, "Error: cd-hit not in PATH"
        
        input_fa = os.path.join(workdir, "input.fa")
        output_fa = os.path.join(workdir, "output.fa")
        clust_fn = os.path.join(workdir, "output.fa.clstr")
        General.write_fasta(seq_dict, input_fa)
        
        cmd = f"{cdhit_bin} -i {input_fa} -o {output_fa} -c {id_cutoff} -aL {cov_cutoff} -aS {cov_cutoff}"
        if verbose:
            print(cmd)
        _ = subprocess.getoutput(cmd)
        if not os.path.exists(clust_fn):
            print(f"Run cd-hit failed, {clust_fn} not exists")
            cluster = None
        else:
            cluster = read_cdhit_clstr(clust_fn)
    else:
        mmseqs_bin = shutil.which('mmseqs')
        assert mmseqs_bin is not None, "Error: mmseqs not in PATH"
        
        input_fa    = os.path.join(workdir, "input.fa")
        output_pref = os.path.join(workdir, "output")
        tmp_dir     = os.path.join(workdir, "tmp")
        General.write_fasta(seq_dict, input_fa)
    
        cmd = f"{mmseqs_bin} easy-cluster --cov-mode 0 --threads 8 --min-seq-id {id_cutoff} -c {cov_cutoff} {input_fa} {output_pref} {tmp_dir}"
        if verbose:
            print(cmd)
        _ = subprocess.getoutput(cmd)
        if not os.path.exists(f"{output_pref}_cluster.tsv"):
            print(f"Run mmseqs failed, {output_pref}_cluster.tsv not exists")
            cluster = None
        else:
            cluster = {}
            for line in open(f"{output_pref}_cluster.tsv"):
                a, b = line.strip().split()
                cluster[a] = cluster.get(a, []) + [b]
            cluster = list(cluster.values())
    
    shutil.rmtree(workdir)
    return cluster

def read_cdhit_clstr(clstr_fn):
    """
    Read sequence clusters from CD-Hit result
    """
    cluster = []
    cur_cluster = []
    for line in open(clstr_fn):
        if line[0] == '>':
            if len(cur_cluster) > 0:
                cluster.append(cur_cluster)
                cur_cluster = []
        else:
            length, seq_id, similarity = re.findall(r"(\d+)aa, >(.*)\.\.\. (.*)$", line)[0]
            length = int(length)
            if similarity == '*':
                similarity = 100.0
            else:
                similarity = float(similarity[3:-1])
            cur_cluster.append( [seq_id, length, similarity] )
    if len(cur_cluster) > 0:
        cluster.append(cur_cluster)
        cur_cluster = []
    return cluster






