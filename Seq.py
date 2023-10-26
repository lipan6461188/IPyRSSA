#-*- coding:utf-8 -*-

import re
from pyliftover import LiftOver
import sys, os
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
### Clustering sequences
###################################

def cluster_sequences(seq_dict, seq_type='prot', cutoff=0.80, verbose=False):
    """
    Parameters
    ------------
    seq_dict: {seq_name: 'AXJIJDIWIJDIWD', ...}
    seq_type: prot or dna
    cutoff: float
        Similarity cutoff
    verbose: bool
        Print the command
    
    Return
    -----------
    cluster: [ [seq_id, seq_length, similarity], ... ]
    """
    assert seq_type in ('prot', )
    import tempfile, shutil, subprocess
    
    cdhit_bin = shutil.which('cd-hit')
    assert cdhit_bin is not None, "Error: cd-hit not in PATH"
    
    workdir = tempfile.mkdtemp(prefix="clusterSeq_")
    input_fa = os.path.join(workdir, "input.fa")
    output_fa = os.path.join(workdir, "output.fa")
    clust_fn = os.path.join(workdir, "output.fa.clstr")
    General.write_fasta(seq_dict, input_fa)
    
    cmd = f"{cdhit_bin} -i {input_fa} -o {output_fa} -c {cutoff}"
    if verbose:
        print(cmd)
    _ = subprocess.getoutput(cmd)
    if not os.path.exists(clust_fn):
        print(f"Run cd-hit failed, {clust_fn} not exists")
        cluster = None
    else:
        cluster = read_cdhit_clstr(clust_fn)
    
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






