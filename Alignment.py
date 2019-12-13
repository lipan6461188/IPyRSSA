#-*- coding:utf-8 -*-

import os, General, random

class BlastNHit(object):
    def __init__(self):
        self.query_id = ""
        self.subject_id = ""
        self.identity = 100
        self.align_len = 0
        self.mismatch = 0
        self.gap_opens = 0
        self.query_start = 0
        self.query_end = 0
        self.subject_start = 0
        self.subject_end = 0
        self.evalue = 0
        self.bit_score = 0
        self.strand = "+"
    def set_sub_range(self, start, end):
        if end<start:
            self.strand = "-"
            start, end = end, start
        self.subject_start, self.subject_end = start, end
    def __repr__(self):
        return "(%s [%d-%d] <==> %s [%d-%d](%c) %.2f%%)" % ( \
            self.query_id, self.query_start, self.query_end, \
            self.subject_id, self.subject_start, self.subject_end, self.strand, \
            self.identity)
    def __str__(self):
        return "(%s [%d-%d] <==> %s [%d-%d](%c) %.2f%%)" % ( \
            self.query_id, self.query_start, self.query_end, \
            self.subject_id, self.subject_start, self.subject_end, self.strand, \
            self.identity)

def blast_seq(query_seq, blastdb, clear=True, verbose=False, perc_identity=90, threads=1):
    """
    Blastn a sequence or sequences and return a list of BlastNHit
    
    query_seq               -- Query sequence, a single seq (str) or a dict of sequences
    blastdb                 -- Balstn db, such as /150T/zhangqf/GenomeAnnotation/INDEX/blast/hg38/hg38
    clear                   -- Clear tempropary files
    verbose                 -- Output command
    perc_identity           -- Minimum percentage of identity
    threads                 -- How many threads to use
    
    Return: [ BlastNHit, BlastNHit,... ]
    
    Require: blastn
    """
    
    blastn = General.require_exec("blastn", "blastn is required")
    if not os.path.exists(blastdb+".nhr"):
        raise RuntimeError("Error: %s.nhr does not exists" % (blastdb))
    if type(query_seq) is not str and type(query_seq) is not dict:
        raise RuntimeError("Error: parameter query_seq should a sequence or a dict of sequence")
    
    rand_id = random.randint(10000,99999)
    seq_fa_file = "/tmp/blast_seq_%s.fa" % (rand_id)
    result_file = "/tmp/blast_result_%s.tabular" % (rand_id)
    
    if type(query_seq) is str:
        General.write_fasta({"query_id":query_seq}, seq_fa_file)
    else:
        General.write_fasta(query_seq, seq_fa_file)
    cmd = "%s -db %s -query %s -out %s -outfmt 7 -num_threads %d -perc_identity %d" % (blastn, blastdb, seq_fa_file, result_file, threads, perc_identity)
    
    if verbose:
        print(cmd)
    os.system(cmd)
    
    Hit_List = []
    for line in open(result_file):
        if line[0]!='#':
            data = line.strip().split()
            hit = BlastNHit()
            hit.query_id = data[0]
            hit.subject_id = data[1]
            hit.identity = float(data[2])
            hit.align_len = int(data[3])
            hit.mismatch = int(data[4])
            hit.gap_opens = int(data[5])
            hit.query_start = int(data[6])
            hit.query_end = int(data[7])
            hit.set_sub_range( int(data[8]), int(data[9]) )
            hit.evalue = float(data[10])
            hit.bit_score = float(data[11])
            Hit_List.append(hit)
    
    if clear:
        os.remove(seq_fa_file)
        os.remove(result_file)
    
    return Hit_List

def annotate_seq(Fasta, genome_blastn_db, Gaper, verbose=True, threads=10):
    """
    Given sequence, blast in genome blastn db and annotate
    
    Fasta                   -- Sequence dict, { tid1:seq1, tid2:seq2, ... }
    genome_blastn_db        -- Balstn db, such as /150T/zhangqf/GenomeAnnotation/INDEX/blast/hg38/hg38
    Gaper                   -- A object of GAP.init()
    verbose                 -- Output information during processing
    threads                 -- How many threads to use
    
    Return: {raw_id=>new_id}, unmaped_ids, unannot_ids
    
    Require: blastn
    """
    
    id_map = {}
    unmapped_ids = []
    unannot_ids = {}
    
    Hits = blast_seq(Fasta, genome_blastn_db, threads=threads)
    Hit_dict = {}
    
    Hit_dict = {}
    for hit in Hits:
        qid = hit.query_id
        try:
            Hit_dict[qid].append(hit)
        except KeyError:
            Hit_dict[qid] = [hit]
    
    for qid in Hit_dict:
        Hit_dict[qid].sort(key=lambda x: x.identity, reverse=True)
    
    for qid in Fasta:
        seq = Fasta[qid]
        if verbose: print("ID:", qid, "Len:", len(seq), sep=" ", end="\t")
        #Hit_list = blast_single_seq(seq, genome_blastn_db, threads=threads)
        if qid not in Hit_dict:
            if verbose: print(" has 0 alignment")
            unmapped_ids.append(qid)
            continue
        
        Hit_list = Hit_dict[qid]
        i = 0
        while i<len(Hit_list):
            hit = Hit_list[i]
            if hit.identity < 90: continue
            trans_loc = Gaper.genomeCoor2transCoor(hit.subject_id, hit.subject_start, hit.subject_end, hit.strand)
            if len(trans_loc)==0: 
                i += 1
            else:
                break
        if i==len(Hit_list):
            if verbose: print(" has no annotation: %s"%(Hit_list[0]))
            hit = Hit_list[0]
            unannot_ids[qid] = "%s:%d-%d(%c)" % (hit.subject_id, hit.subject_start, hit.subject_end, hit.strand)
            continue
        
        gencode_id = trans_loc[0][3]
        trans_region = trans_loc[0][4:6]
        ft = Gaper.getTransFeature(gencode_id)
        annot_name = ft['gene_type']+"-"+ft['gene_name']+"-"+gencode_id+"_%d-%d" % (trans_region[0], trans_region[1])
        if verbose: print("Found", len(Hit_list), "hits", "gene_type:", ft['gene_type'], sep=" ")
        
        id_map[qid] = annot_name
    
    return id_map, unmapped_ids, unannot_ids

