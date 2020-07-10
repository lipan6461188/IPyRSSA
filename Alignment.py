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

def blast_seq(query_seq, blastdb, clear=True, verbose=False, perc_identity=90, evalue=10, maxhit=500, threads=1):
    """
    Blastn a sequence or sequences and return a list of BlastNHit
    
    query_seq               -- Query sequence, a single seq (str) or a dict of sequences
    blastdb                 -- Balstn db, such as /150T/zhangqf/GenomeAnnotation/INDEX/blast/hg38/hg38
    clear                   -- Clear tempropary files
    verbose                 -- Output command
    perc_identity           -- Minimum percentage of identity
    evalue                  -- Maximun evalue
    maxhit                  -- Maximum number of aligned sequences to keep
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
    cmd = f"{blastn} -db {blastdb} -query {seq_fa_file} -out {result_file} -outfmt 7 -num_threads {threads} -perc_identity {perc_identity} -evalue {evalue} -max_target_seqs {maxhit}"
    
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

class AlignedFasta(object):
    def __init__(self, afa_fn, gap_sym="-", verbose=1):
        self.gap_sym = gap_sym
        self.fasta_dict = General.load_fasta(afa_fn)
        self.seq_keys = list(self.fasta_dict.keys())
        if self._check_seq_num() is False:
            raise RuntimeError("Sequence counte should > 1")
        if self._check_len() is False:
            raise RuntimeError("Sequences in file should have same length")
        self.alignLen = len(self.fasta_dict[self.seq_keys[0]])
        if verbose:
            print(f"Total {len(self.seq_keys)} sequences, aligned length: {self.alignLen}")
    
    def _check_seq_num(self):
        if len(self.seq_keys) <= 1:
            return False
        return True
    
    def _check_len(self):
        len_list = [ len(v) for v in self.fasta_dict.values() ]
        if len(set(len_list))!=1:
            return False
        return True
    
    def alignPos2seqPos(self, alignPos, seqID):
        """
        Covert the global aligned position to sequence position
        alignPos                -- Position of alignment. 1-Base
        seqID                   -- Sequence ID
        
        Return:
            Sequence position. 1-Based
        """
        assert 1<=alignPos<=self.alignLen
        assert seqID in self.seq_keys
        pos = 0
        while alignPos != 0:
            if self.fasta_dict[seqID][alignPos-1]!=self.gap_sym:
                pos += 1
            alignPos -= 1
        return pos
    
    def seqPos2alignPos(self, seqPos, seqID):
        """
        Covert the sequence position to global aligned position
        seqPos                  -- Sequence position. 1-Based
        seqID                   -- Sequence ID
        
        Return:
            Position of alignment. 1-Base
        """
        assert 1<=seqPos<=self.alignLen
        assert seqID in self.seq_keys
        pos = 0
        while seqPos != 0:
            pos += 1
            if self.fasta_dict[seqID][pos-1]!=self.gap_sym:
                seqPos -= 1
        return pos
    
    def seqPos2seqPos(self, seqPos, querySeqID, targetSeqID):
        """
        Covert the sequence position to global aligned position
        seqPos                -- Sequence position. 1-Based
        querySeqID            -- Query Sequence ID
        targetSeqID           -- Target Sequence ID
        
        Return:
            Position of alignment. 1-Base
        """
        alignPos = self.seqPos2alignPos(seqPos, querySeqID)
        targetSeqPos = self.alignPos2seqPos(alignPos, targetSeqID)
        return targetSeqPos
    
    def cleanFasta(self, seqID):
        """
        Return sequence without gap
        """
        assert seqID in self.seq_keys
        return self.fasta_dict[seqID].replace(self.gap_sym, "")
    
    def cleanFastaDict(self):
        return { seqID:self.cleanFasta(seqID) for seqID in self.seq_keys}

