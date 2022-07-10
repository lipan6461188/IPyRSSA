#-*- coding:utf-8 -*-

import os, General, random
import tempfile
from xml.dom import minidom
import shutil, Colors

####################################
#### Blast version 1
####################################

class BlastNHit:
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

####################################
#### Blast version 2
####################################

class BlastNHit_V2:
    def __init__(self,bit_score,score,evalue,query_from,query_to,hit_acc,hit_def,hit_from,hit_to,hit_strand, \
                    identity,gap,align_len,match_qseq,match_tseq,align_pat):
        self.bit_score = bit_score
        self.score = score
        self.evalue = evalue
        self.query_from = query_from
        self.query_to = query_to
        self.hit_acc = hit_acc
        self.hit_def = hit_def
        self.hit_from = hit_from
        self.hit_to = hit_to
        self.hit_strand = hit_strand
        self.identity = identity
        self.gap = gap
        self.align_len = align_len
        self.match_qseq = match_qseq
        self.match_tseq = match_tseq
        self.align_pat = align_pat
    
    def __str__(self):
        str1 = "%s\t%d-%d" % (self.match_qseq, self.query_from, self.query_to)
        str2 = "%s\t%s:%d-%d" % (self.match_tseq, self.hit_acc, self.hit_from, self.hit_to)
        return "\n"+str1+"\n"+self.align_pat+"\n"+str2+"\n"
    
    def __repr__(self):
        return self.__str__()

def getChildValue(curElem, childNodeName):
    return curElem.getElementsByTagName(childNodeName)[0].firstChild.nodeValue

def read_blast_result_xml(blastn_xml):
    xmldoc = minidom.parse(blastn_xml)
    iteration = xmldoc.getElementsByTagName('Iteration')
    if len(iteration)>0:
        query_match = iteration[0]
        
        query_id = getChildValue(query_match, 'Iteration_query-def')
        query_len = getChildValue(query_match, 'Iteration_query-len')
        query_len = int(query_len)
        
        hits = query_match.getElementsByTagName("Hit")
        Hit_List = []
        for hit in hits:
            hit_acc = getChildValue(hit, 'Hit_id').split()[0]
            try:
                hit_def = getChildValue(hit, 'Hit_def')
            except:
                hit_def = 'None'
            Hsps = hit.getElementsByTagName("Hsp")
            for hsp in Hsps:
                bit_score = float(getChildValue(hsp, 'Hsp_bit-score'))
                score = float(getChildValue(hsp, 'Hsp_score'))
                evalue = float(getChildValue(hsp, 'Hsp_evalue'))
                query_from = int(getChildValue(hsp, 'Hsp_query-from'))
                query_to = int(getChildValue(hsp, 'Hsp_query-to'))
                hit_from = int(getChildValue(hsp, 'Hsp_hit-from'))
                hit_to = int(getChildValue(hsp, 'Hsp_hit-to'))
                if hit_from>hit_to:
                    hit_strand = '-'
                    hit_from,hit_to = hit_to,hit_from
                else:
                    hit_strand = '+'
                identity = float(getChildValue(hsp, 'Hsp_identity'))
                gap = int(getChildValue(hsp, 'Hsp_gaps'))
                align_len = int(getChildValue(hsp, 'Hsp_align-len'))
                match_qseq = getChildValue(hsp, 'Hsp_qseq')
                match_tseq = getChildValue(hsp, 'Hsp_hseq')
                align_pat = getChildValue(hsp, 'Hsp_midline')
                new_hit = BlastNHit_V2(bit_score,score,evalue,query_from,query_to,hit_acc,hit_def,hit_from,hit_to,hit_strand,identity,gap,align_len,match_qseq,match_tseq,align_pat)
                Hit_List.append(new_hit)
    
    return Hit_List

def blast_sequence_V2(query_seq, blastdb, seqtype='nucl',
    perc_identity=90, strand="both", # Only for blastn
    evalue=10, maxhit=500,  task=None, # For blastn and blastp
    word_size=6, gapopen=10, gapextend=2, clear=True, verbose=False, 
    threads=1, use_LSF=False, LSF_parameters={}, use_SLURM=False, SLURM_parameters={}):
    """
    Blastn a sequence or sequences and return a list of BlastNHit_V2
    
    query_seq               -- Query sequence
    blastdb                 -- Balstn db, such as /150T/zhangqf/GenomeAnnotation/INDEX/blast/hg38/hg38
    seqtype                 -- nucl for DNA, protein for Protein
    perc_identity           -- Minimum percentage of identity (Only for blastn)
    strand                  -- both, minus or plus
    evalue                  -- Maximun evalue
    maxhit                  -- Maximum number of aligned sequences to keep
    task                    -- blastn, blastn-short, dc-megablast, megablast or rmblastn for seqtype='nucl'
                               blastp, blastp-fast, blastp-short for seqtype='protein'
                               default: megablast for seqtype='nucl' and blastp for seqtype='protein'
    clear                   -- Clear tempropary files
    verbose                 -- Output command
    threads                 -- How many threads to use
    use_LSF                 -- Submit to LSF if True
    LSF_parameters          -- { 'queue': 'Z-ZQF', 'cpu': 20, 'job_name': 'cmsearch', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    use_SLURM               -- Submit to SLURM if True
    SLURM_parameters        -- { 'partition': 'normal02', 'job_name': 'blast', 'logFn': '/dev/null', 'errFn': '/dev/null' }
    
    Return: [ BlastNHit_V2, BlastNHit_V2,... ]
    
    Require: blastn or blastp
    """
    
    if seqtype=='nucl':
        blastn = General.require_exec("blastn", "blastn is required")
        if not os.path.exists(blastdb+".nhr"):
            raise RuntimeError("Error: blast index is not valid")
    elif seqtype=='protein':
        blastp = General.require_exec("blastp", "blastn is required")
        if not os.path.exists(blastdb+".pdb"):
            raise RuntimeError("Error: blast index is not valid")
    else:
        raise RuntimeError("seqtype should be nucl or protein")
    
    if type(query_seq) is not str:
        raise RuntimeError("Error: parameter query_seq should a sequence")
    
    if use_LSF or use_SLURM:
        tmpdir = tempfile.mkdtemp(prefix="blast_", dir=os.path.join(os.environ['HOME'],'tmp'))
    else:
        tmpdir = tempfile.mkdtemp(prefix="blast_")
    tmp_query_fa = os.path.join(tmpdir, "query.fa")
    tmp_balst_xml = os.path.join(tmpdir, "result.xml")
    
    with open(tmp_query_fa, 'w') as IN:
        print(">query-seq", file=IN)
        print(query_seq, file=IN)
    
    if task is None:
        task = 'megablast' if seqtype=='nucl' else 'blastp'
    
    if seqtype=='nucl':
        cmd = f"{blastn} -word_size {word_size} -gapopen {gapopen} -gapextend {gapextend} -query {tmp_query_fa} -strand {strand} -task {task} -db {blastdb} -out {tmp_balst_xml} -outfmt 5 -max_hsps 15 -num_threads {threads} -perc_identity {perc_identity} -evalue {evalue} -max_target_seqs {maxhit}"
    if seqtype=='protein':
        cmd = f"{blastp} -word_size {word_size} -gapopen {gapopen} -gapextend {gapextend} -query {tmp_query_fa} -task {task} -db {blastdb} -out {tmp_balst_xml} -outfmt 5 -max_hsps 15 -num_threads {threads} -evalue {evalue} -max_target_seqs {maxhit}"
    
    if use_LSF:
        import Cluster
        job = Cluster.new_job(command=cmd, 
            queue=LSF_parameters.get('queue', 'Z-ZQF'), 
            cpu=LSF_parameters.get('cpu', 20), 
            job_name=LSF_parameters.get('job_name', 'Alignment.blast_sequence_V2'), 
            logFn=LSF_parameters.get('logFn', '/dev/null'),
            errFn=LSF_parameters.get('errFn', '/dev/null'))
        if verbose:
            print(Colors.f(job.get_submit_command(), fc='yellow'))
        job.submit()
        job.wait()
    elif use_SLURM:
        partition = SLURM_parameters.get('partition', 'normal02')
        job_name = SLURM_parameters.get('job_name', 'Alignment.blast_sequence_V2')
        logFn = SLURM_parameters.get('logFn', None)
        errFn = SLURM_parameters.get('errFn', None)
        slurm_cmd = f"srun -p {partition} -J {job_name} -c {threads} "
        if logFn:
            slurm_cmd += "-o {logFn} "
        if errFn:
            slurm_cmd += "-e {logFn} "
        slurm_cmd += cmd
        if verbose:
            print(Colors.f(slurm_cmd, fc='yellow'))
        _ = os.system(slurm_cmd)
        #from multiprocessing import Process
        #job = Process(slurm_cmd)
        #job.start()
        #job.join()
    else:
        if verbose:
            print(Colors.f(cmd, fc='yellow'))
        if clear:
            status = General.Run_catchKI(cmd, [tmpdir])
        else:
            status = os.system(cmd)
        if status != 0:
            raise RuntimeError("Error: Blast has an error")
    
    if not os.path.exists(tmp_balst_xml):
        raise RuntimeError(f"Error: Blast has an error, {tmp_balst_xml} not found")
    hits = read_blast_result_xml(tmp_balst_xml)
    if clear:
        shutil.rmtree(tmpdir)
    
    return hits


####################################
#### Sequence deduplication
####################################

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

def cd_hit_est(id2seq, identity=0.9, global_align=False, band_width=20, memory=800, cpu=0, word_length=10,
    cluster_mode=0, alignment_mode=1, clean=True, verbose=False):
    """
    id2seq                  -- {id1:seq1, ...}
    identity                -- sequence identity threshold. 
                               this is the default cd-hit's "global sequence identity" calculated as:
                               number of identical amino acids or bases in alignment
                               divided by the full length of the shorter sequence
    global_align            -- use global sequence identity
    band_width              -- band_width of alignment
    memory                  -- memory limit (in MB) for the program, 0 for unlimited
    cpu                     -- number of threads, with 0, all CPUs will be used
    word_length             -- word_length, default 10, see user's guide for choosing it
    cluster_mode            -- 1 or 0.
                               by cd-hit's default algorithm, a sequence is clustered to the first
                               cluster that meet the threshold (fast cluster). If set to 1, the program
                               will cluster it into the most similar cluster that meet the threshold
                               (accurate but slow mode)
    alignment_mode          -- 1 or 0, default 1, by default do both +/+ & +/- alignments
                               if set to 0, only +/+ strand alignment
    clean                   -- Delete all tmp files
    verbose                 -- Print command
    
    Remove duplicated sequence from sequence
    
    Return:
        {id1:seq1,...}
    
    Require: cd-hit-est
    """
    import General, Colors
    import shutil, tempfile
    
    exe = General.require_exec("cd-hit-est", exception=True)
    
    ROOT = tempfile.mkdtemp(prefix="cd_hit_est_")
    
    #ROOT = "/tmp/predict_structure_%s/" % (randID, )
    #os.mkdir(ROOT)
    fa_file =  os.path.join(ROOT, "input.fa")
    output_file = os.path.join(ROOT, "output.fa")
    
    General.write_fasta(id2seq, fa_file)
    CMD = f"{exe} -i {fa_file} -o {output_file} -c {identity} -b {band_width} -M {memory} -T {cpu} -n {word_length} -g {cluster_mode} -r {alignment_mode} "
    if global_align:
        CMD += " -G "
    
    if verbose:
        print(Colors.f(CMD, fc='yellow'))
    
    CMD += ' > /dev/null'
    os.system(CMD)
    
    cleaned_fasta = General.load_fasta(output_file)
    if clean:
        shutil.rmtree(ROOT)
    
    return cleaned_fasta


