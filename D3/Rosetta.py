
"""

This module can predict RNA 3D-structure with Rosetta. Python 3.x is required.

from D3 import Rosetta

seq = "CGGGUAGAGAGGGCAGUGGGAGGUAAGAGCUCUUCACCCUUCACCACCUUCUCCACCCAGC"
dot = ".((((.(((((((..(((((.(((..((.....))))).)))))..))))))).))))..."
seq_title = "hsa-mir-197"
outFolder = "/Share/home/zhangqf7/test_pan"
helix_setup_job,rna_denovo_job,extract_job = Rosetta.pred_3D_rosetta(seq, dot, seq_title, outFolder, gen_modellimit=10, topnum=5, verbose=True, queue="Z-ZQF")

helix_setup_job.wait()
rna_denovo_job.wait()
extract_job.wait()

"""

import os, Cluster, sys

if sys.version_info < (3,0):
    raise RuntimeError("Sorry, requires Python 3.x, not Python 2.x")

os.environ['ROSETTA'] = '/Share/home/zhangqf/lipan/usr/rosetta/rosetta_bin_linux_2019.12.60667_bundle'
os.environ['ROSETTA3'] = os.environ['ROSETTA']
os.environ['ROSETTA3_SRC'] = os.environ['ROSETTA']+"/main/source"
os.environ['ROSETTA3_MAIN'] = os.environ['ROSETTA']+"/main"
os.environ['ROSETTA3_DB'] = os.environ['ROSETTA']+"/main/database"

os.environ['RNA_TOOLS'] = os.environ['ROSETTA']+"/tools/rna_tools"
os.environ['PATH'] = os.environ['RNA_TOOLS']+'/bin/:'+os.environ['ROSETTA']+"/main/source/bin:"+os.environ['PATH']
os.environ['PYTHONPATH'] = os.environ['RNA_TOOLS']+'/bin/:'+os.environ['PYTHONPATH']

helix_setup_CMD = "cd %s;helix_setup.py -fasta %s -secstruct %s -out_prefix %s -out_cmdlines %s"
rna_denovo_CMD = "cd %s;rna_denovo -fasta %s -secstruct_file %s -tag %s -nstruct %s -fixed_stems -s %s -include_neighbor_base_stacks -working_res 1-%s -minimize_rna true"
extract_CMD = "cd %s;extract_lowscore_decoys.py %s %s"

def pred_3D_rosetta(seq, dot, seq_title, outFolder, gen_modellimit=1000, topnum=50, verbose=False, cpu=1, queue="Z-ZQF"):
    """
    seq                 -- RNA sequence
    dot                 -- Dotbracket structure
    seq_title           -- Sequence title: the final output files will be seq_title.out.1.pdb, seq_title.out.2.pdb, ...
    gen_modellimit      -- How many structures to produce
    topnum              -- Get top models from all models
    verbose             -- Show commands
    """
    outFolder = outFolder.rstrip('/') + '/'
    
    seq = seq.lower().replace('t','u')
    assert len(seq) == len(dot)
    assert gen_modellimit > topnum
    
    faFn = outFolder+"seq.fa"
    strFn = outFolder+"dot.str"
    cmdFn = outFolder+"commands"
    errFn = outFolder+"error"
    logFn = outFolder+"log"
    helix_prefix = outFolder+"helix"
    denovoFn = outFolder+seq_title
    
    open(faFn, "w").writelines( ">"+seq_title+"\n"+seq+"\n" )
    open(strFn, "w").writelines( dot+"\n" )
    CMD1 = helix_setup_CMD % (outFolder, faFn, strFn, helix_prefix, cmdFn)
    helix_setup_job = Cluster.new_job(command=CMD1, queue=queue, cpu=cpu, job_name=seq_title+"_helix_setup_1", logFn=logFn,  errFn=errFn)
    
    CMD2 = rna_denovo_CMD % (outFolder, faFn, strFn, denovoFn, gen_modellimit, helix_prefix+"*.pdb", len(seq))
    rna_denovo_job = Cluster.new_job(command=CMD2, queue=queue, cpu=cpu, job_name=seq_title+"_rna_denovo_2", logFn=logFn,  errFn=errFn)
    rna_denovo_job.set_job_depends([helix_setup_job])
    
    CMD3 = extract_CMD % (outFolder, denovoFn+".out", topnum)
    extract_job = Cluster.new_job(command=CMD3, queue=queue, cpu=cpu, job_name=seq_title+"_extract_3", logFn=logFn,  errFn=errFn)
    extract_job.set_job_depends([rna_denovo_job])
    
    helix_setup_job.submit(verbose)
    rna_denovo_job.submit(verbose)
    extract_job.submit(verbose)
    
    return helix_setup_job, rna_denovo_job, extract_job



