#!/bin/env python
import subprocess
import argparse

parser = argparse.ArgumentParser(
    description='Setup scripts to preassemble helices for FARFAR.')
parser.add_argument(
    '-secstruct', required=True, help='Secondary structure file.')
parser.add_argument('-fasta', required=True, help='Fasta file.')
parser.add_argument(
    '-offset', type=int, default=0, help='Residue number offset')
parser.add_argument(
    '-out_prefix', default='helix', help='Prefix of the output scripts.')
parser.add_argument(
    '-out_cmdlines', default='HELIX_CMDLINES',
    help='Output file for the cmdlines for running the jobs.')
args = parser.parse_args()

def read_fasta(fasta_file):
    with open(fasta_file) as f:
        fasta = f.readlines()[1].strip()
    return fasta

# Get the helix stems
def get_helix_stems(secstruct_file):
    def read_sectruct(secstruct_file):
        with open(secstruct_file) as f:
            secstruct = f.readline().strip()
        return secstruct   
    def get_bp_list(secstruct):
        left_base_in_bp = []
        bp = []
        for i, char in enumerate(secstruct):
            if char == '(':
                left_base_in_bp.append(i + 1)
            elif char == ')':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        for i, char in enumerate(secstruct):
            if char == '[':
                left_base_in_bp.append(i + 1)
            elif char == ']':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        return bp

    def get_stems(bp):
        stems = []
        stem = []
        for i, j in bp:
            if (not stem) or (i - 1 == stem[-1][0] and j + 1 == stem[-1][1]):
                stem.append((i, j))
            else:
                stems.append(stem)
                stem = [(i, j)]
        if stem:
            stems.append(stem)
        return stems

    secstruct = read_sectruct(secstruct_file)
    bp = sorted(get_bp_list(secstruct))
    return get_stems(bp)

fasta = read_fasta(args.fasta)
helix_stems = get_helix_stems(args.secstruct)
offset = args.offset

out = open(args.out_cmdlines, 'w')
out.writelines("# helix_setup cmdlines:\n")
stem_idx = 0
min_helix_len = 2
for stem in helix_stems:
    if len(stem) < min_helix_len:
        continue
    seq1 = ''
    seq2 = ''
    for (i, j) in stem:
        seq1 += fasta[i-1]
        seq2 = fasta[j-1] + seq2
    res_num = (
        '%d-%d ' % (stem[0][0] + offset, stem[-1][0] + offset) +
        '%d-%d' % (stem[-1][1] + offset, stem[0][1] + offset)
    )

    tag = "%s%d" % (args.out_prefix, stem_idx)
    out_script = "%s.RUN" % tag
    out_pdb = '%s.pdb ' % tag
    stem_idx += 1

    #cmdline = 'python /Share/home/zhangqf7/yanqiu/rna_denovo/rna_tools/bin/rna_helix.py '
    cmdline = 'rna_helix.py '
    cmdline += '-o %s ' % out_pdb
    cmdline += "-resnum %s " % res_num
    #cmdline += "-out_script %s " % out_script
    #cmdline += "-tag %s " % tag
    cmdline += '-seq %s %s ' % (seq1, seq2)
    out.writelines(cmdline+"\n")
    subprocess.check_call([str(elem) for elem in cmdline.split()])
