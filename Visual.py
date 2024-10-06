#-*- coding:utf-8 -*-

"""

This module use VARNA (http://varna.lri.fr) to visualize the RNA secondary structure
Set the PATH of VARNAv3-93.jar in variable VARNA
export VARNA=[PATH to VARNAv3-93.jar]

"""

import os, sys, random, platform, shutil, tempfile, subprocess
from . import Colors

if 'VARNA' in os.environ:
    VARNAProg = os.environ['VARNA']
else:
    from . import Colors
    sys.stderr.writelines( Colors.f("Warning: VARNA variable not found, please specify",fc="red",bc="white",ft="blink") + "\n" )
    VARNAProg = "/Users/lee/Documents/VARNAv3-93.jar"

Gradient_Colors = Colors.linear_gradient("#0b6b38", "#fffbb8", 50)['hex'] + Colors.linear_gradient("#fffbb8", "#aa0626", 50)['hex']

def __highlight_region_cmd(region_list):
    """
    The input must be [ (1, 10, '#FF0000'), (90, 100, '#00FF00') ]...
    1-based coordination system
    """
    assert( (isinstance(region_list, list) or isinstance(region_list, tuple)) and len(region_list) >= 1 )
    assert( isinstance(region_list[0], list) or isinstance(region_list[0], tuple) )
    CMD = '-highlightRegion "'
    for region in region_list:
        if len(region) == 2:
            CMD += "%s-%s:fill=%s,outline=#FFFFFF,radius=15;" % (region[0], region[1], "00FF00")
        else:
            CMD += "%s-%s:fill=%s,outline=#FFFFFF,radius=15;" % (region[0], region[1], region[2])
    CMD += '"'
    return CMD

def __base_color_shape_cmd(shape_list, cutofflist=[0.3,0.5,0.7]):
    
    assert len(cutofflist)==3
    
    Level1 = ""
    Level2 = ""
    Level3 = ""
    Level4 = ""
    NoData = ""
    for idx in range(len(shape_list)):
        if shape_list[idx] == 'NULL':
            NoData += str(idx+1) if NoData == "" else ','+str(idx+1)
        elif float(shape_list[idx]) > cutofflist[2]:
            Level1 += str(idx+1) if Level1 == "" else ','+str(idx+1)
        elif float(shape_list[idx]) > cutofflist[1]:
            Level2 += str(idx+1) if Level2 == "" else ','+str(idx+1)
        elif float(shape_list[idx]) > cutofflist[0]:
            Level3 += str(idx+1) if Level3 == "" else ','+str(idx+1)
        else:
            Level4 += str(idx+1) if Level4 == "" else ','+str(idx+1)
    CMD = ""
    if Level1: CMD += "-applyBasesStyle1on \"%s\" " % (Level1, )
    if Level2: CMD += "-applyBasesStyle2on \"%s\" " % (Level2, )
    if Level3: CMD += "-applyBasesStyle3on \"%s\" " % (Level3, )
    if Level4: CMD += "-applyBasesStyle4on \"%s\" " % (Level4, )
    if NoData: CMD += "-applyBasesStyle5on \"%s\" " % (NoData, )
    return CMD

def __dot_match_bpprob(dot, bpprob, warning=True):
    """
    dot                 -- Dot-bracket
    bpprob              -- [(base1, base2, prob), ...]
    """
    import Structure
    bpprob.sort(key=lambda x: x[0])
    dot_bp = sorted(Structure.dot2ct(dot), key=lambda x: x[0])
    
    i, j = 0, 0
    final_bpprob = []
    while i<len(bpprob) and j<len(dot_bp):
        if bpprob[i][:2]==dot_bp[j][:2]:
            final_bpprob.append( bpprob[i] )
            i += 1
            j += 1
        elif bpprob[i][0]>dot_bp[j][0]:
            final_bpprob.append( list(dot_bp[j])+[None] ) # default
            j += 1
        else:
            x,y = bpprob[i][:2]
            if x<=y and warning:
                print(f"Warning: ({x},{y}) is not a base pair in structure. Skip it.")
            i += 1
    while j<len(dot_bp):
        final_bpprob.append( list(dot_bp[j])+[None] ) # default
        j += 1
    return final_bpprob

def __basepair_bpprob_cmd(bpprob, cutofflist=[0.6,0.8,0.95], mode='color'):
    """
    bpprob             -- [(base1, base2, prob), ...]
    cutofflist          -- The cutoff of values
    mode                -- color/thickness/both
    """
    if len(cutofflist)!=3:
        raise RuntimeError("Error: cutofflist should have length of 3")
    
    Level1 = []
    Level2 = []
    Level3 = []
    Level4 = []
    NoData = []
    for b1,b2,prob in bpprob:
        if prob is None:
            NoData.append( (b1,b2,prob) )
        elif prob > cutofflist[2]:
            Level1.append( (b1,b2,prob) )
        elif prob > cutofflist[1]:
            Level2.append( (b1,b2,prob) )
        elif prob > cutofflist[0]:
            Level3.append( (b1,b2,prob) )
        else:
            Level4.append( (b1,b2,prob) )
    CMD = ""
    if Level1:
        for b1,b2,prob in Level1:
            if mode=='color':
                CMD += f"({b1},{b2}):color=#2306f7;"
            elif mode=='both':
                CMD += f"({b1},{b2}):color=#2306f7,Thickness=4;"
            elif mode=='thickness':
                CMD += f"({b1},{b2}):Thickness=4;"
    if Level2:
        for b1,b2,prob in Level2:
            if mode=='color':
                CMD += f"({b1},{b2}):color=#7167f9;"
            elif mode=='both':
                CMD += f"({b1},{b2}):color=#7167f9,Thickness=3;"
            elif mode=='thickness':
                CMD += f"({b1},{b2}):Thickness=3;"
    if Level3:
        for b1,b2,prob in Level3:
            if mode=='color':
                CMD += f"({b1},{b2}):color=#b7b2f7;"
            elif mode=='both':
                CMD += f"({b1},{b2}):color=#b7b2f7,Thickness=2;"
            elif mode=='thickness':
                CMD += f"({b1},{b2}):Thickness=2;"
    if Level4:
        for b1,b2,prob in Level4:
            if mode=='color':
                CMD += f"({b1},{b2}):color=#e4e3fc;"
            elif mode=='both':
                CMD += f"({b1},{b2}):color=#e4e3fc,Thickness=1;"
            elif mode=='thickness':
                CMD += f"({b1},{b2}):Thickness=1;"
    if NoData:
        for b1,b2,prob in NoData:
            if mode=='color':
                CMD += f"({b1},{b2}):color=#aeaeaf;"
            elif mode=='both':
                CMD += f"({b1},{b2}):color=#aeaeaf,Thickness=1;"
            elif mode=='thickness':
                CMD += f"({b1},{b2}):Thickness=1;"
    
    CMD = f"-auxBPs \"{CMD}\""
    return CMD

def __manual_period(seqlen, first_base_pos, period, peroid_color='#ff5722'):
    """
    seqlen                      -- The sequence length
    first_base_pos              -- The position of the first base
    period                      -- Numbering period
    peroid_color                -- Number color
    """
    cmd = "-periodNum 0 "
    j = 0
    annotation_cmd = ""
    for i in range(first_base_pos, first_base_pos+seqlen):
        j += 1
        if i%period == 0:
            annotation_cmd += f"{i}:type=B,anchor={j},size=8,color={peroid_color};"
    if annotation_cmd:
        cmd += f"-annotations \"{annotation_cmd}\""
    return cmd

def __base_color_heatmap_cmd(shape_list):
    shape_str = ""
    null_base_idx = []
    for idx,shape in enumerate(shape_list):
        if shape == 'NULL':
            shape_str += '0.0;'
            null_base_idx.append(str(idx+1))
        else:
            shape_str += '%s;' % (shape, )
    CMD = ""
    if null_base_idx:
        CMD += "-basesStyle1 \"label=#828282\" -applyBasesStyle1on \"%s\" " % ( ",".join(null_base_idx), )
    CMD += "-colorMap \"%s\" " % (shape_str[:-1], )
    return CMD

def __annotation_cmd(annotation_list):
    """
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    annotation_str = ""
    for annot in annotation_list:
        assert 'text' in annot
        assert 'anchor' in annot
        annotation_str += ';%s:anchor=%d,size=%d,color=%s,type=%s' % (annot['text'], annot['anchor'], annot.get('size',7), annot.get('color','#000000'), annot.get('type','B'))
    return "-annotations \"" + annotation_str[1:] + "\""

def Plot_RNAStructure_Shape(sequence, dot, shape_list, 
    mode='label', correctT=True, scaling=0.8, 
    highlight_region=[], annotation=[], cutofflist=[0.3,0.5,0.7], 
    bpprob=[], bpprob_cutofflist=[0.6,0.8,0.95], bpprob_mode='color', bpwarning=True,
    period=10, first_base_pos=1, peroid_color='#000000',
    title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    cutofflist          -- The color cutoff
    bpprob              -- Base pairing probability, only provide base pairs in the structure
    bpprob_cutofflist   -- Base pairing color/thickness cutoff
    bpprob_mode         -- color/thickness/both
    bpwarning           -- Base pairing warning when provide base pairs not in the structure
    period              -- Numbering the base for how many bases as period
    first_base_pos      -- The number of first base
    peroid_color        -- The period color
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot) == len(shape_list)
    assert mode in ('label', 'fill', 'heatmap')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -drawBases false -bpStyle simple " % (sequence, dot)
    
    if mode == 'label':
        CMD += "-basesStyle1 \"label=#B61D22\" " + "-basesStyle2 \"label=#ED9616\" " + "-basesStyle3 \"label=#194399\" " + "-basesStyle4 \"label=#040000\" " + "-basesStyle5 \"label=#828282\" "
        CMD += __base_color_shape_cmd(shape_list, cutofflist=cutofflist)
    elif mode == 'fill':
        CMD += "-basesStyle1 \"fill=#B61D22\" " + "-basesStyle2 \"fill=#ED9616\" " + "-basesStyle3 \"fill=#194399\" " + "-basesStyle4 \"fill=#040000\" " + "-basesStyle5 \"fill=#828282\" "
        CMD += __base_color_shape_cmd(shape_list, cutofflist=cutofflist)
    else:
        CMD += __base_color_heatmap_cmd(shape_list)
    
    if highlight_region:
        CMD += " " + __highlight_region_cmd(highlight_region)
    
    if annotation:
        CMD += " " + __annotation_cmd(annotation)
    
    if scaling:
        CMD += " -spaceBetweenBases \"%s\"" % (scaling, )
    
    if bpprob:
        new_bpprob = __dot_match_bpprob(dot, bpprob, bpwarning)
        CMD += " " + __basepair_bpprob_cmd(new_bpprob, bpprob_cutofflist, bpprob_mode)
    
    if first_base_pos==1 and peroid_color=='#000000':
        CMD += f" -period {period}"
    else:
        CMD += " " + __manual_period(len(sequence), first_base_pos, period, peroid_color)
    
    if title:
        CMD += " -title \"%s\"" % (title, )
    
    if not wait:
        CMD = "nohup " + CMD + " &"
    
    return CMD

def __base_color_base_cmd(seq):
    A = ""
    T = ""
    C = ""
    G = ""
    for idx in range(len(seq)):
        if seq[idx] == 'A':
            A += str(idx+1) if A == "" else ','+str(idx+1)
        if seq[idx] == 'C':
            C += str(idx+1) if C == "" else ','+str(idx+1)
        if seq[idx] == 'T' or seq[idx] == 'U':
            T += str(idx+1) if T == "" else ','+str(idx+1)
        if seq[idx] == 'G':
            G += str(idx+1) if G == "" else ','+str(idx+1)
    CMD = ""
    if A: CMD += "-applyBasesStyle1on \"%s\" " % (A, )
    if T: CMD += "-applyBasesStyle2on \"%s\" " % (T, )
    if C: CMD += "-applyBasesStyle3on \"%s\" " % (C, )
    if G: CMD += "-applyBasesStyle4on \"%s\" " % (G, )
    return CMD

def Plot_RNAStructure_Base(sequence, dot, mode='fill', correctT=True, scaling=0.8, 
    highlight_region=[], annotation=[], 
    bpprob=[], bpprob_cutofflist=[0.6,0.8,0.95], bpprob_mode='color', bpwarning=True,
    period=10, first_base_pos=1, peroid_color='#000000',
    title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    bpprob              -- Base pairing probability, only provide base pairs in the structure
    bpprob_cutofflist   -- Base pairing color/thickness cutoff
    bpprob_mode         -- color/thickness/both
    bpwarning           -- Base pairing warning when provide base pairs not in the structure
    period              -- Numbering the base for how many bases as period
    first_base_pos      -- The number of first base
    peroid_color        -- The period color
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot)
    assert mode in ('label', 'fill')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -drawBases false -bpStyle simple " % (sequence, dot)
    
    if mode == 'label':
        CMD += "-basesStyle1 \"label=#CCFFCC\" " + "-basesStyle2 \"label=#CCFFFF\" " + "-basesStyle3 \"label=#FFFFCC\" " + "-basesStyle4 \"label=#FFCCFF\" "
    else:
        CMD += "-basesStyle1 \"fill=#CCFFCC\" " + "-basesStyle2 \"fill=#CCFFFF\" " + "-basesStyle3 \"fill=#FFFFCC\" " + "-basesStyle4 \"fill=#FFCCFF\" "
    
    CMD += __base_color_base_cmd(sequence)
    
    if highlight_region:
        CMD += " " + __highlight_region_cmd(highlight_region)
    
    if annotation:
        CMD += " " + __annotation_cmd(annotation)

    if scaling:
        CMD += " -spaceBetweenBases \"%s\"" % (scaling, )
    
    if bpprob:
        new_bpprob = __dot_match_bpprob(dot, bpprob, bpwarning)
        CMD += " " + __basepair_bpprob_cmd(new_bpprob, bpprob_cutofflist, bpprob_mode)
    
    if first_base_pos==1 and peroid_color=='#000000':
        CMD += f" -period {period}"
    else:
        CMD += " " + __manual_period(len(sequence), first_base_pos, period, peroid_color)

    if title:
        CMD += " -title \"%s\"" % (title, )
    
    if not wait:
        CMD = "nohup " + CMD + " &"
    
    return CMD

def Plot_RNAStructure_highlight(sequence, dot, hg_base_list=[], mode='fill', correctT=True, 
    scaling=0.8, highlight_region=[], annotation=[], 
    bpprob=[], bpprob_cutofflist=[0.6,0.8,0.95], bpprob_mode='color', bpwarning=True,
    period=10, first_base_pos=1, peroid_color='#000000',
    title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    hg_base_list        -- Regions of bases to highlight in circle
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    bpprob              -- Base pairing probability, only provide base pairs in the structure
    bpprob_cutofflist   -- Base pairing color/thickness cutoff
    bpprob_mode         -- color/thickness/both
    bpwarning           -- Base pairing warning when provide base pairs not in the structure
    period              -- Numbering the base for how many bases as period
    first_base_pos      -- The number of first base
    peroid_color        -- The period color
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot)
    assert mode in ('label', 'fill')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -drawBases false -bpStyle simple " % (sequence, dot)
    
    if hg_base_list:
        if mode == 'label':
            CMD += "-basesStyle1 \"label=#FF0000\" "
        else:
            CMD += "-basesStyle1 \"fill=#FF0000\" "
        hg_base_list = [str(item) for item in hg_base_list]
        CMD += "-applyBasesStyle1on \"%s\" " % (",".join(hg_base_list), )
    
    if highlight_region:
        CMD += " " + __highlight_region_cmd(highlight_region)
    
    if annotation:
        CMD += " " + __annotation_cmd(annotation)
    
    if scaling:
        CMD += " -spaceBetweenBases \"%s\"" % (scaling, )
    
    if bpprob:
        new_bpprob = __dot_match_bpprob(dot, bpprob, bpwarning)
        CMD += " " + __basepair_bpprob_cmd(new_bpprob, bpprob_cutofflist, bpprob_mode)
    
    if first_base_pos==1 and peroid_color=='#000000':
        CMD += f" -period {period}"
    else:
        CMD += " " + __manual_period(len(sequence), first_base_pos, period, peroid_color)
    
    if title:
        CMD += " -title \"%s\"" % (title, )
    
    if not wait:
        CMD = "nohup " + CMD + " &"
    
    return CMD


##################################
####    rRNA structure ps
##################################

PS_PATH = os.path.dirname(os.path.abspath(__file__))+"/ps/"
PS_FILE = { 
    'yeast_small':      'yeast_small.ps',
    'yeast_large':      ['yeast_large_5.ps', 'yeast_large_3.ps'],
    'yeast_5S':         'yeast_5S.ps',
    'yeast_smallMito':  'yeast_smallMito.ps',

    'human_small':      'human_small.ps',
    'human_5S':         'human_5S.ps',
    'human_smallMito':  'human_smallMito.ps',

    'mouse_small':      'mouse_small.ps',
    'mouse_5S':         'mouse_5S.ps',
    'mouse_smallMito':  'mouse_smallMito.ps',

    'arabidopsis_small':    'arabidopsis_small.ps',
    'arabidopsis_large':    ['arabidopsis_large_5.ps', 'arabidopsis_large_3.ps']
}

def load_ps(psFn):
    import re
    header_lines = []
    tail_lines = []
    base_lines = []
    base_through = False
    seq = ""
    for line in open(psFn):
        if re.match(r"^\([AUTCG]\)", line):
            base_lines.append(line)
            base_through = True
            seq += line[1]
        elif base_through:
            tail_lines.append(line)
        else:
            header_lines.append(line)
    
    return header_lines, base_lines, tail_lines, seq

def Map_rRNA_Shape(sequence, shape_list, target, outPrex, title=None, mode='label', cutofflist=[0.3,0.5,0.7]):
    """
    sequence            -- Raw sequence is used to double-check
    shape_list          -- A list of SHAPE scores
    target              -- rRNA target
                           yeast_small, yeast_large, yeast_5S, yeast_smallMito,
                           human_small, human_5S, human_smallMito,
                           mouse_small, mouse_5S, mouse_smallMito,
                           arabidopsis_small, arabidopsis_large
    outPrex             -- Output file prefix
    title               -- Title for PDF file
    mode                -- Color mode: [label|fill|heatmap]
    cutofflist          -- The color cutoff
    """
    assert target in PS_FILE
    assert mode in ('label', 'fill', 'heatmap')
    assert len(sequence) == len(shape_list)
    
    title_map = {   'yeast_small':'Yeast small subunit rRNA',
                    'yeast_large':'Yeast large subunit rRNA',
                    'yeast_5S': 'Yeast 5S rRNA',
                    'yeast_smallMito': 'Yeast mitochodria small subunit rRNA',
                    'human_small':'Human small subunit rRNA',
                    'human_5S': 'Human 5S rRNA',
                    'human_smallMito': 'Human mitochodria small subunit rRNA',
                    'mouse_small':'Mouse small subunit rRNA',
                    'mouse_5S': 'Mouse 5S rRNA',
                    'mouse_smallMito': 'Mouse mitochodria small subunit rRNA',
                    'arabidopsis_small': 'Arabidopsis small subunit rRNA',
                    'arabidopsis_large': 'Arabidopsis large subunit rRNA' }

    if title is None:
        title = title_map[target]
    
    psfiles = PS_FILE[target]
    if not isinstance(psfiles, list):
        psfiles = [psfiles]
    start=0
    for psfile in psfiles:
        full_path = PS_PATH + psfile
        header_lines, base_lines, tail_lines, seq = load_ps(full_path)
        cur_len = len(seq)
        #print( len(sequence[start:start+cur_len]) )
        #print( len(seq) )
        if sequence[start:start+cur_len].replace('T','U') != seq.replace('T','U'):
            #sys.stderr.writelines(sequence[start:start+cur_len].replace('T','U')+"|\n")
            #sys.stderr.writelines(seq.replace('T','U')+"|\n")
            sys.stderr.writelines("Error: Different Sequence!!\n")
            return
        if '_5.ps' in psfile:
            outFn = outPrex + '_5p.ps'
        elif '_3.ps' in psfile:
            outFn = outPrex + '_3p.ps'
        else:
            outFn = outPrex + '.ps'
        #if mode == 'label':
        print("Write file: "+outFn+"...")
        assert len(base_lines)==cur_len
        #print(base_lines[:100], cur_len)
        write_label_ps(header_lines, base_lines, tail_lines, shape_list[start:start+cur_len], title, outFn, cutofflist=cutofflist, mode=mode)
        #elif mode == 'fill':
        #    sys.stderr.writelines("Sorry: Not applicant now\n")
        #else:
        #    sys.stderr.writelines("Sorry: Not applicant now\n")
        start += cur_len

def get_rRNA_refseq(target):
    """
    Return reference sequence
    target              -- rRNA target
                           yeast_small, yeast_large, yeast_5S, yeast_smallMito
                           human_small, human_5S, human_smallMito,
                           mouse_small, mouse_5S, mouse_smallMito,
                           arabidopsis_small, arabidopsis_large
    """
    assert target in PS_FILE
    
    fullseq = ""
    psfiles = PS_FILE[target]
    if not isinstance(psfiles, list):
        psfiles = [psfiles]
    for psfile in psfiles:
        full_path = PS_PATH + psfile
        header_lines, base_lines, tail_lines, seq = load_ps(full_path)
        fullseq += seq
    return fullseq

def _color_command_segmented(base_score, cutofflist=[0.3,0.5,0.7]):
    if base_score == 'NULL':
        return "0.51 0.51 0.51 setrgbcolor"
    elif float(base_score) < cutofflist[0]:
        return "0.00 0.00 0.00 setrgbcolor"
    elif float(base_score) < cutofflist[1]:
        return "0.10 0.26 0.60 setrgbcolor"
    elif float(base_score) < cutofflist[2]:
        return "0.93 0.59 0.09 setrgbcolor"
    else:
        return "0.71 0.11 0.13 setrgbcolor"

def _color_command_heatmap(base_score, gradient_list, min_score=0, max_score=1):
    import numpy as np
    if base_score == 'NULL':
        return "0.51 0.51 0.51 setrgbcolor"
    else:
        v = np.clip(float(base_score), min_score, max_score)
        ind = int( (v-min_score)/(max_score-min_score) * len(gradient_list) )
        #print(ind)
        ind = int(np.clip(ind, 0, len(gradient_list)-1))
        color = gradient_list[ind]
        #print(color)
        r,g,b = Colors._hex_to_RGB(color)
        r,g,b = r/255,g/255,b/255
        #print(r,g,b)
        return f"{r:.2f} {g:.2f} {b:.2f} setrgbcolor"


def write_label_ps(header_lines, base_lines, tail_lines, shape_list, title, outFn, cutofflist=[0.3,0.5,0.7], mode='fill'):
    """
    header_lines            -- produced by load_ps
    base_lines              -- produced by load_ps
    tail_lines              -- produced by load_ps
    shape_list              -- A list of SHAPE scores
    title                   -- Title
    outFn                   -- Output file name
    cutofflist              -- The color cutoff
    """
    OUT = open(outFn, "w")
    for header_line in header_lines:
        if r'{title}' in header_line:
            header_line = header_line.format(title=title)
        OUT.writelines(header_line)
    #print(len(shape_list), len())
    for shape,base_line in zip(shape_list,base_lines):
        if mode=='label':
            OUT.writelines( _color_command_segmented(shape, cutofflist)+"\n" )
        elif mode=='heatmap':
            OUT.writelines( _color_command_heatmap(shape, Gradient_Colors, 0, 1)+"\n" )
        else:
            raise RuntimeError("Sorry: mode='fill' Not applicant now")
        OUT.writelines(base_line)
    for tail_line in tail_lines:
        OUT.writelines(tail_line)
    OUT.close()


def visual_structure_entropy(sequence, shape_list=None, si=-0.6, sm=1.8, md=None, figsize=(15,6)):
    """
    Plot the shannon entropy and the pairing probability
    
    sequence                -- Sequence
    shape_list              -- Shape list
    """
    import Structure, Figures
    import matplotlib.pyplot as plt
    
    if shape_list:
        assert len(sequence) == len(shape_list)
    Len = len(sequence)
    
    probList = Structure.partition(sequence, shape_list=shape_list, si=si, sm=sm, md=md, return_pfs=False)
    shannon = Structure.calcShannon(probList, Len)
    
    # Plot figure
    fig, axs = plt.subplots(2,1,figsize=figsize)
    axs[0].bar(range(1, Len+1), shannon, width=1)
    axs[0].set_ylim(0, 0.4)
    axs[0].set_xlim(1, Len)
    axs[0].set_ylabel("Shannon Entropy")
    Figures.rainbowPlot(probList, axs[1], length=Len, lw=1)
    return fig, axs



#####################################
### Protein Structure Visualization
#####################################


class Chimera:
    """
    Run UCSF Chimera with python binding
    
    Example:
    
    os.environ['LD_LIBRARY_PATH'] += ':/nfs_beijing/kubeflow-user/lipan/miniconda3/envs/torch2/lib'

    chimera = Chimera()
    chimera_command_list = [
        'background solid white',
        'open AF-P54219-F1-model_v4.pdb',
        'rangecolor bfactor, 50 #f08253 70 #fada4d 90 #7ec9ef 100 #1b57ce',
        'copy file AF-P54219-F1-model_v4.png png width 800 height 800'
    ]
    chimera.run(chimera_command_list, nogui=True)
    """
    def __init__(self, chimera_bin=None, work_dir=None, delete_workdir=None):
        if chimera_bin is None:
            chimera_bin = shutil.which('chimera')
        self.chimera_bin = chimera_bin
        assert os.path.exists(self.chimera_bin), self.chimera_bin
        if work_dir is None:
            work_dir = tempfile.mkdtemp(prefix='chimera_')
            if delete_workdir is None:
                delete_workdir = True
        elif delete_workdir is None:
            delete_workdir = False
        self.work_dir         = work_dir
        self.delete_workdir   = delete_workdir
        self.pre_command_list = ['background solid white']
        self.post_command_list= []
    
    def set_coordinate_system(self, length=100, x_color='red', y_color='green', z_color='blue'):
        text = f".color {x_color}\n.arrow 0 0 0 {length} 0 0\n.color {y_color}\n.arrow 0 0 0 0 {length} 0\n.color {z_color}\n.arrow 0 0 0 0 0 {length}"
        coordinate_file = os.path.join(self.work_dir, 'coordinates.bild')
        print(text, file=open(coordinate_file, 'w'))
        self.post_command_list.append(f"open {coordinate_file}")
    
    def run(self, chimera_command_list, nogui=True, stop_at_end=True):
        chimera_py_script = tempfile.mktemp(suffix='.py', prefix='chimera_', dir=self.work_dir)
        with open(chimera_py_script, 'w') as OUT:
            print("from chimera import runCommand as rc", file=OUT)
            for cmd in self.pre_command_list:
                print(f"rc('{cmd}')", file=OUT)
            for cmd in chimera_command_list:
                print(f"rc('{cmd}')", file=OUT)
            for cmd in self.post_command_list:
                print(f"rc('{cmd}')", file=OUT)
            if stop_at_end:
                print("rc('stop now')", file=OUT)
        cmd = f"{self.chimera_bin} --silent --script {chimera_py_script}"
        if nogui:
            cmd += ' --nogui'
        status, output = subprocess.getstatusoutput(cmd)
        if status != 0 or len(output) > 0:
            print("===============UCSF Chimera Error information===============")
            print(output)
    
    def clear(self):
        import os
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
    
    def __del__(self):
        if self.delete_workdir:
            self.clear()

