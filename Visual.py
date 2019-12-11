#-*- coding:utf-8 -*-

"""

This module use VARNA (http://varna.lri.fr) to visualize the RNA secondary structure
Set the PATH of VARNAv3-93.jar in variable VARNA
export VARNA=[PATH to VARNAv3-93.jar]

"""

import os, sys, random, platform

if 'VARNA' in os.environ:
    VARNAProg = os.environ['VARNA']
else:
    import Colors
    sys.stderr.writelines( Colors.f("Warning: VARNA variable not found, please specify",fc="red",bc="white",ft="blink") + "\n" )
    VARNAProg = "/Users/lee/Documents/VARNAv3-93.jar"

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

def Plot_RNAStructure_Shape(sequence, dot, shape_list, mode='label', correctT=True, scaling=0.8, highlight_region=[], annotation=[], cutofflist=[0.3,0.5,0.7], title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    shape_list          -- A list of SHAPE scores
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    cutofflist          -- The color cutoff
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot) == len(shape_list)
    assert mode in ('label', 'fill', 'heatmap')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -bpStyle simple " % (sequence, dot)
    
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

def Plot_RNAStructure_Base(sequence, dot, mode='fill', correctT=True, scaling=0.8, highlight_region=[], annotation=[], title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot)
    assert mode in ('label', 'fill')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -bpStyle simple " % (sequence, dot)
    
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
    
    if title:
        CMD += " -title \"%s\"" % (title, )
    
    if not wait:
        CMD = "nohup " + CMD + " &"
    
    return CMD

def Plot_RNAStructure_highlight(sequence, dot, hg_base_list=[], mode='fill', correctT=True, scaling=0.8, highlight_region=[], annotation=[], title="", wait=True, VARNAProg=VARNAProg):
    """
    sequence            -- Raw sequence
    dot                 -- Dotbracket structure
    hg_base_list        -- Regions of bases to highlight in circle
    mode                -- Color mode: [label|fill|heatmap]
    correctT            -- Covert T to U
    highlight_region    -- Regions to highlight
    title               -- Title of plot
    wait                -- Nohup the command
    VARNAProg           -- Path of VARNA
    
    annotation_list    -- [ {'text': 'loop1', 'anchor':10, 'color': '#ff9800', 'size': 10, 'type':'B'},... ]
    """
    
    assert len(sequence) == len(dot)
    assert mode in ('label', 'fill')
    
    if correctT:
        sequence = sequence.replace('T', 'U')
    
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -drawBackbone false -bpStyle simple " % (sequence, dot)
    
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
    'mouse_smallMito':  'mouse_smallMito.ps'
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

def Map_rRNA_Shape(sequence, shape_list, target, outPrex, mode='label', cutofflist=[0.3,0.5,0.7]):
    """
    sequence            -- Raw sequence is used to double-check
    shape_list          -- A list of SHAPE scores
    target              -- rRNA target
                           yeast_small, yeast_large, yeast_5S, yeast_smallMito,
                           human_small, human_5S, human_smallMito,
                           mouse_small, mouse_5S, mouse_smallMito
    outPrex             -- Output file prefix
    mode                -- Color mode: [label|fill|heatmap]
    cutofflist          -- The color cutoff
    """
    assert target in PS_FILE
    assert mode in ('label', 'fill', 'heatmap')
    assert len(sequence) == len(shape_list)
    
    psfiles = PS_FILE[target]
    if not isinstance(psfiles, list):
        psfiles = [psfiles]
    start=0
    for psfile in psfiles:
        full_path = PS_PATH + psfile
        header_lines, base_lines, tail_lines, seq = load_ps(full_path)
        cur_len = len(seq)
        if sequence[start:start+cur_len].replace('T','U') != seq:
            sys.stderr.writelines("Error: Different Sequence!!\n")
            return
        start += cur_len
        if '_5.ps' in psfile:
            outFn = outPrex + '_5p.ps'
        elif '_3.ps' in psfile:
            outFn = outPrex + '_3p.ps'
        else:
            outFn = outPrex + '.ps'
        if mode == 'label':
            print("Write file: "+outFn+"...")
            write_label_ps(header_lines, base_lines, tail_lines, shape_list, outFn, cutofflist=cutofflist)
        elif mode == 'fill':
            sys.stderr.writelines("Sorry: Not applicant now\n")
        else:
            sys.stderr.writelines("Sorry: Not applicant now\n")

def get_rRNA_refseq(target):
    """
    Return reference sequence
    target              -- rRNA target
                           yeast_small, yeast_large, yeast_5S, yeast_smallMito
                           human_small, human_5S, human_smallMito,
                           mouse_small, mouse_5S, mouse_smallMito
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

def write_label_ps(header_lines, base_lines, tail_lines, shape_list, outFn, cutofflist=[0.3,0.5,0.7]):
    """
    header_lines            -- produced by load_ps
    base_lines              -- produced by load_ps
    tail_lines              -- produced by load_ps
    shape_list              -- A list of SHAPE scores
    outFn                   -- Output file name
    cutofflist              -- The color cutoff
    """
    OUT = open(outFn, "w")
    for header_line in header_lines:
        OUT.writelines(header_line)
    for shape,base_line in zip(shape_list,base_lines):
        if shape == 'NULL':
            OUT.writelines("0.51 0.51 0.51 setrgbcolor\n")
        elif float(shape) < cutofflist[0]:
            OUT.writelines("0.00 0.00 0.00 setrgbcolor\n")
        elif float(shape) < cutofflist[1]:
            OUT.writelines("0.10 0.26 0.60 setrgbcolor\n")
        elif float(shape) < cutofflist[2]:
            OUT.writelines("0.93 0.59 0.09 setrgbcolor\n")
        else:
            OUT.writelines("0.71 0.11 0.13 setrgbcolor\n")
        OUT.writelines(base_line)
    for tail_line in tail_lines:
        OUT.writelines(tail_line)
    OUT.close()




