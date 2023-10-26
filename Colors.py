#-*- coding:utf-8 -*-

fontColors = {
    'black': 30,
    'red': 31,
    'green': 32,
    'yellow': 33,
    'blue': 34,
    'magenta': 35,
    'cyan': 36,
    'lightgray': 37,
    'default': 39,
    'darkgray': 90,
    'lightred': 91,
    'lightgreen': 92,
    'lightyellow': 93,
    'lightblue': 94,
    'lightmagenta': 95,
    'lightcyan': 96,
    'white': 97
}

backgroundColors = {
    'black': 40,
    'red': 41,
    'green': 42,
    'yellow': 43,
    'blue': 44,
    'magenta': 45,
    'cyan': 46,
    'lightgray': 47,
    'default': 49,
    'darkgray': 100,
    'lightred': 101,
    'lightgreen': 102,
    'lightyellow': 103,
    'lightblue': 104,
    'lightmagenta': 105,
    'lightcyan': 106,
    'white': 107
}

formatting = {
    'normal': 0,
    'bold': 1,
    'light': 2,
    'italic': 3,
    'underline': 4,
    'blink': 5,
    'reverse': 7
}

RGB = {
    'red'    : "#f44336",
    'pink'   : "#e91e63",
    'purple' : "#9c27b0",
    'indigo' : "#3f51b5",
    'blue'   : "#2196f3",
    'cyan'   : "#00bcd4",
    'aqua'   : "#00ffff",
    'teal'   : "#009688",
    'green'  : "#4caf50",
    'khaki'  : "#f0e68c",
    'yellow' : "#ffeb3b",
    'amber'  : "#ffc107",
    'orange' : "#ff9800",
    'brown'  : "#795548",
    'gray'   : "#9e9e9e",
    'lime'   : "#cddc39",
    
    'deep_purple' : "#673ab7",
    'light_blue'  : "#87ceeb",
    'light_green' : "#8bc34a",
    'deep_orange' : "#ff5722",
    'blue_gray'   : "#607d8b",
    'dark_gray'   : "#616161",
    'light_gray'  : "#f1f1f1",
    'pale_yellow' : "#ffffcc",
    'pale_red'    : "#ffdddd",
    'pale_blue'   : "#ddffff",
    'pale_green'  : "#ddffdd",
}

import sys

def format(text, fc='red', bc='default', ft='normal'):
    """
    fc                      -- font color
    bc                      -- background color
    ft                      -- formatting
    
    Colors: 
        blue, lightgray, lightyellow, default, darkgray, yellow, lightred, lightcyan, black, lightmagenta, lightblue, cyan, green, magenta, lightgreen, white, red
    Formattings:
        bold, normal, light, blink, italic, underline, reverse
    """
    code = "%d;%d;%d" % (formatting[ft], fontColors[fc], backgroundColors[bc])
    return "\x1b["+code+"m"+text+"\x1b[0m"

def f(text, fc='red', bc='default', ft='normal'):
    """
    fc                      -- font color
    bc                      -- background color
    ft                      -- formatting
    
    Colors: 
        blue, lightgray, lightyellow, default, darkgray, yellow, lightred, lightcyan, black, lightmagenta, lightblue, cyan, green, magenta, lightgreen, white, red
    Formattings:
        bold, normal, light, blink, italic, underline, reverse
    """
    return format(text, fc=fc, bc=bc, ft=ft)

def color_SHAPE(shape_list, cutoff=[0.3, 0.5, 0.7]):
    """
    shape_list              -- A list of SHAPE scores
    cutoff                  -- Cutoff of SHAPE color boundaries.
    
    Transform SHAPE values to color blocks
    """
    color_blocks = ""
    
    for value in shape_list:
        if value == 'NULL':
            color_blocks += f(" ", bc='lightgray')
        else:
            shape = float(value)
            if shape < cutoff[0]:
                color_blocks += f(" ", bc='blue')
            elif shape < cutoff[1]:
                color_blocks += f(" ", bc='cyan')
            elif shape < cutoff[2]:
                color_blocks += f(" ", bc='green')
            else:
                color_blocks += f(" ", bc='red')
    
    return color_blocks

def color_Seq_SHAPE(sequence, shape_list, cutoff=[0.3, 0.5, 0.7]):
    """
    sequence                -- Raw sequence
    shape_list              -- A list of SHAPE scores
    cutoff                  -- Cutoff of SHAPE color boundaries.
    
    Transform seuquence to colorful sequence according to their shape values
    """
    
    assert len(sequence) == len(shape_list)
    
    color_seq = ""
    
    for base, value in zip(sequence,shape_list):
        if value == 'NULL':
            color_seq += f(base, fc='lightgray')
        else:
            shape = float(value)
            if shape < cutoff[0]:
                color_seq += f(base, fc='blue')
            elif shape < cutoff[1]:
                color_seq += f(base, fc='cyan')
            elif shape < cutoff[2]:
                color_seq += f(base, fc='green')
            else:
                color_seq += f(base, fc='red')
    
    return color_seq

def browse_shape(sequence, shape_list_list, linelen=200, dot="", shape_title_list=[], colorcutoff=[0.3, 0.5, 0.7], OUT=sys.stdout):
    """
    sequence                -- Sequence
    shape_list_list         -- [ [shape_list1, shape_list2, shape_list3, ...], [], []... ]
    linelen                 -- Number of bases for each line
    dot                     -- Dot structure
    shape_title_list        -- Title for each shape list
    colorcutoff             -- Cutoff for colorblock
    
    Print/compare shape scores in screen
    """
    
    import General
    
    ### Check the correctness of parameters
    for shape_list in shape_list_list:
        assert len(sequence) == len(shape_list)
    if shape_title_list:
        assert len(shape_list_list) == len(shape_title_list)
    else:
        shape_title_list = [ "" ] * len(shape_list_list)
    if dot:
        assert len(sequence) == len(dot)
    
    ### Print lengend
    red_legend = format("  ", bc="red")+" >%s " % (colorcutoff[2], )
    green_legend = format("  ", bc="green")+" %s-%s " % (colorcutoff[1], colorcutoff[2])
    cyan_legend = format("  ", bc="cyan")+" %s-%s " % (colorcutoff[0], colorcutoff[1])
    blue_legend = format("  ", bc="blue")+" <%s " % (colorcutoff[0], )
    null_legend = format("  ", bc="lightgray")+" NULL "
    OUT.writelines("\n#### Legend\n")
    OUT.writelines(" "*5 + red_legend + green_legend + cyan_legend + blue_legend + null_legend + "\n\n")
    
    ### Calculate AUC
    if dot:
        OUT.writelines("#### AUC\n")
        for head, shape_list in zip(shape_title_list, shape_list_list):
            roc = General.calc_shape_structure_ROC(dot, shape_list, start=0.0, step=0.01, stop=1.0)
            auc = round(General.calc_AUC(roc), 3)
            OUT.writelines("     "+head+"\t"+str(auc)+"\n")
        OUT.writelines("\n")
    
    ### Estimate min head length
    max_title_len = max([len(title) for title in shape_title_list])
    max_seqnum_len = 2*len(str(len(sequence)))+1
    min_head_len = max(max_seqnum_len, max_title_len) + 2
    
    ### Print sequence, structure and shape
    i = 0
    while i<len(sequence):
        end = min(i+linelen, len(sequence))
        head = "%s-%s" % (i+1, end)
        head += " "*(min_head_len-len(head))
        OUT.writelines(head+sequence[i:end]+"\n")
        if dot:
            OUT.writelines(" "*min_head_len+dot[i:end]+"\n")
        for head, shape_list in zip(shape_title_list, shape_list_list):
            head += " "*(min_head_len-len(head))
            OUT.writelines(head+color_SHAPE(shape_list[i:end], cutoff=colorcutoff)+"\n")
        OUT.writelines("\n")
        i += linelen

def browse_multi_shape(sequence_list, shape_list_list, linelen=200, dot="", shape_title_list=[], OUT=sys.stdout):
    """
    sequence_list           -- Sequence list
    shape_list_list         -- [ [shape_list1, shape_list2, shape_list3, ...], [], []... ]
    linelen                 -- Number of bases for each line
    dot                     -- Dot structure of the first sequence
    shape_title_list        -- Title for each sequence/shape
    
    Align and print/compare shape scores in screen
    """
    
    from IPyRSSA import General, Structure
    
    ### Check the correctness of parameters
    assert len(sequence_list) == len(shape_list_list)
    if shape_title_list:
        assert len(sequence_list) == len(shape_title_list)
    else:
        shape_title_list = [ "" ] * len(shape_list_list)
    for sequence, shape_list in zip(sequence_list, shape_list_list):
        assert len(sequence) == len(shape_list)
    if dot:
        assert len(sequence_list[0]) == len(dot)
    
    aligned_seq_list = Structure.multi_alignment(sequence_list, clean=True, verbose=False)
    aligned_shape_list = [ Structure.shape_to_alignSHAPE(raw_shape,aligned_seq) for raw_shape,aligned_seq in zip(shape_list_list, aligned_seq_list) ]
    aligned_dot = ""
    if dot:
        aligned_dot = Structure.dot_to_alignDot(dot, aligned_seq_list[0])
    
    ### Estimate min head length
    max_title_len = max([len(title) for title in shape_title_list])
    max_seqnum_len = max( [2*len(str(len(raw_seq)))+1 for raw_seq in sequence_list] )
    min_head_len = max(max_seqnum_len, max_title_len) + 2
    
    ### Print lengend
    red_legend = format("  ", bc="red")+" >0.7 "
    green_legend = format("  ", bc="green")+" 0.5-0.7 "
    cyan_legend = format("  ", bc="cyan")+" 0.3-0.5 "
    blue_legend = format("  ", bc="blue")+" <0.3 "
    null_legend = format("  ", bc="lightgray")+" NULL "
    OUT.writelines("\n#### Legend\n")
    OUT.writelines(" "*5 + red_legend + green_legend + cyan_legend + blue_legend + null_legend + "\n\n")
    
    ### Print sequence, structure and shape
    i = 0
    aligned_seq_len = len(aligned_seq_list[0])
    while i<aligned_seq_len:
        end = min(i+linelen, aligned_seq_len)
        index = 0
        for title,aligned_seq,aligned_shape in zip(shape_title_list,aligned_seq_list,aligned_shape_list):
            raw_start = len(aligned_seq[:i].replace("-",""))+1
            raw_end = len(aligned_seq[:end].replace("-",""))
            head = "%s-%s" % (raw_start, raw_end)
            head += " "*(min_head_len-len(head))
            OUT.writelines(head+aligned_seq[i:end]+"\n")
            if index == 0:
                if aligned_dot:
                    OUT.writelines(" "*min_head_len+aligned_dot[i:end]+"\n")
            head = title
            head += " "*(min_head_len-len(head))
            OUT.writelines(head+color_SHAPE(aligned_shape[i:end], cutoff=[0.3, 0.5, 0.7])+"\n")
            index += 1
        i += linelen
        OUT.writelines("\n")


def _hex_to_RGB(hex):
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i+2], 16) for i in range(1,6,2)]

def _RGB_to_hex(RGB):
    ''' [255,255,255] -> "#FFFFFF" '''
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])

def _color_dict(gradient):
    ''' Takes in a list of RGB sub-lists and returns dictionary of
        colors in RGB and hex form for use in a graphing function
        defined later on '''
    return {
        "hex":[_RGB_to_hex(RGB) for RGB in gradient],
        "r":[RGB[0] for RGB in gradient],
        "g":[RGB[1] for RGB in gradient],
        "b":[RGB[2] for RGB in gradient]}

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
    # Starting and ending colors in RGB form
    s = _hex_to_RGB(start_hex)
    f = _hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [
            int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
            for j in range(3)
        ]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)
    
    return _color_dict(RGB_list)

