
##########################
##  Using Python 2.7 and IPyRSSA to visualize icSHAPE in terminal
##########################

## Add PYTHONPATH to IPyRSSA at first
## https://github.com/lipan6461188/IPyRSSA
## export PYTHONPATH=[PATH of IPyRSSA]:$PYTHONPATH

## It is recommended to set the background of your terminal to dark black or white.

## Install muscle (https://www.drive5.com/muscle/downloads.htm)

import General, Colors

human_cy_shape = General.load_shape("hek293_cy_icshape.out")['28S']
human_wc_shape = General.load_shape("hek293_wc_icshape.out")['28S']
mouse_cy_shape = General.load_shape("mes_cy_icshape.out")['28S']

human_28S = General.load_fasta("human_28S.fa")['28S']
mouse_28S = General.load_fasta("mouse_28S.fa")['NR_003279.1']

human_28S_dot = General.load_dot("human_28S.dot")['human_28S'][1]

## Stretch the window of the terminal to a sufficient width

linelen = 200 # recommend window width of the terminal > 220

## 1. Visual human icSHAPE
Colors.browse_shape(human_28S, [human_cy_shape, human_wc_shape], dot=human_28S_dot, linelen=200, shape_title_list=['HEK293_28S_cy', 'HEK293_28S_wc'])

## 2. Visual mouse icSHAPE
Colors.browse_shape(mouse_28S, [mouse_cy_shape], linelen=200, shape_title_list=['mES_28S_cy'])

## 3. Align human/mouse sequence and visualize
Colors.browse_multi_shape([human_28S, mouse_28S], [human_cy_shape, mouse_cy_shape], dot=human_28S_dot, linelen=200, shape_title_list=['HEK293_28S_cy', 'mES_28S_cy'])

######## To save the result as PDF file in mac:
##	1) Type command+K to clean the terminal screen;
## 	2) Run Colors.browse_shape or Colors.browse_multi_shape
##  3) Type command+P and save as PDF

######## To save the result as txt file in mac or linux:

OUT = open("human_28S_icSHAPE.txt", "w")
Colors.browse_shape(human_28S, [human_cy_shape, human_wc_shape], dot=human_28S_dot, linelen=200, shape_title_list=['HEK293_28S_cy', 'HEK293_28S_wc'], OUT=OUT)
OUT.close()

## to view: less -SR human_28S_icSHAPE.txt
