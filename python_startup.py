# scp -P 20000 python_startup.py zhangqf8@166.111.152.116:/Share2/home/zhangqf/lipan/Load_Env
# ssh -tt -p 20000 zhangqf7@166.111.152.116 'scp /Share2/home/zhangqf/lipan/Load_Env/python_startup.py lipan@10.10.91.12:/home/lipan/usr/Load_Env'
# Lipan's .startup.py file

def warn_f(text):
    code = "%d;%d;%d" % (0, 31, 49)
    return "\x1b["+code+"m"+text+"\x1b[0m"
 
print(f"(run {__file__})")
import os, sys, time, re, random, pickle, copy, gzip, io, yaml, logging, configparser, math, shutil, pathlib, tempfile, hashlib, argparse, json, inspect, urllib, collections, subprocess, requests, platform, multiprocessing, importlib, string, code, warnings, concurrent, gc, functools, types, traceback, base64, bz2, ctypes, tarfile, shlex
from queue import PriorityQueue, Queue, deque, LifoQueue
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool
from os.path import exists, join, getsize, isfile, isdir, abspath, basename, realpath, dirname
from os import listdir
from typing import Dict, Union, Optional, List, Tuple, Mapping
from functools import partial
from datetime import datetime
import numpy  as np
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm, trange
from tabulate import tabulate
from argparse import Namespace
# import getopt, random, os, math, re, sys, time, copy, datetime, importlib, tempfile, collections, pickle, io, gzip, json
def import_module(model_name, from_import_name=None):
    try:
        module = importlib.import_module(model_name)
    except ModuleNotFoundError:
        print( warn_f(f"{model_name} not found, import failed") )
        return None
    except ImportError:
        if warn: print( warn_f(f"{model_name} import error, import failed") )
        return None
    if from_import_name is not None:
        try:
            return module.__dict__[from_import_name]
        except KeyError:
            err_code = 1
        try:
            return getattr(module, from_import_name)
        except AttributeError:
            err_code = 2
        print( warn_f(f"{from_import_name} not in {model_name}") )
        return None
    else:
        return module

global Colors, Structure, Visual, General, Alignment, Cluster, GAP, FFindex, Covariation
Colors = import_module("Colors")
Structure = import_module("Structure")
Visual = import_module("Visual")
General = import_module("General")
Alignment = import_module("Alignment")
Cluster = import_module("Cluster")
if Colors is None:
    Colors = import_module("IPyRSSA", "Colors")
    Structure = import_module("IPyRSSA", "Structure")
    Visual = import_module("IPyRSSA", "Visual")
    General = import_module("IPyRSSA", "General")
    Alignment = import_module("IPyRSSA", "Alignment")
    Cluster = import_module("IPyRSSA", "Cluster")


f = Colors.f
print_tensor_info, print_dict = General.print_tensor_info, General.print_dict
GAP = import_module("GAP")
reload = import_module("imp", "reload")
FFindex = import_module("FFindex", "FFindex")

if sys.version_info.major == 2:
    import commands
    print("(import commands)")
    getoutput = commands.getoutput
else:
    #print(Colors.f("commands not exists","yellow"))
    import subprocess
    getoutput = subprocess.getoutput

print("(import os, sys, time, re, random, pickle, copy, gzip, io, configparser, math, shutil, pathlib, tempfile, hashlib, argparse, json, inspect, urllib, collections, subprocess, requests, platform, multiprocessing, importlib)")
print("(import Colors, Structure, Visual, General, Alignment, Covariation, GAP)")
print("(from imp import reload)")

#random.seed(1216)

from os import environ as env
HOME = env['HOME']
TEMP = tempfile.gettempdir()
STDOUT = sys.stdout
STDERR = sys.stderr

if not hasattr(np, 'object'):
    np.object = np.object_

if not hasattr(np, 'int'):
    np.int = np.int64

if not hasattr(np, 'float'):
    np.float = np.float64

if not hasattr(np, 'bool'):
    np.bool = np.bool_

print("importCommon() to import more common modules")
def importCommon():
    ### Standard library
    global getopt
    global random
    global os
    global math
    global re
    global sys
    global time
    global tempfile
    global copy
    global datetime
    global subprocess
    #global Pool
    global partial
    global reload
    global _pickle
    global IPython
    global clear_output
    global h5py
    global flush
    global plot_decision_regions
    
    getopt                  = import_module("getopt")
    random                  = import_module("random")
    os                      = import_module("os")
    math                    = import_module("math")
    re                      = import_module("re")
    sys                     = import_module("sys")
    time                    = import_module("time")
    tempfile                = import_module("tempfile")
    copy                    = import_module("copy")
    datetime                = import_module("datetime")
    subprocess              = import_module("subprocess")
    #Pool                    = import_module("multiprocessing", "Pool")
    partial                 = import_module("functools", "partial")
    reload                  = import_module("imp", "reload")
    _pickle                 = import_module("_pickle")
    IPython                 = import_module("IPython")
    clear_output            = import_module("IPython.display", "clear_output")
    h5py                    = import_module("h5py")
    flush                   = sys.stdout.flush
    plot_decision_regions   = import_module("mlxtend.plotting", "plot_decision_regions")

    ### Image
    global sns
    global matplotlib
    global Patch
    global Line2D
    global plt
    global PdfPages
    global venn2, venn3
    
    sns         = import_module("seaborn")
    matplotlib  = import_module("matplotlib")
    Patch       = import_module("matplotlib.patches", "Patch")
    Line2D      = import_module("matplotlib.lines", "Line2D")
    plt         = import_module("matplotlib.pyplot")
    PdfPages    = import_module("matplotlib.backends.backend_pdf", "PdfPages")
    venn2       = import_module("matplotlib_venn", 'venn2')
    venn3       = import_module("matplotlib_venn", 'venn3')
    
    # plt.rc('font', family='Helvetica')
    plt.rcParams['pdf.fonttype'] = 42
    
    ### Data Science
    global pd, tabulate
    global numpy, np
    global scipy
    global statsmodels, multicomp
    
    pd          = import_module("pandas")
    tabulate    = import_module("tabulate", "tabulate")
    numpy       = import_module("numpy")
    np          = import_module("numpy")
    scipy       = import_module("scipy")
    statsmodels = import_module("statsmodels")
    multicomp   = import_module("statsmodels.sandbox.stats.multicomp")
    
    ### IPySSSA
    #global Colors
    #global Structure
    #global Visual
    #global General
    #global Seq
    #global Cluster
    #global Figures
    #global Rosetta
    
    #Colors = import_module("Colors")
    #Structure = import_module("Structure")
    #Visual = import_module("Visual")
    #General = import_module("General")
    #Seq = import_module("Seq")
    #Cluster = import_module("Cluster")
    #Figures = import_module("Figures")
    #try:
    #    Rosetta = import_module("D3", "Rosetta")
    #except RuntimeError:
    #    print("import Rosetta failed")
    #    Rosetta = None
    
    ### GAP
    global GAP
    
    GAP = import_module("GAP")
    
    ### sklearn
    global sklearn
    global OneHotEncoder
    global svm
    global roc_curve
    global auc
    global precision_recall_curve
    global biweight_midcorrelation
    global ensemble
    global metrics
    
    sklearn = import_module("sklearn")
    OneHotEncoder = import_module("sklearn.preprocessing", "OneHotEncoder")
    svm = import_module("sklearn.svm")
    roc_curve = import_module("sklearn.metrics", "roc_curve")
    auc = import_module("sklearn.metrics", "auc")
    precision_recall_curve = import_module("sklearn.metrics", "precision_recall_curve")
    ensemble = import_module("sklearn.ensemble")
    biweight_midcorrelation = import_module("astropy.stats", "biweight_midcorrelation")
    metrics = import_module("sklearn.metrics")
    
    ### Biopackage
    global pysam
    global pyBigWig
    global cmap_parse
    global cmap_write
    global mrcfile
    global SVDSuperimposer
    
    pysam      = import_module("pysam")
    pyBigWig   = import_module("pyBigWig")
    cmap_parse = import_module("cmapPy.pandasGEXpress.parse", "parse")
    cmap_write = import_module("cmapPy.pandasGEXpress.write_gctx", "write")
    mrcfile    = import_module('mrcfile')
    SVDSuperimposer = import_module("Bio.SVDSuperimposer", 'SVDSuperimposer')

    ### Random Seed
    #np.random.seed(1216)
    #random.seed(1216)

print("import_tf2() to import tensorflow2")
def import_tf2(set_visible_gpu=False, visible_gpu_id=-1):
    global tf
    global keras
    global layers
    global optimizers
    global losses
    global Sequential
    global tfa
    global v1
    
    tf = import_module("tensorflow")
    keras = import_module("tensorflow", "keras")
    layers = import_module("tensorflow.keras.layers")
    optimizers = import_module("tensorflow.keras.optimizers")
    losses = import_module("tensorflow.keras.losses")
    Sequential = import_module("tensorflow.keras", 'Sequential')
    tfa = import_module("tensorflow_addons")
    v1 = import_module("tensorflow.compat.v1")

    gpus = tf.config.experimental.list_physical_devices("GPU")
    print( warn_f(f"list_physical_devices GPUs: {gpus}") )
    if set_visible_gpu:
        tf.config.experimental.set_visible_devices([gpus[visible_gpu_id]], "GPU")
        tf.config.experimental.set_memory_growth(gpus[visible_gpu_id], True)
    
    tf.random.set_seed(1216)

print("it() to import pytorch")
def it():
    global torch
    global nn
    global F
    global pytree
    global Dataset
    global DataLoader
    global checkpoint
    
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    import torch.utils._pytree as pytree
    from torch.utils.checkpoint import checkpoint
    from torch.utils.data import Dataset, DataLoader

print("it() to import huggingface")
def hf():
    global AutoModelForCausalLM
    global AutoTokenizer
    global AutoConfig
    global safetensors
    global diffusers
    
    from transformers import AutoModelForCausalLM, AutoTokenizer, AutoConfig
    import safetensors
    import diffusers

print("import_pdb() to import PDB related functions")
def import_pdb():
    global pdb_data
    global pdb_features
    global af2_features
    global DensityInfo
    global AF2_Prot_Module
    global AF2_Protein
    global rc
    global Kalign
    
    pdb_data        = import_module("pdb_data")
    pdb_features    = import_module("pdb_features")
    af2_features    = import_module("af2_features")
    
    if pdb_data is not None:
        pdb_data.init_metadata()
    
    DensityInfo     = import_module("cryonet.libs.structure.density", "DensityInfo")
    AF2_Prot_Module = import_module("alphafold.common", "protein")
    AF2_Protein     = import_module("alphafold.common.protein")
    rc              = import_module("alphafold.common", "residue_constants")
    Kalign          = import_module("alphafold.data.tools", "kalign")


def str2date(date_str, fmt='%Y-%m-%d'):
    from datetime import datetime
    return datetime.strptime(date_str, fmt)

def date2str(date, fmt='%Y-%m-%d'):
    return date.strftime(fmt)

