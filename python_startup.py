# scp -P 20000 python_startup.py zhangqf8@166.111.152.116:/Share2/home/zhangqf/lipan/Load_Env
# ssh -tt -p 20000 zhangqf7@166.111.152.116 'scp /Share2/home/zhangqf/lipan/Load_Env/python_startup.py lipan@10.10.91.12:/home/lipan/usr/Load_Env'
# Lipan's .startup.py file


print(f"(run {__file__})")
import getopt, random, os, math, re, sys, time, copy, datetime, importlib, tempfile
def import_module(model_name, from_import_name=None):
    try:
        module = importlib.import_module(model_name)
    except ModuleNotFoundError:
        print( Colors.f(f"{model_name} not found, import failed", fc='red') )
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
        print( Colors.f(f"{from_import_name} not in {model_name}", fc='red') )
        return None
    else:
        return module

global Colors, Structure, Visual, General, Alignment, Covariation, GAP
Colors = import_module("Colors")
Structure = import_module("Structure")
Visual = import_module("Visual")
General = import_module("General")
Alignment = import_module("Alignment")
Covariation = import_module("Covariation")
GAP = import_module("GAP")
reload = import_module("imp", "reload")

if sys.version_info.major == 2:
    import commands
    print("(import commands)")
    getoutput = commands.getoutput
else:
    #print(Colors.f("commands not exists","yellow"))
    import subprocess
    getoutput = subprocess.getoutput

print("(import getopt, random, os, math, re, sys, time, copy, datetime, importlib, tempfile)")
print("(import Colors, Structure, Visual, General, Alignment, Covariation, GAP)")
print("(from imp import reload)")

#random.seed(1216)

from os import environ as env
HOME = env['HOME']
TEMP = tempfile.gettempdir()
global join
join = os.path.join

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
    global Pool
    global partial
    global reload
    global tqdm
    global _pickle
    global IPython
    global clear_output
    global h5py
    global flush
    global plot_decision_regions
    
    getopt = import_module("getopt")
    random = import_module("random")
    os = import_module("os")
    math = import_module("math")
    re = import_module("re")
    sys = import_module("sys")
    time = import_module("time")
    tempfile = import_module("tempfile")
    copy = import_module("copy")
    datetime = import_module("datetime")
    subprocess = import_module("subprocess")
    Pool = import_module("multiprocessing", "Pool")
    partial = import_module("functools", "partial")
    reload = import_module("imp", "reload")
    tqdm = import_module("tqdm", "tqdm")
    _pickle = import_module("_pickle")
    IPython = import_module("IPython")
    clear_output = import_module("IPython.display", "clear_output")
    h5py = import_module("h5py")
    flush = sys.stdout.flush
    plot_decision_regions = import_module("mlxtend.plotting", "plot_decision_regions")

    ### Image
    global sns
    global matplotlib
    global Patch
    global Line2D
    global plt
    global PdfPages
    global venn2, venn3
    
    sns = import_module("seaborn")
    matplotlib = import_module("matplotlib")
    Patch = import_module("matplotlib.patches", "Patch")
    Line2D = import_module("matplotlib.lines", "Line2D")
    plt = import_module("matplotlib.pyplot")
    PdfPages = import_module("matplotlib.backends.backend_pdf", "PdfPages")
    venn2 = import_module("matplotlib_venn", 'venn2')
    venn3 = import_module("matplotlib_venn", 'venn3')
    
    ### Data Science
    global pd
    global numpy, np
    global scipy
    global statsmodels, multicomp
    
    pd = import_module("pandas")
    numpy = import_module("numpy")
    np = import_module("numpy")
    scipy = import_module("scipy")
    statsmodels = import_module("statsmodels")
    multicomp = import_module("statsmodels.sandbox.stats.multicomp")
    
    ### IPySSSA
    global Colors
    global Structure
    global Visual
    global General
    global Seq
    global Cluster
    global Figures
    global Rosetta
    
    Colors = import_module("Colors")
    Structure = import_module("Structure")
    Visual = import_module("Visual")
    General = import_module("General")
    Seq = import_module("Seq")
    Cluster = import_module("Cluster")
    Figures = import_module("Figures")
    try:
        Rosetta = import_module("D3", "Rosetta")
    except RuntimeError:
        print("import Rosetta failed")
        Rosetta = None
    
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
    
    pysam = import_module("pysam")
    pyBigWig = import_module("pyBigWig")

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
    print( Colors.f(f"list_physical_devices GPUs: {gpus}", fc='yellow') )
    if set_visible_gpu:
        tf.config.experimental.set_visible_devices([gpus[visible_gpu_id]], "GPU")
        tf.config.experimental.set_memory_growth(gpus[visible_gpu_id], True)
    
    tf.random.set_seed(1216)






# print("importCommon() to import more common modules")
# def importCommon():
#     ### Standard library
#     global getopt
#     global random
#     global os
#     global math
#     global re
#     global sys
#     global time
#     global tempfile
#     global copy
#     global datetime
#     global subprocess
#     global Pool
#     global partial
#     global reload
#     global tqdm
#     global _pickle
#     global IPython
#     global clear_output
    
#     import getopt
#     import random
#     import os
#     import math
#     import re
#     import sys
#     import time
#     import tempfile
#     import copy
#     import datetime
#     import subprocess
#     from multiprocessing import Pool
#     from functools import partial
#     from imp import reload
#     from tqdm import tqdm
#     import _pickle
#     import IPython
#     from IPython.display import clear_output
    
#     ### Image
#     global sns
#     global matplotlib
#     global plt
#     global PdfPages
#     global venn2, venn3
    
#     import seaborn as sns
#     import matplotlib
#     import matplotlib.pyplot as plt
#     import matplotlib.backends.backend_pdf
#     from matplotlib.backends.backend_pdf import PdfPages
#     from matplotlib_venn import venn2, venn3
    
#     ### Data Science
#     global pd
#     global numpy, np
#     global scipy
#     global statsmodels, multicomp
    
#     import pandas as pd
#     import numpy
#     import numpy as np
#     import scipy
#     import statsmodels
#     import statsmodels.sandbox.stats.multicomp as multicomp
    
#     ### IPySSSA
#     global Colors
#     global Structure
#     global Visual
#     global General
#     global Seq
#     global Cluster
#     global Figures
#     global Rosetta
    
#     import Colors
#     import Structure
#     import Visual
#     import General
#     import Seq
#     import Cluster
#     import Figures
#     try:
#         from D3 import Rosetta
#     except RuntimeError:
#         print("import Rosetta failed")
    
#     ### GAP
#     global GAP
    
#     import GAP
    
#     ### sklearn
#     global sklearn
#     global OneHotEncoder
#     global svm
#     global roc_curve
#     global auc
#     global precision_recall_curve
#     global biweight_midcorrelation
#     global ensemble
#     global metrics
    
#     import sklearn
#     from sklearn.preprocessing import OneHotEncoder
#     from sklearn import svm
#     from sklearn.metrics import roc_curve, auc, precision_recall_curve
#     import sklearn.ensemble as ensemble
#     from astropy.stats import biweight_midcorrelation
#     import sklearn.metrics as metrics
    
#     np.random.seed(1216)
#     random.seed(1216)

# print("import_tf2() to import tensorflow2")
# def import_tf2(set_visible_gpu=False, visible_gpu_id=-1):
#     global tf
#     global keras
#     global layers
#     global v1
#     global optimizers
#     global losses
#     global Sequential
#     global tfa
    
#     import tensorflow as tf
#     from tensorflow import keras
#     from tensorflow.keras import layers
#     import tensorflow.compat.v1 as v1
#     from tensorflow.keras import optimizers
#     from tensorflow.keras import losses
#     from tensorflow.keras import Sequential
#     import tensorflow_addons as tfa
    
#     gpus = tf.config.experimental.list_physical_devices("GPU")
#     print(f"list_physical_devices GPUs: {gpus}")
#     if set_visible_gpu:
#         tf.config.experimental.set_visible_devices([gpus[visible_gpu_id]], "GPU")
#         tf.config.experimental.set_memory_growth(gpus[visible_gpu_id], True)
    
#     tf.random.set_seed(1216)

