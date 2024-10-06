#-*- coding:utf-8 -*-

import os, sys, subprocess, random, time, re
from IPyRSSA import General

if 'getstatusoutput' in dir(subprocess):
    from subprocess import getstatusoutput
else:
    from commands import getstatusoutput

class GPUProcess(object):
    def __init__(self, gpuid, pid, ptype, pname, memo):
        super(GPUProcess, self).__init__()
        self.gpuid = gpuid
        self.pid = pid
        self.ptype = ptype
        self.pname = pname
        self.memo = memo
    def __str__(self):
        return "__GPU(%s)__PID(%s)__Proc(%s)__" % (self.gpuid, self.gpuid, self.pname)
    def __repr__(self):
        return "__GPU(%s)__PID(%s)__Proc(%s)__" % (self.gpuid, self.gpuid, self.pname)

def get_gpu_processes():
    """
    Get a list of GPUProcess objects
    
    Require: nvidia-smi
    
    """
    NVIDIASMI = General.require_exec("nvidia-smi", exception=False)
    
    lines = getstatusoutput(NVIDIASMI)[1].strip().split('\n')
    sl = 0
    while "Processes" not in lines[sl]:
        sl += 1
    GPUProc_List = []
    for line in lines[sl+3:-1]:
        if "No running processes found" in line: break
        data = line.strip('|').strip().split()
        gpu_id = int(data[0])
        pid = int(data[1])
        ptype = data[2]
        memo = int(data[-1].rstrip('MiB'))
        pname = " ".join(data[3:-1])
        proc = GPUProcess(gpu_id, pid, ptype, pname, memo)
        GPUProc_List.append(proc)
    return GPUProc_List

def get_gpu_list():
    """
    Get a list of available gpu
    
    Require: nvidia-smi
    
    """
    NVIDIASMI = General.require_exec("nvidia-smi", exception=False)
    
    lines = getstatusoutput(NVIDIASMI+" -L")[1].strip().split('\n')
    gpuids = [ int(line.split()[1].rstrip(':')) for line in lines ]
    return gpuids

def get_free_gpus():
    """
    Get a list of GPU id without process run on it
    example:
        [3,4]
    
    """
    gpu_list = set(get_gpu_list())
    gpu_processes = get_gpu_processes()
    busy_gpus = set([ i.gpuid for i in gpu_processes ])
    free_gpus = list( gpu_list - busy_gpus )
    return sorted( free_gpus )

