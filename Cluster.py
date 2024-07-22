#-*- coding:utf-8 -*-
"""

This module can manage LSF task submission in python environment.

########### Example 1 -- Quick start

command="sleep 10"
job = new_job(command)
job.submit()
job.wait() # wait this job

job.kill() # kill this job

job.job_status() # see job status

########### Example 2 -- A complex command

#####   give raw string here   ######
command = r"echo \\"Hello World\\" | awk -F \\" \\" '{print \\$1\\"\\$\\$\\$\\"\\$2}' > ~/output.txt"

job2 = new_job(command)
job2.get_submit_command()

job2.submit()
job.wait() # wait this job

########### Example 3 -- Submit to Z-BNODE and give a job name and log/error file

command = r"echo \\"Hello World\\"; sleep 10"

job = new_job(command=command, queue="Z-BNODE", cpu=1, job_name="My job", log_file="~/log",  error_file="~/error")
job.get_submit_command()

job.submit()
job.wait()

########### Example 4 -- job2 will be excute when job1 done

job1 = Cluster.new_job(command="sleep 5", job_name="job1")
job2 = Cluster.new_job(command="sleep 2", job_name="job1")
job2.set_job_depends(job1) # you can use list for multiple process
job1.submit() # job1 must be sumbit before job2
job2.submit()


"""

import os, sys, random, re, time, subprocess, shutil, tempfile, gc, signal
from . import General, Colors
from typing import List, Optional
from multiprocessing import Pool
from tqdm.auto import trange, tqdm
from os.path import abspath, join, exists, realpath
from pathlib import Path

if 'getstatusoutput' in dir(subprocess):
    from subprocess import getstatusoutput
else:
    from commands import getstatusoutput

try:
    bkill = General.require_exec("bkill", warning="", exception=False)
    bsub = General.require_exec("bsub", warning="", exception=False)
    bjobs = General.require_exec("bjobs", warning="", exception=False)
except NameError:
    #print(Colors.f("Error: bkill not found in PATH",fc='red'))
    pass

try:
    srun = General.require_exec("srun", warning="", exception=False)
except NameError:
    #print(Colors.f("Error: srun not found in PATH",fc='red'))
    pass

def host_name():
    return os.environ['HOSTNAME']

def get_bsub_id(output_info):
    return int(re.findall("\\<(\\d+)\\>", output_info)[0])

def bjob_finish(job_id):
    """
    job_id              -- Job id
    
    Return True if the job had finished.
    """
    if job_status(job_id) in ('Not_Found', 'DONE', 'EXIT'):
        return True
    else:
        return False

def kill_job(job_id):
    """
    job_id              -- Job id
    
    Kill the job
    """
    CMD = bkill + " %s > /dev/null" % (job_id)
    os.system(CMD)

def job_status(job_id):
    """
    job_id              -- Job id
    
    Return one of Not_Found, DONE, RUN, PEND, EXIT
    """
    CMD = bjobs + " %s" % (job_id, )
    out_info = getstatusoutput(CMD)[1] 
    if "is not found" in out_info:
        return "Not_Found"
    return out_info.split('\n')[1].split()[2]

def compile_depence_string(job_depends):
    """
    job_depends             -- Job id/JOB_HANDLE object list to wait beform running
    """
    depence_string = ""
    for job in job_depends:
        if isinstance(job, JOB_HANDLE):
            job_id = job.job_id
        elif isinstance(job, int):
            job_id = job
        else:
            raise RuntimeError("job_depends should contain JOB_HANDLE objects and ints")
        depence_string += "&&done(%s)" % (job_id, )
    
    return "\""+depence_string.lstrip('&&')+"\""

class JOB_HANDLE:
    def __init__(self, queue, cpu=1, job_name=""):
        """
        queue               -- Queue name: Z-ZQF, Z-BNODE, Z-HNODE
        cpu                 -- CPU number
        job_name            -- Give a name for job
        
        Return a handle of a new job
        """
        assert queue in ('Z-ZQF', 'Z-BNODE', 'Z-HNODE')
        
        self.queue = queue
        self.cpu = cpu
        self.job_name = job_name
        self.has_submit = False
    
    def set_job_command(self, command):
        """
        command             -- Shell command to submit
        
        Set the command to submit
        """
        self.command = command
    
    def set_job_depends(self, job_depends):
        """
        job_depends         -- Job id/JOB_HANDLE object list to wait beform running
        """
        if type(job_depends) == list:
            self.job_depends = job_depends
        elif type(job_depends) == JOB_HANDLE:
            self.job_depends = [ job_depends ]
        elif type(job_depends) == int:
            self.job_depends = [ job_depends ]
        else:
            raise RuntimeError("job_depends should be list or JOB_HANDLE or int")
    
    def set_log_file(self, logFn):
        """
        logFn            -- A file to save the content from standard output
        """
        self.logFn = logFn
    
    def set_error_file(self, errFn):
        """
        errFn            -- A file to save the content from standard error output
        """
        self.errFn = errFn
    
    def ready(self):
        """
        Return True if the job have ready
        """
        if host_name() not in ('loginview02', 'mgt01', 'loginview03', 'ZIO01'):
            sys.stderr.writelines("Please submit a job in loginview02, loginview03 or mgt01 or ZIO01\n")
            return False
        
        if not hasattr(self, 'command'):
            sys.stderr.writelines("Please set command using set_job_command\n")
            return False
        
        return True
    
    def get_submit_command(self):
        """
        Return the Shell command if it is ready
        """
        if not self.ready():
            return ""
        
        CMD = bsub + " -q %s -n %s -J \"%s\"" % (self.queue, self.cpu, self.job_name)
        
        if hasattr(self, 'logFn'):
            CMD += " -o %s" % (self.logFn, )
        
        if hasattr(self, 'errFn'):
            CMD += " -e %s" % (self.errFn, )
        
        if hasattr(self, 'job_depends'):
            CMD += " -w " + compile_depence_string(self.job_depends)
        
        CMD += """ "%s" """ % ( self.command,  )
        
        return CMD
    
    def submit(self, verbose=False):
        """
        Submit the command if it is ready
        Return the job id
        """
        
        if self.has_submit:
            raise NameError("Error: this job has submited")
        
        if not self.ready():
            raise NameError("Error: JOB not ready")
        
        CMD_string = self.get_submit_command()
        if verbose:
            tqdm.write(CMD_string)
        self.job_id = int( get_bsub_id(getstatusoutput(CMD_string)[1]) )
        self.has_submit = True
        
        return self.job_id
    
    def has_finish(self):
        """
        Return True if the job has finished
        """
        if not self.has_submit:
            raise NameError("Error: Please submit the job at first")
        
        return bjob_finish(self.job_id)
    
    def job_status(self):
        """
        Return one of Not_Found, DONE, RUN, PEND, EXIT
        """
        if not self.has_submit:
            raise NameError("Error: Please submit the job at first")
        
        return job_status(self.job_id)
    
    def wait(self, interval=1):
        """
        interval            -- Seconds for each search
        
        Wait the job to finish
        """
        
        if not self.has_submit:
            raise NameError("Error: Please submit the job at first")
        
        while not bjob_finish(self.job_id):
            time.sleep(interval)
    
    def kill(self):
        """
        Kill the job
        """
        kill_job(self.job_id);
    
    def __del__(self):
        if self.has_submit:
            self.kill()

def new_job(command, queue="Z-ZQF", cpu=1, job_name="", logFn="", errFn=""):
    """
    command             -- Shell command to submit
    queue               -- Queue name: Z-ZQF, Z-BNODE, Z-HNODE
    cpu                 -- CPU number
    job_name            -- Give a name for job
    logFn               -- A file to save the content from standard output
    errFn               -- A file to save the content from standard error output
    
    Return a JOB_HANDLE object
    """
    if queue not in ('Z-ZQF', 'Z-BNODE', 'Z-HNODE'):
        raise NameError("Error: queue must be one of Z-ZQF/Z-BNODE/Z-HNODE")
    
    if job_name == "":
        job_name = command
    
    job = JOB_HANDLE(queue, cpu, job_name)
    job.set_job_command(command)
    if logFn != "":
        job.set_log_file(logFn)
    
    if errFn != "":
        job.set_error_file(errFn)
    
    return job


def slurm_build_dep_chain(cmd_list: List[str], 
                    name_list: Optional[List[str]] = None, 
                    num_chain: int = 1, 
                    cpu: int = 1,
                    gpu: int = 0,
                    partition: Optional[str] = None,
                    node: Optional[str] = None,
                    job_pref: Optional[str] = None, 
                    sleep: float = 2, 
                    shuffle: bool = False):
    """
    Run commands with dependence chain
    
    Parameters
    -------------------
    cmd_list: List of bash commands
    name_list: List of job name postfix
    num_chain: Number of dependence chain
    cpu: Number of CPU
    gpu: Number of GPU
    partition: Partition name
    node: Node name
    job_pref: Job name prefix, random str if not specified
    sleep: Sleep time
    shuffle: If shuffle all jobs

    Return
    -------------------
    chain_job_id_list: List[List[int]]
        Slurm job ID
    """
    import os, sys, subprocess, time, re, random, string, shutil
    import numpy as np
    
    assert shutil.which("squeue") is not None

    num_jobs = len(cmd_list)
    assert num_chain >= 1
    assert num_chain <= num_jobs
    if name_list is not None:
        assert num_jobs == len(name_list)
    else:
        name_list = [ str(idx) for idx in np.arange(num_jobs) ]
    
    if job_pref is None:
        job_pref = "".join([ random.choice(string.ascii_letters) for _ in range(5) ])
    
    def get_last_pid():
        USER = os.environ['USER']
        cmd = f"squeue -o \"%.18i %.9P %.28j\" -u {USER} --states=all -S V | grep {job_pref} | tail -n 1 | awk '{{print $1}}'"
        pid = subprocess.getoutput(cmd)
        if len(pid) > 0:
            pid = int(pid)
            return pid
        else:
            return None
    
    if shuffle:
        index = np.arange(num_jobs)
        np.random.shuffle(index)
        cmd_list  = [ cmd_list[d] for d in index ]
        name_list = [ name_list[d] for d in index ]
    
    chain_size = int(np.ceil(num_jobs / num_chain))
    
    all_chain_pid_list = []
    chain_pid_list = []
    
    pid = None
    for idx in range(num_jobs):
        cmd  = cmd_list[idx]
        name = name_list[idx]
        srun_cmd = f"srun -J {job_pref}_{name} -c {cpu} --gres=gpu:{gpu}"
        if partition:
            srun_cmd += f" -p {partition}"
        if node:
            srun_cmd += f" -w {node}"
        if pid:
            srun_cmd += f" -d afterany:{pid} {cmd} &"
        else:
            srun_cmd += f" {cmd} &"
        tqdm.write(srun_cmd)
        os.system(srun_cmd)
        if pid is None:
            time.sleep(sleep*2)
        else:
            time.sleep(sleep)
        pid = get_last_pid()
        chain_pid_list.append(pid)
        
        if (idx + 1) % chain_size == 0:
            pid = None
            all_chain_pid_list.append(chain_pid_list)
            chain_pid_list = []
    
    if len(chain_pid_list) > 0:
        all_chain_pid_list.append(chain_pid_list)
     
    return all_chain_pid_list

def print_slurm_dep_chain():
    """
    Print depedence chain
    """
    assert shutil.which('squeue') is not None
    import os, subprocess
    
    cmd = f"squeue -u {os.environ['USER']} | sed '1,1d' | awk '{{print $1}}'"
    job_id_list = [ int(job_id) for job_id in subprocess.getoutput(cmd).strip().split() ]
    
    def get_job_dep(job_id):
        job_id = int(job_id)
        dep = None
        for line in subprocess.getoutput(f"scontrol show job {job_id}").split('\n'):
            if 'Dependency=' in line:
                line = line.strip()
                idx = line.index('Dependency=')
                dep = line[idx+11:].strip()
                if len(dep) == 0 or dep == '(null)':
                    return None
                assert ':' in dep, dep
                dep = dep.split(':')[1]
                break
        if dep is not None:
            dep = int(dep)
        return dep
    
    job_dep_list = [ get_job_dep(job_id) for job_id in job_id_list ]
    
    chains = []
    for job_idx, job_id in enumerate(job_id_list):
        if job_id not in job_dep_list:
            chain = [job_id]
            while job_dep_list[job_idx] is not None:
                chain.append(job_dep_list[job_idx])
                job_idx = job_id_list.index(job_dep_list[job_idx])
            chains.append(chain[::-1])
    
    for ch_idx, chain in enumerate(chains):
        print( f"Chain {ch_idx}:", "->".join([str(ch) for ch in chain]) )

def get_slurm_node_resource(node):
    output = subprocess.getoutput(f"scontrol show nodes {node}")
    for line in output.split('\n'):
        line = line.strip()
        if line.startswith('CfgTRES'):
            total_cpu_num = int(re.findall(r"CfgTRES=cpu=(\d+)", line)[0])
            try:
                total_gpu_num = int(re.findall(r"gres/gpu\S*=(\d+)", line)[0])
            except:
                total_gpu_num = 0
        elif line.startswith('AllocTRES'):
            try:
                used_cpu_num = int(re.findall(r"AllocTRES=cpu=(\d+)", line)[0])
            except:
                used_cpu_num = 0
            try:
                used_gpu_num = int(re.findall(r"gres/gpu\S*=(\d+)", line)[0])
            except:
                used_gpu_num = 0
    return {
        'total_cpu_num': total_cpu_num,
        'total_gpu_num': total_gpu_num,
        'used_cpu_num': used_cpu_num,
        'used_gpu_num': used_gpu_num
    }

def get_slurm_resource():
    """
    Get SLURM Resource
    
    Return
    ------------
    partition2nodeResource: Dict of partition -> node -> resource dict
    """
    assert shutil.which('scontrol') is not None
    partition2nodes = {}
    for line in os.popen("scontrol show node"):
        if line.startswith("NodeName="):
            node = re.findall(r"NodeName=([\d\w]+)", line)[0]
        if 'Partitions=' in line:
            partition = re.findall(r"Partitions=([\d\w]+)", line)[0]
            partition2nodes[partition] = partition2nodes.get(partition, []) + [node]
    partition2nodeResource = {}
    for partition, nodes in partition2nodes.items():
        partition2nodeResource[partition] = {}
        resource = {}
        for node in nodes:
            node_resource = get_slurm_node_resource(node)
            #for it in resource_new:
            #    resource[it] = resource.get(it, 0) + resource_new[it]
            partition2nodeResource[partition][node] = node_resource
    return partition2nodeResource

def get_slurm_used_gpu_dict():
    """
    Return dict of used GPUs --  { 'gpu01': [1,2,3], ... }
    """
    info = subprocess.getoutput( "scontrol show job -d | grep \"gpu(IDX:\"" )
    info = [ line.strip() for line in info.split('\n') ]
    
    def parse_line_info(line):
        """ 
        Input can be 
            Nodes=gpu05 CPU_IDs=2-3 Mem=0 GRES_IDX=gpu(IDX:1-4)
        Or
            Nodes=gpu[04-07] CPU_IDs=0-47 Mem=0 GRES=gpu(IDX:0-7)
        Or
            Nodes=gpu[01-02,04,06] CPU_IDs=0-55 Mem=0 GRES=gpu(IDX:0-7)
        """
        gpu_id = re.findall(r"Nodes=([\w\d\[\]\-\,]+)", line)[0]
        gpu_devices = re.findall(r"gpu\(IDX:(.*)\)", line)[0]
        devices = []
        for device_reg in gpu_devices.split(','):
            if '-' in device_reg:
                start, end = device_reg.split('-')
                start, end = int(start), int(end)
                for i in range(start,end+1):
                    devices.append(i)
            else:
                devices.append(int(device_reg))
        gpu_id_list = [  ]
        if '[' in gpu_id or ']' in gpu_id:
            assert ']' in gpu_id and ']' in gpu_id
            s = gpu_id.find('[')
            e = gpu_id.find(']')
            prefix = gpu_id[ : s ]
            gpu_range_list = gpu_id[ s+1: e ].split(',')
            for gpu_range in gpu_range_list:
                if '-' in gpu_range:
                    gpu_s_str, gpu_e_str = gpu_range.split('-')
                    gpu_s, gpu_e = int(gpu_s_str), int(gpu_e_str)
                    for i in range(gpu_s, gpu_e+1):
                        i = str(i)
                        single_gpu_id = prefix + i.zfill( len(gpu_s_str) )
                        gpu_id_list.append( single_gpu_id )
                else:
                    single_gpu_id = prefix + gpu_range
                    gpu_id_list.append( single_gpu_id )
        else:
            gpu_id_list = [ gpu_id ]
        return gpu_id_list, devices
    
    gpu_used_devs = {}
    for line in info:
        if line == "":
            continue
        gpu_id_list, gpu_devices = parse_line_info(line)
        for gpu_id in gpu_id_list:
            gpu_used_devs[gpu_id] = gpu_used_devs.get(gpu_id, []) + gpu_devices
    
    for gpu_id in gpu_used_devs:
        gpu_used_devs[gpu_id].sort()
    
    return gpu_used_devs

def print_slurm_resource():
    """
    Print SLURM resouce to screen
    """
    assert shutil.which('scontrol') is not None
    partition2nodeResource = get_slurm_resource()
    # print(f"{'Partition':<19s}{'GPU Free/Total'}{' '*20}{'CPU Free/Total'}{' '*20}Nodes")
    print(f"{'Partition':<19s}{'Node':<19s}{'GPU Free/Total'}{' '*20}{'CPU Free/Total'}{' '*20}")
    for partition in partition2nodeResource:
        nodeResource = partition2nodeResource[partition]
        for node in nodeResource:
            res = nodeResource[node]
            cpu_tol, cpu_used = res['total_cpu_num'], res['used_cpu_num']
            gpu_tol, gpu_used = res['total_gpu_num'], res['used_gpu_num']
            print(f"{partition:<19s}{node:<19s}{gpu_tol-gpu_used:7d}/{gpu_tol:<7d}{' '*20}{cpu_tol-cpu_used:7d}/{cpu_tol:<7d}{' '*20}")
            
            # res = partition2resource[partition]
            # nodes = partition2nodes[partition]
            # nodes = ",".join(nodes)
            # cpu_tol, cpu_used = res['total_cpu_num'], res['used_cpu_num']
            # gpu_tol, gpu_used = res['total_gpu_num'], res['used_gpu_num']
            # print(f"{partition:<19s}{gpu_tol-gpu_used:7d}/{gpu_tol:<7d}{' '*20}{cpu_tol-cpu_used:7d}/{cpu_tol:<7d}{' '*20}{nodes}")

##################################
#### Multi-GPU for Shell Commands
##################################

def run_cmd_auto_gpu(cmd, gpu_list, state_dir_name):
    from tqdm.auto import trange, tqdm
    gpu_dev = None
    gpu_exists_file = None
    time.sleep( round(random.random(), 3)*2 )
    for gpu_idx in gpu_list:
        gpu_exists_file = os.path.join(state_dir_name, str(gpu_idx))
        if not os.path.exists(gpu_exists_file):
            open(gpu_exists_file, 'a').close()
            while not os.path.exists(gpu_exists_file):
                time.sleep(0.5)
            gpu_dev = gpu_idx
            break
    assert gpu_dev is not None, f"Expect one of gpu devices {gpu_list} exists in {state_dir_name}"
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_dev)
    # cmd = f"CUDA_VISIBLE_DEVICES={gpu_dev} "+cmd
    tqdm.write(f'CUDA_VISIBLE_DEVICES: {gpu_dev} {cmd}')
    os.system(cmd)
    if os.path.exists(gpu_exists_file):
        os.remove(gpu_exists_file)
    else:
        tqdm.write(f"Warning: {gpu_exists_file} not exists")

class Multi_GPU_CMD:
    """
    cmd_list = ['sleep 10'] * 10
    runner = Multi_GPU_CMD([1,2,3], cmd_list)
    runner.run()
    """
    def __init__(self, gpu_list, cmd_list, state_dir_name=None):
        from tqdm.auto import trange, tqdm
        self.gpu_list = gpu_list
        if state_dir_name is None:
            self.state_dir_name = tempfile.mkdtemp(prefix='multigpu_')
            self.state_dir_type = 'generated'
        else:
            self.state_dir_name = state_dir_name
            self.state_dir_type = 'specified'
        tqdm.write(f"state_dir_name: {self.state_dir_name}")
        self.cmd_list = cmd_list
    
    def run(self, disable_tqdm=True):
        assert isinstance(self.cmd_list, (list, tuple)), self.cmd_list
        for cmd in self.cmd_list:
            assert isinstance(cmd, str), cmd
        assert isinstance(self.gpu_list, list), self.gpu_list
        for gpu in self.gpu_list:
            assert isinstance(gpu, int), gpu
        num_proc = len(self.gpu_list)
        with Pool(num_proc) as pool:
            h_list = []
            for cmd in self.cmd_list:
                h = pool.apply_async(run_cmd_auto_gpu, (cmd, self.gpu_list, self.state_dir_name))
                h_list.append([cmd, h])
            for cmd, h in tqdm(h_list, disable=disable_tqdm, dynamic_ncols=True):
                print(cmd)
                _ = h.get()
    
    def __del__(self):
        if self.state_dir_type == 'generated':
            shutil.rmtree(self.state_dir_name)

##################################
#### Multi-GPU for Python Function
##################################

def run_pyfunc_auto_gpu(pyfunc, args, kwargs, gpu_param_name, gpu_list, state_dir_name):
    
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}
    if gpu_param_name is not None:
        assert gpu_param_name not in kwargs, f"{gpu_param_name} should not in kwargs"
    
    gpu_dev = None
    gpu_exists_file = None
    time.sleep( round(random.random(), 3)*2 )
    for gpu_idx in gpu_list:
        gpu_exists_file = os.path.join(state_dir_name, str(gpu_idx))
        if not os.path.exists(gpu_exists_file):
            open(gpu_exists_file, 'a').close()
            while not os.path.exists(gpu_exists_file):
                time.sleep(0.5)
            gpu_dev = gpu_idx
            break
    assert gpu_dev is not None, f"Expect one of gpu devices {gpu_list} exists in {state_dir_name}"
    if gpu_param_name is not None:
        kwargs[gpu_param_name] = gpu_dev
    else:
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_dev)
    output = pyfunc(*args, **kwargs)
    if os.path.exists(gpu_exists_file):
        os.remove(gpu_exists_file)
    else:
        tqdm.write(f"Warning: {gpu_exists_file} not exists")
    return output

class Multi_GPU_PYFUNC:
    """
    args_list = [[10]] * 10
    kwargs_list = None
    runner = Multi_GPU_PYFUNC([1,2,3], time.sleep, None, args_list, kwargs_list)
    runner.run()
    """
    def __init__(self, gpu_list, pyfunc, gpu_param_name, args_list, kwargs_list, state_dir_name=None):
        from tqdm.auto import trange, tqdm
        #print('args_list:', args_list)
        #print('kwargs_list:', kwargs_list)
        if args_list is None and kwargs_list is None:
            raise RuntimeError("At least one of args_list and kwargs_list should be sepcified")
        if args_list is None:
            args_list = [[] for _ in range(len(kwargs_list))]
        if kwargs_list is None:
            kwargs_list = [{} for _ in range(len(args_list))]
        assert len(args_list) == len(kwargs_list), f"Expect same length, but got {len(args_list)} and {len(kwargs_list)}"
        
        if state_dir_name is None:
            self.state_dir_name = tempfile.mkdtemp(prefix='multigpu_')
            self.state_dir_type = 'generated'
        else:
            self.state_dir_name = state_dir_name
            self.state_dir_type = 'specified'
        tqdm.write(f"state_dir_name: {self.state_dir_name}")
        self.pyfunc      = pyfunc
        self.gpu_param_name = gpu_param_name
        self.gpu_list    = gpu_list
        self.args_list   = args_list
        self.kwargs_list = kwargs_list
        assert isinstance(self.args_list, (list, tuple)), self.args_list
        assert isinstance(self.kwargs_list, (list, tuple)), self.kwargs_list
        assert isinstance(self.gpu_list, (list, tuple)), self.gpu_list
        for gpu in self.gpu_list:
            assert isinstance(gpu, int), gpu
    
    def run(self, disable_tqdm=True):
        num_proc = len(self.gpu_list)
        output_list = []
        with Pool(num_proc) as pool:
            h_list = []
            for args, kwargs in zip(self.args_list, self.kwargs_list):
                h = pool.apply_async(run_pyfunc_auto_gpu, (self.pyfunc, args, kwargs, self.gpu_param_name, self.gpu_list, self.state_dir_name))
                h_list.append(h)
            for h in tqdm(h_list, disable=disable_tqdm, dynamic_ncols=True):
                output_list.append(h.get())
        return output_list
    
    def __del__(self):
        if self.state_dir_type == 'generated':
            shutil.rmtree(self.state_dir_name)


###################################               
### Get resource limitations in docker environment
###################################               
                                                  
def get_docker_limit():
    """
    Get Docker limits

    Return
    -----------
    limits: dict
        -- CPU: Number of CPU
        -- MEM: Numver of Mem
    """
    limits = {}

    ### CPU
    from multiprocessing import cpu_count
    if os.path.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") and os.path.exists("/sys/fs/cgroup/cpu/cpu.cfs_period_us"):
        cfs_quota_us   = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
        cfs_period_us  = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
        container_cpus = cfs_quota_us // cfs_period_us
        # For physical machine, the `cfs_quota_us` could be '-1'
        limits['CPU'] = cpu_count() if container_cpus < 1 else container_cpus
    else:
        limits['CPU'] = cpu_count()

    ### MEM
    limits['MEM'] = int(open('/sys/fs/cgroup/memory/memory.limit_in_bytes').read()) // 1024 ** 3
    limits['MEM_Usage'] = int(open('/sys/fs/cgroup/memory/memory.usage_in_bytes').read()) // 1024 ** 3

    return limits

###################################               
### Lock programs to permit only one process can run
###################################     

class ProgLock:
    """
    Semaphore to control only one process can enter the environment
    """
    def __init__(self, tagfile, wait_second=10, verbose=True):
        """
        tagfile: a tag file to label the Semaphore condition. this file must be create before run the program
        wait_second: seconds to wait when tagfile not found
        verbose: print information when tagfile not found
        """
        self.tagfile = tagfile
        self.wait_second = wait_second
        self.verbose = verbose
    
    def __enter__(self):
        from tqdm.auto import trange, tqdm
        while True:
            while not os.path.exists(self.tagfile):
                if self.verbose:
                    tqdm.write(f"Tag file {self.tagfile} not found. wait {self.wait_second}s...")
                time.sleep(self.wait_second)
            try:
                os.remove(self.tagfile)
                break
            except FileNotFoundError:
                tqdm.write(f"Delete file {self.tagfile} failed. wait {self.wait_second}s...")
                time.sleep(self.wait_second)
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        from pathlib import Path
        Path(self.tagfile).touch()


class ParallelTaskManager:
    def __init__(self, input_dir, file_suffix, output_dir, process_fn, lock_file, filter_fn=None):
        """
        input_dir: directory to contain candidate files
        file_suffix: files with this surfix will be runned
        output_dir: output directory
        process_fn: a callable object with two parameters (in_file, out_prefix)
        lock_file: a lock file to prevent conflics
        filter_fn: filter function with parameters (file_name), return True to exec and False to skip
        """
        self.input_dir   = input_dir
        self.file_suffix = file_suffix
        self.output_dir  = output_dir
        self.process_fn  = process_fn
        self.lock_file   = lock_file
        self.filter_fn   = filter_fn if filter_fn is not None else (lambda x: True)
        
        # assert realpath(input_dir) != realpath(output_dir), f"input_dir must be diffrent from output_dir"
        def sigint_handler(sig, frame):
            raise RuntimeError("SIGINT")
        signal.signal(signal.SIGINT, sigint_handler) # 2
    
    def get_candidates(self, verbose=True):
        from tqdm.auto import trange, tqdm
        candidate_files = [f for f in os.listdir(self.input_dir) if f.endswith(self.file_suffix) and self.filter_fn(f)]
        running_files, done_files, err_files, to_run_files = [], [], [], []
        for f in candidate_files:
            f = f[:-len(self.file_suffix)]
            if exists(join(self.output_dir,   f"{f}.done")):
                done_files.append(f)
            elif exists(join(self.output_dir, f"{f}.running")):
                running_files.append(f)
            elif exists(join(self.output_dir, f"{f}.error")):
                err_files.append(f)
            else:
                to_run_files.append(f)
        
        if verbose:
            tqdm.write(Colors.f(f"Total {len(candidate_files)} candidate files, {len(running_files)} running, {len(done_files)} done, {len(err_files)} error, {len(to_run_files)} to run.", 'yellow'))
        
        return running_files, done_files, err_files, to_run_files
    
    def run(self, verbose=True):
        from tqdm.auto import trange, tqdm
        while True:
            gc.collect()
            with ProgLock(self.lock_file):
                
                running_files, done_files, err_files, to_run_files = self.get_candidates(verbose=verbose)
                
                if len(to_run_files) == 0 and len(err_files) == 0:
                    break
                
                if len(err_files) > 0:
                    name     = err_files[0]
                else:
                    name     = to_run_files[0]
                
                in_file          = join(self.input_dir,  name+self.file_suffix)
                out_prefix       = join(self.output_dir, name)
                running_tag_file = join(self.output_dir, f"{name}.running")
                done_tag_file    = join(self.output_dir, f"{name}.done")
                err_tag_file     = join(self.output_dir, f"{name}.error")
                
                if exists(running_tag_file):
                    continue
                tqdm.write(Colors.f(f"Create {running_tag_file}", 'blue'))
                Path(running_tag_file).touch()
                if exists(err_tag_file):
                    os.remove(err_tag_file)
            try:
                # run_main(in_file, out_prefix, n_process=n_process, n_threads=n_threads, jackhmmer_ncpu=jackhmmer_ncpu, min_msa=min_msa, disable_tqdm=disable_tqdm)
                start_time = time.time()
                self.process_fn(in_file, out_prefix)
                exec_time = time.time() - start_time
                if verbose:
                    tqdm.write(Colors.f(f"{name} finished: exec {exec_time:.1f}s", 'green'))
            except Exception as e:
                tqdm.write(Colors.f("===============ERROR===============", 'red'))
                tqdm.write(e)
                os.remove(running_tag_file)
                Path(err_tag_file).touch()
                return -1
            
            os.remove(running_tag_file)
            Path(done_tag_file).touch()
        
        return 0


###################################               
### Manage Torch DDP environment
###################################     

def get_rank():
    """
    Get current rank. 
    Priority:
    - torch.distributed.get_rank()
    - OMPI_COMM_WORLD_RANK
    - RANK
    """
    import torch
    if torch.distributed.is_initialized():
        return torch.distributed.get_rank()
    if 'OMPI_COMM_WORLD_RANK' in os.environ:
        return int(os.environ['OMPI_COMM_WORLD_RANK'])
    if 'RANK' in os.environ:
        return int(os.environ['RANK'])
    return 0

def get_world_size():
    """
    Get current rank. 
    Priority:
    - torch.distributed.get_world_size()
    - OMPI_COMM_WORLD_SIZE
    - WORLD_SIZE
    - torch.cuda.device_count()
    """
    import torch
    if torch.distributed.is_initialized():
        return torch.distributed.get_world_size()
    if 'OMPI_COMM_WORLD_SIZE' in os.environ:
        return int(os.environ['OMPI_COMM_WORLD_SIZE'])
    if 'WORLD_SIZE' in os.environ:
        return int(os.environ['WORLD_SIZE'])
    return torch.cuda.device_count()

def print_rank_0(*args, **kwargs):
    """
    Print only in rank_0
    """
    if get_rank() == 0:
        print(*args, **kwargs)

def init_distributed(verbose=True):
    """
    Initialize distributed environment
    """
    import torch, os
    if 'WORLD_SIZE' in os.environ and 'LOCAL_WORLD_SIZE' in os.environ:
        world_size, local_world_size = int(os.environ['WORLD_SIZE']), int(os.environ['LOCAL_WORLD_SIZE'])
        if verbose:
            print_rank_0(f">>>>>>>>> WORLD_SIZE={world_size}, LOCAL_WORLD_SIZE={local_world_size} <<<<<<<<<<<")
        assert world_size % local_world_size == 0
        num_nodes = world_size // local_world_size
        devices = local_world_size
    elif 'OMPI_COMM_WORLD_SIZE' in os.environ:
        import deepspeed
        deepspeed.init_distributed()
        world_size, local_world_size = int(os.environ['OMPI_COMM_WORLD_SIZE']), int(os.environ['OMPI_COMM_WORLD_LOCAL_SIZE'])
        if verbose:
            print_rank_0(f">>>>>>>>> WORLD_SIZE={world_size}, LOCAL_WORLD_SIZE={local_world_size} <<<<<<<<<<<")
        assert world_size % local_world_size == 0
        num_nodes = world_size // local_world_size
        devices = local_world_size
    else:
        num_nodes = 1
        devices = torch.cuda.device_count()
        world_size = torch.cuda.device_count()
    return num_nodes, devices, world_size

