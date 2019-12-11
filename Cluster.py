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

import os, sys, random, re, time
import General, subprocess, Colors

if 'getstatusoutput' in dir(subprocess):
    from subprocess import getstatusoutput
else:
    from commands import getstatusoutput

try:
    bkill = General.require_exec("bkill", warning="", exception=True)
    bsub = General.require_exec("bsub", warning="", exception=True)
    bjobs = General.require_exec("bjobs", warning="", exception=True)
except NameError:
    print(Colors.f("Error: bkill not found in PATH, Cluster module cannot be used",fc='red'))

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
            print(CMD_string)
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






