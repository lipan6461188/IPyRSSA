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

"""

import os, sys, commands, random

import General

bkill = General.require_exec("bkill", warning="", exception=True)
bsub = General.require_exec("bsub", warning="", exception=True)
bjobs = General.require_exec("bjobs", warning="", exception=True)

def host_name():
    import os
    return os.environ['HOSTNAME']

def get_bsub_id(output_info):
    import re
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
    import os
    CMD = bjobs + " %s" % (job_id, )
    out_info = commands.getstatusoutput(CMD)[1] 
    if "is not found" in out_info:
        return "Not_Found"
    return out_info.split('\n')[1].split()[2]

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
        if host_name() not in ('loginview02', 'mgt01', 'loginview03'):
            print >>sys.stderr, "Please submit a job in loginview02, loginview03 or mgt01"
            return False
        
        if not hasattr(self, 'command'):
            print >>sys.stderr, "Please set command using set_job_command"
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
        
        CMD += """ "%s" """ % ( self.command,  )
        
        return CMD
    
    def submit(self):
        """
        Submit the command if it is ready
        Return the job id
        """
        import commands
        
        if self.has_submit:
            raise NameError("Error: this job has submited")
        
        if not self.ready():
            raise NameError("Error: JOB not ready")
        
        CMD_string = self.get_submit_command()
        self.job_id = int( get_bsub_id(commands.getstatusoutput(CMD_string)[1]) )
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
        import time
        
        if not self.has_submit:
            raise NameError("Error: Please submit the job at first")
        
        while not bjob_finish(self.job_id):
            time.sleep(interval)
    
    def kill(self):
        """
        Kill the job
        """
        if not self.has_finish():
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








