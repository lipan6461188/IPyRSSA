#########
#########   Test Cluster.py
#########

sys.path.append("/Share/home/zhangqf8/lipan/python_utils/PyPsBL")

import Cluster

#####################
#  __host_name()
#####################

print Cluster.host_name()

#####################
#  new_job(command, queue="Z-ZQF", cpu=1, job_name="", logFn="", errFn="")
#####################

##### Example 1

command="sleep 10"
job = Cluster.new_job(command)
job.submit()
job.wait()
job.kill()
job.job_status()

##### Example 2

command = r"echo \"Hello World\" | awk -F \" \" '{print \$1\"\$\$\$\"\$2}' > ~/output.txt"

job2 = Cluster.new_job(command)
job2.submit()
job.wait()

########### Example 3 -- Submit to Z-BNODE and give a job name and log/error file

command = """echo "Hello World"; sleep 10"""

job3 = Cluster.new_job(command=command, queue="Z-BNODE", cpu=1, job_name="My job", logFn="~/log",  errFn="~/error")
job3.get_submit_command()

job3.submit()
job3.wait()


