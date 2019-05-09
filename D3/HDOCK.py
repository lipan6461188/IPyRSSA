
"""

This module can dock protein and DNA/RNA with HDOCK

from D3 import HDOCK

## 1. Upload a task
proteinPDBFn = "/tmp/Dicer/dicer.pdb"
RNAPDBFn = "/tmp/Dicer/miRNA.pdb"
jobname = "dicer222"
pResi= "324:A, 327:A, 328:A, 331:A, 455:A, 459:A, 463:A, 479:A, 480:A, 481:A, 482:A, 483:A, 788:A, 936:A, 937:A, 954:A, 958:A, 959:A, 960:A, 961:A, 962:A, 963:A, 968:A, 971:A, 972:A, 975:A, 976:A, 1014:A, 1016:A, 1019:A, 1020:A, 1023:A, 1024:A, 1026:A, 1027:A, 1028:A, 1030:A, 1031:A, 1032:A, 1033:A, 1886:A, 1887:A, 1888:A, 1889:A, 1890:A, 1911:A, 1913:A"
dist = "1032:A 73:C 25, 445:A 25:C 35"

jobID, res = HDOCK.upload_HDOCK_job(proteinPDBFn, RNAPDBFn, jobname=jobname, pResi=pResi, dist=dist)

## 2. Wait a task -- very long time
status = HDOCK.get_HDOCK_status(jobID)
while status == 'RUNNING':
    status = HDOCK.get_HDOCK_status(jobID)
    timeleft = HDOCK.guess_HDOCK_time_left(jobID)
    print(status)+"\t"+str(timeleft)+" seconds left"
    time.sleep(5)

## 3. Download all PDB
HDOCK.fetch_HDOCK_results(jobID, "/tmp/all.tar.gz")
HDOCK.fetch_HDOCK_top10_results(jobID, "/tmp/top10.tar.gz")

"""


import os, sys, random, time, re
import threading 
from bs4 import BeautifulSoup
import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

if sys.version_info.major == 2:
    import urllib2
else:
    import urllib.request as urllib2

def upload_HDOCK_job(proteinPDBFn, RNAPDBFn, email=None, jobname=None, pResi=None, nResi=None, dist=None, verbose=True):
    """
    proteinPDBFn            -- Protein PDB filename
    RNAPDBFn                -- RNA or DNA PDB filename
    email                   -- Email address
    jobname                 -- Job name
    pResi                   -- Protein interface residues constraints
                                such as: 324:A, 327:A, 328:A, 331:A
    nResi                   -- Nucleoacid interface residues constraints
                                such as: 324:C, 327:C, 328:C, 331:C
    dist                    -- Maximun distance constraints
                                such as: 1032:A 73:C 25, 445:A 25:C 35
    verbose                 -- Print some information
    
    Return:
        jobID, res
        jobID           -- A ID string
        res             -- Requests Responds object
    
    Exception:
        RuntimeError
    """
    
    files = { 'pdbfile1':open(proteinPDBFn, 'rb'), 'pdbfile2':open(RNAPDBFn, 'rb') }
    param = { 
                'upload': 'Submit',
                'email': "",
                'jobname': "",
                'sitenum1': "",
                'sitenum2': "",
                'restrnum': ""
            }
    
    if email:
        param['email'] = email
    if jobname:
        param['jobname'] = jobname
    if pResi:
        param['sitenum1'] = pResi
    if nResi:
        param['sitenum2'] = nResi
    if dist:
        param['restrnum'] = dist
    
    headers = {'User-Agent': 'Mozilla/5.0'}
    res = requests.post("http://hdock.phys.hust.edu.cn/runhdock.php", data=param, files=files, headers=headers)
    res.raise_for_status()
    soup = BeautifulSoup(res.text,'html.parser')
    jobID = soup.title.string.split()[-1]
    
    if verbose:
        print("JobID: "+str(jobID))
        print("http://hdock.phys.hust.edu.cn/data/"+str(jobID)+"/")
    
    return jobID, res

def is_HDOCK_exists(HDOCK_jobID):
    
    res = requests.get("http://hdock.phys.hust.edu.cn/data/%s/" % (HDOCK_jobID, ))
    if res.status_code == 404:
        return False
    if res.status_code == 200:
        return True
    sys.stderr.writelines("Error: Unknown error in is_HDOCK_exists: "+str(res.status_code)+"\n")
    return False

def is_HDOCK_complete(HDOCK_jobID):
    res = requests.get("http://hdock.phys.hust.edu.cn/data/%s/" % (HDOCK_jobID, ))
    res.raise_for_status()
    soup = BeautifulSoup(res.text, 'html.parser')
    isRunning = bool( soup.find_all(string=re.compile("is RUNNING")) )
    return not isRunning

def guess_HDOCK_time_left(HDOCK_jobID):
    res = requests.get("http://hdock.phys.hust.edu.cn/data/%s/" % (HDOCK_jobID, ))
    res.raise_for_status()
    soup = BeautifulSoup(res.text, 'html.parser')
    match = soup.find_all(string=re.compile("Your job will be finished in about .* seconds.", re.DOTALL))
    if len(match) == 0:
        return -1
    
    match = match[0].replace("\n", "").encode("ascii")
    seconds_left = float(match.split()[-2])
    return seconds_left

def get_HDOCK_status(HDOCK_jobID, test_exists=True):
    """
    Get HDOCK running status
    
    Return: NOT_EXISTS COMPLETE RUNNING
    """
    if test_exists:
        isExists = is_HDOCK_exists(HDOCK_jobID)
        if not isExists:
            return "NOT_EXISTS"
    
    isComplete = is_HDOCK_complete(HDOCK_jobID)
    if isComplete:
        return "COMPLETE"
    
    return "RUNNING"

def fetch_HDOCK_results(HDOCK_jobID, outTarGZFn, verbose=False):
    """
    Download all PDB files and energy file
    
    out_folder              -- A directory to contain all files
    outTarGZFn              -- Out filename (.tar.gz)
    verbose                 -- Output some information
    
    """
    
    url = "http://hdock.phys.hust.edu.cn/data/%s/all_results.tar.gz"%(HDOCK_jobID,)
    if verbose:
        os.system("wget --no-verbose "+url+" -O "+outTarGZFn)
    else:
        os.system("wget --quiet "+url+" -O "+outTarGZFn)

def fetch_HDOCK_top10_results(HDOCK_jobID, outTarGZFn, verbose=False):
    """
    Download top 10 PDB files
    
    out_folder              -- A directory to contain all files
    outTarGZFn              -- Out filename (.tar.gz)
    verbose                 -- Output some information
    
    """
    url = "http://hdock.phys.hust.edu.cn/data/%s/top10_models.tar.gz"%(HDOCK_jobID,)
    if verbose:
        os.system("wget --no-verbose "+url+" -O "+outTarGZFn)
    else:
        os.system("wget --quiet "+url+" -O "+outTarGZFn)

