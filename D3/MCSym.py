
"""

This module can predict RNA 3D-structure with MCSym

from D3 import MCSym

## 1. Upload a task
seq = "UGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUACC"
dot = "((((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))).."
jobID, res = MCSym.upload_MCSym_job(seq, dot, seq_title="let-7", time_limit="48h", gen_modellimit=10, verbose=True)

## 2. Wait a task -- very long time
status = MCSym.get_MCSym_status(jobID)
while status == 'RUNNING':
    status = MCSym.get_MCSym_status(jobID)
    print(status)
    time.sleep(2)

## 3. Minimize a task and wait -- very long time
thread_handle = MCSym.minimize_MCSym_newThread(jobID, log_out="/tmp/logout")
# You can see /tmp/logout about log information
thread_handle.join()

## 4. Score and sort -- short time
thread_handle = MCSym.score_MCSym_newThread(jobID, log_out="/tmp/logout")
# You can see /tmp/logout about log information
thread_handle.join()

## 5. Download all PDB
thread_handle = MCSym.fetch_top_MCSym_pdb_newThread(jobID, "/tmp/outFolder", topnum=50)
thread_handle.join()

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

def check_valid_structure(dot):
    if 1.0*dot.count('.')/len(dot) > 0.39:
        return False
    return True

def upload_MCSym_job(seq, dot, seq_title, time_limit="48h", gen_modellimit=1000, verbose=True, proxies=None):
    """
    seq                 -- RNA sequence
    dot                 -- Dotbracket structure
    seq_title           -- Sequence title
    time_limit          -- Maximun time to run: 30m, 1h, 12h, 24h, 48h, none
    gen_modellimit      -- How many models to produce
    verbose             -- Output some information
    proxies             -- Proxy
                            such as: {'http': 'http://username:password@ip:port', 'https': 'http://usernamepassword@ip:port'}
    
    Return:
        jobID, res
        jobID           -- A ID string
        res             -- Requests Responds object
    
    Exception:
        RuntimeError
    """
    assert time_limit in ('30m','1h','12h','24h','48h','none')
    seq = seq.upper()
    assert len(seq) == len(dot)
    
    if not check_valid_structure(dot):
        raise RuntimeError("Error: The number of unpaired bases should not exceed 40%")
    
    input_structure = ">%s\n%s\n%s" % (seq_title, seq, dot)
    param = {   
        "scriptgen":input_structure,
        "action":"Submit",
        "gen_fragRMSD": "0.1",
        "gen_theo": "1",
        "gen_mergeRMSD": "1.5",
        "gen_clash": "1.5",
        "gen_ribosethrs": "2.0",
        "gen_ribosemthd": "ccm",
        "gen_method": "probabilistic",
        "gen_modellimit": str(gen_modellimit),
        "gen_timelimit": time_limit,
        "gen_modeldvsty": "1.0"
    }
    
    if verbose:
        print(input_structure+"\n")
    
    headers = {'User-Agent': 'Mozilla/5.0'}
    if proxies:
        res = requests.post("https://www.major.iric.ca/cgi-bin/MC-Sym/mcsym.cgi", data=param, headers=headers, proxies=proxies)
    else:
        res = requests.post("https://www.major.iric.ca/cgi-bin/MC-Sym/mcsym.cgi", data=param, headers=headers)
    res.raise_for_status()
    jobID = re.findall("<a HREF=\"http://www.major.iric.ca/MC-Sym/Work/(\\w+)/\" target=\"_blank\">HERE</a>", res.text)
    if len(jobID) == 0:
        if verbose:
            print("Failed: maybe you have exceeded the maximum limit")
            print("Hint: print(res.text) to see")
        return -1, res
    jobID = jobID[0]
    if verbose:
        print("JobID: "+str(jobID))
        print("https://major.iric.ca/MC-Sym/Work/"+str(jobID)+"/")
    
    return jobID, res

def is_MCSym_exists(MCSym_jobID):
    try:
        mcsym_out = urllib2.urlopen("https://major.iric.ca/MC-Sym/Work/%s/mcsym.out" % (MCSym_jobID, ))
    except urllib2.HTTPError:
        return False
    return True

def is_MCSym_error(MCSym_jobID):
    mcsym_err = urllib2.urlopen("https://major.iric.ca/MC-Sym/Work/%s/mcsym.err" % (MCSym_jobID, ))
    lines = mcsym_err.readlines()
    mcsym_err.close()
    if len(lines)==0:
        return False
    if 'error' in lines[-1]:
        return True
    else:
        return False

def is_MCSym_complete(MCSym_jobID):
    mcsym_out = urllib2.urlopen("https://major.iric.ca/MC-Sym/Work/%s/mcsym.out" % (MCSym_jobID, ))
    lines = mcsym_out.readlines()
    mcsym_out.close()
    if len(lines)==0:
        return False
    if "report for structure saved in log file" in lines[-1]:
        return True
    else:
        return False

def is_MCSym_minimizing(MCSym_jobID):
    url = "https://major.iric.ca/MC-Sym/Work/%s/"%(MCSym_jobID,)
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    for node in soup.find_all('a'):
        if node.get('href').endswith('.xyz'):
            return True
    return False

def get_MCSym_status(MCSym_jobID, test_exists=True):
    """
    Get MCSym running status
    
    Return: NOT_EXISTS COMPLETE ERROR RUNNING
    """
    if test_exists:
        isExists = is_MCSym_exists(MCSym_jobID)
        if not isExists:
            return "NOT_EXISTS"
    
    isComplete = is_MCSym_complete(MCSym_jobID)
    if isComplete:
        return "COMPLETE"
    
    isError = is_MCSym_error(MCSym_jobID)
    if isError:
        return "ERROR"
        
    return "RUNNING"

def minimize_MCSym(MCSym_jobID, minitype="refine", log_out="/dev/null"):
    """
    Minimize raw models
    
    minitype                -- Method to minimize: relieve, refine, brushup, anneal
    log_out                 -- Output information
    """
    assert minitype in ('relieve', 'refine', 'brushup', 'anneal')
    param = { "workkey": MCSym_jobID, "minimize": minitype }
    headers = {'User-Agent': 'Mozilla/5.0'}
    res = requests.post("https://www.major.iric.ca/cgi-bin/MC-Sym/mcsym.cgi", data=param, headers=headers, stream=True)
    fd = open(log_out, 'wb')
    for chunk in res.iter_content(chunk_size=128):
        fd.write(chunk)
        fd.flush()
    fd.close()

def minimize_MCSym_newThread(MCSym_jobID, minitype="refine", log_out="/dev/null"):
    """
    Minimize raw models with a new thread
    
    minitype                -- Method to minimize: relieve, refine, brushup, anneal
    log_out                 -- Output information
    
    Return Thread object
    """
    thd = threading.Thread(target=minimize_MCSym, args=(MCSym_jobID,minitype,log_out))
    thd.start() 
    return thd

def score_MCSym(MCSym_jobID, log_out="/dev/null"):
    """
    Score models
    
    log_out                 -- Output information
    """
    param = { "workkey": MCSym_jobID, "score": '1' }
    headers = {'User-Agent': 'Mozilla/5.0'}
    res = requests.post("https://www.major.iric.ca/cgi-bin/MC-Sym/mcsym.cgi", data=param, headers=headers, stream=True)
    fd = open(log_out, 'wb')
    for chunk in res.iter_content(chunk_size=128):
        fd.write(chunk)
        fd.flush()
    fd.close()

def score_MCSym_newThread(MCSym_jobID, log_out="/dev/null"):
    """
    Score models with a new thread
    
    log_out                 -- Output information
    
    Return Thread object
    """
    thd = threading.Thread(target=score_MCSym, args=(MCSym_jobID,log_out))
    thd.start() 
    return thd

def read_MCSym_pdb_energy(MCSym_jobID):
    """
    Read MCSym energy file: energy.sorted.data. score_MCSym must be runned
    
    Return [ (pdb_file_name, energy), (pdb_file_name, energy), (pdb_file_name, energy)... ]
    """
    energy_list = []
    mcsym_energy = urllib2.urlopen("https://major.iric.ca/MC-Sym/Work/%s/energy.sorted.data" % (MCSym_jobID, ))
    for line in mcsym_energy.readlines():
        fn,energy = re.findall("Total energy of `(\\S+)` is:\\s*(\\S+) Kcal/mol.", line)[0]
        energy = float(energy)
        energy_list.append( (fn, energy) )
    return energy_list

def fetch_MCSym_pdb(MCSym_jobID, out_folder, energy_fn=None, verbose=False):
    """
    Download all PDB files and energy file
    
    out_folder              -- A directory to contain all files
    energy_fn               -- Energy file name
    verbose                 -- Output some information
    
    """
    out_folder = out_folder.rstrip('/') + '/'
    
    pdb_energy = read_MCSym_pdb_energy(MCSym_jobID)
    if energy_fn:
        OUT = open(out_folder+energy_fn, "w")
    else:
        OUT = open(out_folder+"energy.txt", "w")
    OUT.writelines( "\n".join([ it[0]+"\t"+str(it[1]) for it in pdb_energy ]) )
    OUT.close()
    
    url = "https://major.iric.ca/MC-Sym/Work/%s/"%(MCSym_jobID,)
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    for node in soup.find_all('a'):
        if node.get('href').endswith('.pdb.gz'):
            file = node.get('href')
            if verbose:
                os.system("wget --no-verbose "+url+file+" -O "+out_folder+file)
            else:
                os.system("wget --quiet "+url+file+" -O "+out_folder+file)

def fetch_top_MCSym_pdb(MCSym_jobID, out_folder, topnum=50, energy_fn=None, verbose=False):
    """
    Download top minimum energy PDB files and energy file
    
    out_folder              -- A directory to contain all files
    topnum                  -- Number of structures with minimun energy
    energy_fn               -- Energy file name
    verbose                 -- Output some information
    
    """
    out_folder = out_folder.rstrip('/') + '/'
    
    pdb_energy = read_MCSym_pdb_energy(MCSym_jobID)
    if energy_fn:
        OUT = open(out_folder+energy_fn, "w")
    else:
        OUT = open(out_folder+"energy.txt", "w")
    OUT.writelines( "\n".join([ it[0]+"\t"+str(it[1]) for it in pdb_energy ]) )
    OUT.close()
    
    url = "https://major.iric.ca/MC-Sym/Work/%s/"%(MCSym_jobID,)
    #page = requests.get(url).text
    #soup = BeautifulSoup(page, 'html.parser')
    for raw_file,energy in pdb_energy[:topnum]:
        file = raw_file + ".gz"
        if verbose:
            os.system("wget --no-verbose "+url+file+" -O "+out_folder+file)
        else:
            os.system("wget --quiet "+url+file+" -O "+out_folder+file)

def fetch_top_MCSym_pdb_newThread(MCSym_jobID, out_folder, topnum=50, energy_fn=None):
    """
    Download top minimum energy PDB files and energy file with a new thread
    
    out_folder              -- A directory to contain all files
    topnum                  -- Number of structures with minimun energy
    energy_fn               -- Energy file name
    verbose                 -- Output some information
    
    Return Thread object
    """
    thd = threading.Thread(target=fetch_top_MCSym_pdb, args=(MCSym_jobID, out_folder, topnum, energy_fn))
    thd.start()
    return thd














