
import Bio
from Bio.PDB import PDBList
import tempfile, shutil, os, Colors
import json, General

def ssecFromPDB(pdb, clean=True, verbose=False):
    """
    pdb                     -- A PDB id or PDB file, for example, 5CD1, /tmp/5CD1.cif
                               If you provide a PDB ID, it will download pdf file automaticly
    clean                   -- Clean the directory
    verbose                 -- Print the command to the screen
    
    Return:
        { 'm1_chain_L5': {'bseq': 'ACTAGTCAGCTAGC', 'sstr': '....((((..))))...'}, ... }
    """
    import tempfile, shutil, os, Colors
    import json, General
    
    tmpdir = tempfile.mkdtemp(prefix="PDB_", dir='/tmp/')
    tmp_result_json = os.path.join(tmpdir, "result.json")
    
    if len(pdb)==4:
        import Bio
        from Bio.PDB import PDBList
        pdbl = PDBList(verbose=False)
        pdbFn = pdbl.retrieve_pdb_file(pdb, pdir=tmpdir, overwrite=True)
    else:
        pdbFn = pdb
    
    cmd = f'x3dna-dssr -i={pdbFn} -o={tmp_result_json} --json 2>/dev/null'
    if verbose:
        print( Colors.f(cmd, 'yellow') )
    os.system(cmd)
    
    x3dna_result = json.load(open(tmp_result_json))
    ssec = {}
    for key in x3dna_result['dbn']:
        if key != 'all_chains':
            ssec[key] = { 'bseq':x3dna_result['dbn'][key]['bseq'], 'sstr':x3dna_result['dbn'][key]['sstr'] }
    
    if clean:
        shutil.rmtree(tmpdir)
    
    return ssec


