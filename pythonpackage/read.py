import urllib3


def read_pdbid(pdbid):
    """
    Read a pdb file from the PDB. Need an internet connection.

    Parameters
    ----------
    pdbid: Name of the pdb file. (ex: 1dpx or 1dpx.pdb)

    Return
    ------
    The content of the pdb file
    """
    # URL of the PDB
    url = 'http://files.rcsb.org/download/'
    if pdbid[-4:] == '.pdb':
        file_name = pdbid
    else:
        file_name = pdbid + '.pdb'
    http = urllib3.PoolManager()
    url = url + file_name
    req = http.request('GET', url, preload_content = False)
    data = req.read()
    req.release_conn()
    return data
