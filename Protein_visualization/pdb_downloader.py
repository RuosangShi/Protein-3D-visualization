import os
import requests
import pandas as pd


def download_pdb_structure(uniprot_id, alphafold=False):

    """
    First checks pdb database, then alphafold db
    ______________
    Args:
        uniprot_id: Uniprot identifier of the protein
        alphafold: To only download from alphafold db
    ______________
    Returns:
        str: Path to downloaded structure file or None if no structure found
    """

    os.makedirs('pdb_structures', exist_ok=True)
    os.makedirs('alphafold_structures', exist_ok=True)

    if not alphafold:
        # First, check pdb 
        pdb_mapping_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
        try:
            response = requests.get(pdb_mapping_url)
            if response.status_code == 200:
                pdb_entries = []
                for line in response.text.split('\n'):
                    if line.startswith('DR   PDB;'):
                        pdb_id = line.split(';')[1].strip()
                        pdb_entries.append(pdb_id)
                
                if pdb_entries:
                    pdb_id = pdb_entries[0]
                    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    pdb_path = os.path.join('pdb_structures', f"{pdb_id}_{uniprot_id}.pdb")
                    
                    pdb_response = requests.get(pdb_url)
                    if pdb_response.status_code == 200:
                        with open(pdb_path, 'wb') as f:
                            f.write(pdb_response.content)
                        print(f"Downloaded PDB structure for {uniprot_id}: {pdb_path}")
                        return pdb_path
        except Exception as e:
            print(f"Error checking PDB for {uniprot_id}: {e}")

    # Try alphafold db
    alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    alphafold_path = os.path.join('alphafold_structures', f"{uniprot_id}_alphafold.pdb")
    
    try:
        response = requests.get(alphafold_url)
        if response.status_code == 200:
            with open(alphafold_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded AlphaFold structure for {uniprot_id}: {alphafold_path}")
            return alphafold_path
    except Exception as e:
        print(f"Error downloading from AlphaFold for {uniprot_id}: {e}")

    print(f"No structure found for {uniprot_id}")
    return None

def batch_download_structures(uniprot_ids):

    """
    Download structures for a list of uniprot ids
    ______________
    Args:
        uniprot_ids: List of uniprot ids
    ______________
    Returns:
        dict: Mapping of uniprot id to downloaded structure path
    """

    download_results = {}
    
    if isinstance(uniprot_ids, pd.Series):
        uniprot_ids = uniprot_ids.tolist()
    
    for uniprot_id in uniprot_ids:
        structure_path = download_pdb_structure(uniprot_id)
        download_results[uniprot_id] = structure_path
    
    return download_results

