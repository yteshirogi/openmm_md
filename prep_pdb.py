from pathlib import Path

from rdkit import Chem
import requests


def get_pdb_each(pdbid, ligand_name, DATA=Path("./data")):
    """
    Download pdb file and ligand sdf file from rcsb.org
    """
    DATA.mkdir(exist_ok=True)
    pdb_path = DATA / f"{pdbid}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    ligand_path = DATA / f"{ligand_name}.sdf"
    ligand_url = f"https://files.rcsb.org/ligands/download/{ligand_name}_ideal.sdf"
    for path, url in zip([pdb_path, ligand_path], [pdb_url, ligand_url]):
        r = requests.get(url)
        r.raise_for_status()
        with open(path, "wb") as f:
            f.write(r.content)

    lig_smiles = converter(str(ligand_path))

    return pdb_path, lig_smiles

def get_pdb(pdbid, DATA=Path("./data")):
    """
    Download pdb file from rcsb.org
    """
    DATA.mkdir(exist_ok=True)
    pdb_path = DATA / f"{pdbid}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    r = requests.get(pdb_url)
    r.raise_for_status()
    with open(pdb_path, "wb") as f:
        f.write(r.content)

    return pdb_path


def converter(sdf_file, DATA=Path("./data")):
    supplier = Chem.SDMolSupplier(sdf_file)
    
    smiles_file = DATA / "smiles.txt"
    with open(smiles_file, "w") as f:
        for mol in supplier:
            if mol is not None:             # avoiding compounds that cannot be loaded.
                smiles = Chem.MolToSmiles(mol)
                f.write(f"{smiles}\n")
    
    return smiles


# get current directory when called
# HERE = Path(__file__).parent.resolve()
# create data directory if not exists
# DATA = HERE / "data"
# DATA.mkdir(exist_ok=True)
