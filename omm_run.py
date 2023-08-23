import fix_ligand
import fix_protein
import gen_ff
import gen_system
import merge_complex
import prep_pdb
import rd2omm
import run_md
import sys


def main(pdbid, ligand_name):
    """
    Run a simulation.
    """
    # download pdb file and ligand sdf file from rcsb.org
    # prep_pdb.get_pdb_each(pdbid, ligand_name)

    # download pdb file from rcsb.org
    pdb_path, lig_smiles = prep_pdb.get_pdb_each(pdbid, ligand_name)

    # fix protein
    protein = fix_protein.prepare_protein(pdb_path)

    # fix ligand
    ligand = fix_ligand.prepare_ligand(
        pdb_path,
        resname=ligand_name,
        smiles=lig_smiles,
        depict=False
    )

    # convert rdkit mol to openmm topology and positions
    omm_ligand= rd2omm.rdkit_to_openmm(ligand, name=ligand_name)
    # merge protein and ligand
    complex_topology, complex_positions = merge_complex.merge_protein_and_ligand(
        protein, omm_ligand
    )

    # generate forcefield
    forcefield = gen_ff.generate_forcefield(
        rdkit_mol=ligand,
        protein_ff="amber14-all.xml",
        solvent_ff="amber14/tip3pfb.xml"
    )

    # generate system
    system = gen_system.generate_system(
        complex_topology, forcefield, complex_positions
    )

    # run simulation
    run_md.run_simulation(system, steps=50000)


if __name__ == "__main__":
    # get pdbid and ligand name from command line
    pdbid = sys.argv[1]
    ligand_name = sys.argv[2]
    # smiles = "CC(C)(O)CC(=O)NCCn1ccc2ncnc(Nc3ccc(Oc4cccc(c4)C(F)(F)F)c(Cl)c3)c12"
    main(pdbid, ligand_name)
