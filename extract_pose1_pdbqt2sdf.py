# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:24:18 2025

@author: ismail
"""
import os
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd

def extract_affinity(pdbqt_file):
    with open(pdbqt_file, 'r') as file:
        lines = file.readlines()
    
    pose_lines = []
    in_pose = False
    affinity = None
    
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            affinity = float(line.split()[3])
        if line.startswith("MODEL 1"):
            in_pose = True
        if in_pose:
            pose_lines.append(line)
        if line.startswith("TORSDOF"):
            pose_lines.append("ENDMDL")
            break
    return affinity

#print(extract_affinity('204_out.pdbqt'))

def pdbqt2sdf(pdbqt_file):
    # Extract binding affinities from the PDBQT file
    affinities = extract_affinity(pdbqt_file)

    # Load the PDBQT file
    molecules = pybel.readfile("pdbqt", pdbqt_file)

    # Write to SDF file with binding affinity and model number
    output = pybel.Outputfile("sdf", pdbqt_file.split('.')[0]+'_pose1'+'.sdf', overwrite=True)
    for idx, mol in enumerate(molecules):
        # Append model number to the molecule title
        mol.title = pdbqt_file.split('.')[0]
        # Add binding affinity as a property
        mol.data["BindingAffinity"] = affinities
        output.write(mol)
        if idx == 0 :
            break
    output.close()
    print("Conversion complete: sdf_file created with binding affinities and model numbers.")

def concatenate_files():
    input_files = os.listdir(os.getcwd())
    output = pybel.Outputfile("sdf", 'docked_ligands_concatenated.sdf', overwrite=True)
    for file in input_files:
        if file.endswith('.sdf'):
            molecules = pybel.readfile("sdf", file)
            for mol in molecules:
                output.write(mol)
    output.close()
    return output

def rearrange_sdf(sdf_file):

    supplier = Chem.SDMolSupplier(sdf_file)

    # Convert to a Pandas DataFrame
    mols = [mol for mol in supplier if mol is not None]
    df = PandasTools.LoadSDF(sdf_file)

    # Specify the column to sort by (e.g., 'MolecularWeight')
    sort_column = 'BindingAffinity'

    # Ensure the column exists
    if sort_column not in df.columns:
        raise ValueError(f"Column '{sort_column}' not found in SDF file.")

    # Sort the DataFrame by the specified column
    df_sorted = df.sort_values(by=sort_column)

    # Write the sorted molecules to a new SDF file
    output_sdf = 'docked_ligands_sorted.sdf'
    writer = Chem.SDWriter(output_sdf)

    for mol in df_sorted['ROMol']:
    #    mol.SetProp("BindingAffinity", str(mol.GetProp("BindingAffinity")))
        writer.write(mol)

    writer.close()
    print(f"Sorted molecules have been written to {output_sdf}")


def main():
    for pdbqt_file in os.listdir(os.getcwd()):
        pdbqt2sdf(pdbqt_file)
    concatenate_files()
    #rearrange_sdf('docked_ligands_concatenated.sdf')

main()