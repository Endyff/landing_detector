import argparse
from pymol import cmd
from pathlib import Path
from tqdm import tqdm

import pandas as pd
import numpy as np
import os
import logging

POCKET_KEYWORD='pocket'
LIGAND_KEYWORD='ligand'
COMPLEX_KEYWORD='complex'

# SET UP CUSTOM LIGAND_PROTEIN PATH WHEN USING --data_dir
# SET UP FOR PLINDER DATASET
def protein_path(system_path):
    return system_path / 'receptor.pdb'

def ligand_path(system_path):
    return list((system_path / 'ligand_files').iterdir())[0]

def get_paths(system_path):
    systems, protein_paths, ligand_paths = [],[],[]
    for path in Path(system_path).iterdir():
        systems.append(path.name)
        protein_paths.append(protein_path(path))
        ligand_paths.append(ligand_path(path))
    return pd.DataFrame({'system_id': systems, 'protein_path': protein_paths, 'ligand_path': ligand_paths})



def setup_parser():
    parser = argparse.ArgumentParser(description='Buried ratio calculator')

    # input output
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--data_dir', type=str, help='Directory containing system subdirectories')
    group.add_argument('--path_csv', type=str, help='CSV file containing paths to systems')
    parser.add_argument('--result_csv', type=str, default='buried_ratio.csv', help='Path to the result CSV file')
    
    # log
    parser.add_argument('--log_file', type=str, default='buried_ratio_annotation.log', help='Path to the log file')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logging level')

    # pymol_settings
    parser.add_argument('--dot_density', type=int, default=2, help='Density of the dots around the atoms')
    parser.add_argument('--dot_solvent', type=int, default=1, help='Boolean switch between using vdW radii computation of the surface (0) or the solvent accessible area (1)')
    
    return parser

def setup_logging(log_filename, logging_level):
    logging.basicConfig(level=logging_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', filename=log_filename)
    return logging.getLogger(__name__)


def get_buried_ratio(protein_path, ligand_path, dot_solvent, dot_density):
    # method by Petr Kouba 
    # https://github.com/KoubaPetr/binding-classifier

    cmd.reinitialize()
    cmd.load(str(protein_path), POCKET_KEYWORD)
    cmd.load(str(ligand_path), LIGAND_KEYWORD)
    cmd.h_add()
    cmd.create(COMPLEX_KEYWORD, f'{LIGAND_KEYWORD} or {POCKET_KEYWORD}')

    cmd.flag(flag='ignore', selection='all', action='clear')
    cmd.flag(flag='ignore', selection='solvent', action='set') 

    cmd.set('dot_solvent', dot_solvent) # this is boolean switch between using vdW radii computation of the surface (0) or the solvent accessible area (1)
    cmd.set('dot_density', dot_density) # this takes integers 0 to 4 - to set the density of the dots around the atoms

    # calculate the surface of the lipid and the pocket
    lipid_area = cmd.get_area(selection=LIGAND_KEYWORD)
    protein_area = cmd.get_area(selection=POCKET_KEYWORD)
    complex_area = cmd.get_area(selection=COMPLEX_KEYWORD)
    log.debug(f'areas: \n\t{lipid_area=}\n\t{protein_area=}\n\t{complex_area=}')

    buried_ratio = (lipid_area + protein_area - complex_area) / (2*lipid_area) if lipid_area != 0 else np.nan
    control_ratio = complex_area / protein_area if protein_area != 0 else np.nan
    if buried_ratio == np.nan:
        log.error(f"In case of {protein_path} {ligand_path} the lipid area is 0 !")
    if control_ratio == np.nan:
        log.error(f"In case of {protein_path} {ligand_path} the protein area is 0 !")
    return buried_ratio, control_ratio


def process_complexes(args, log):

    processed_systems = set(pd.read_csv(args.result_csv).system_id) if Path(args.result_csv).exists()  else {}
    if args.data_dir:
        paths_df = get_paths(args.data_dir)
    elif args.path_csv:
        paths_df = pd.read_csv(args.path_csv)
    else:
        print('No data source provided.')
        log.error("No data source provided.")
        return
    log.info(f"Processing {len(paths_df)} systems")
    for _, row in paths_df.iterrows():
        system_id = row.system_id
        protein_path = Path(row.protein_path)
        ligand_path = Path(row.ligand_path)
        if system_id in processed_systems and not args.overwrite:
            log.debug(f"Skipping {system_id}")
            continue
        if protein_path.exists() and ligand_path.exists():
            try:
                buried_lipid_ratio, control_ratio = get_buried_ratio(protein_path, ligand_path, args.dot_solvent, args.dot_density)
                result_dict = {
                    'system_id': system_id,
                    'buried_lipid_ratio': buried_lipid_ratio,
                    'control_ratio': control_ratio
                }
                pd.DataFrame([result_dict]).to_csv(args.result_csv, index=False, mode='a', header=not os.path.exists(args.result_csv))

            except Exception as e:
                log.warning(f"Error processing {system_id}: {e}")
        else:
            log.warning(f"Files for {system_id} not found.")
        
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    log = setup_logging(args.log_file, args.log_level)
    process_complexes(args, log)