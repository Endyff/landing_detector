import argparse
import os
import shutil
from pathlib import Path

from plinder.core import PlinderSystem, get_config
from plinder.core.scores import query_index
import logging




parser = argparse.ArgumentParser(description='Download complexes from plinder')
parser.add_argument('--filter', type=str, default=None, help='Filtering list')
parser.add_argument('--system_id', type=str, default=None, help='System ID list')
parser.add_argument('--remove_others_once', action='store_true', help='Remove other files')
parser.add_argument('--remove_other_every', action='store_true', help='Remove other files')
parser.add_argument('--num_complexes', type=int, default=3, help='Number of complexes to download')
parser.add_argument('--ligand_chains', type=int, default=1, help='Number of ligand chains')
parser.add_argument('--protein_chains', type=int, default=1, help='Number of protein chains')
parser.add_argument('--resolution', type=float, default=3.0, help='Resolution of the complex')
# parser.add_argument('--split', type=str, default=None, help='Split of the complex')
parser.add_argument('--lipinski', action='store_true', help='Apply Lipinski filter')

parser.add_argument('--filter_file', type=str, default='filters.txt', help='Filter file')
parser.add_argument('--system_list_file', type=str, default=None, help='System list file')
parser.add_argument('--plinder_dir', type=str, default='/Users/bouceond/phd/landing_detector/data', help='Plinder directory')




args = parser.parse_args()


def download_complex(system_id):
    plinder_system = PlinderSystem(system_id=system_id).ligands
    plinder_system.ligands

def clear_other_files(system_id_list, data_dir):

    for path in Path(f'{str(data_dir)}/systems').iterdir():
        # remove unwanted complexes
        if os.path.isdir(str(path)) and path.name not in system_id_list:
            shutil.rmtree(path)
        # zip bucked that contained the system
        if path.suffix == '.zip':
            os.remove(path)
        # remove check that bucket was processed
        if path.name[-4:] == 'done':
            os.remove(path)


def parse_filter_file(filter_file):
    filters = []
    with open(filter_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            filters.append(tuple(line.split(' ')))
    return filters

def get_filters(args):
    filters = []
    if args.ligand_chains is not None:
        filters.append(("system_num_ligand_chains", "<=", args.ligand_chains))
    if args.protein_chains is not None:
        filters.append(("system_num_protein_chains", "<=", args.protein_chains))
    if args.resolution is not None:
        filters.append(("system_resolution", "<=", args.resolution))
    if args.lipinski:
        filters.append(("ligand_is_lipinski", "==", True))

    return filters

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', filename='plinder_downloader.log')
    log = logging.getLogger(__name__)

    agrs = parser.parse_args()
    os.makedirs(args.plinder_dir, exist_ok=True)
    os.environ["PLINDER_MOUNT"] = args.plinder_dir
    
    config = get_config()
    plinder_dir = config.data.plinder_dir
    log.info(f"Plinder directory: {plinder_dir}")

    if args.system_list_file is not None:
        log.info(f"Reading system list from {args.system_list_file}")
        with open(args.system_list_file, 'r') as f:
            system_id_list = f.readlines()
    else:
        log.info("Querying system list")
        if args.filter_file:
            filters = parse_filter_file(args.filter_file)
        else:
            filters = get_filters(args)
        print(filters)
        filters_str = "\n\t".join(map(repr, filters))
        log.info(f"Filters: \n\t{filters_str}")
        # add split
        df = query_index(filters=filters)
        system_id_list = df['system_id'].tolist()
        log.info(f"Found {len(system_id_list)} systems")

    system_id_complexes = system_id_list[:args.num_complexes] if args.num_complexes else system_id_list
    print(system_id_complexes)
    log.info(f"Downloading {len(system_id_complexes)} complexes")
    for system_id in system_id_complexes:

        log.debug(f"Downloading {system_id}")
        PlinderSystem(system_id=system_id).ligands

        if args.remove_other_every:
            clear_other_files(system_id_complexes, plinder_dir)

    if args.remove_others_once:
        clear_other_files(system_id_complexes, plinder_dir)

    log.info("Download complete")

if __name__ == '__main__':
    main()