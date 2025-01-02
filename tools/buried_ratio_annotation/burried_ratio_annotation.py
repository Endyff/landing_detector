import freesasa
import argparse

def get_buried_ratio(pdb_file, sdf_file):
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    total = result.totalArea()
    chains = result.chains()
    buried = 0
    for chain in chains:
        buried += chain.totalArea()
    return buried / total