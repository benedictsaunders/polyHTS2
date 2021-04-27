from protocol import *
import os, argparse

# runScreen(name, monomer_list, style, length, solvent):
homedir = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='filename', default='monomers.txt')
parser.add_argument('-l', dest='length', default=4)
parser.add_argument('-r', dest='style', default='AB')
parser.add_argument('-p', dest='parallel', default=4)
parser.add_argument('-n', dest='name', default='screen')
targs = parser.parse_args()

ms = [
    'BrC1=CC2=CC3=C(C=C(Br)N3)C=C2N1',
    'BrC1=CC2=CC3=C(C=C(Br)O3)C=C2O1',
    'BrC1=CC=C(Br)C2=NON=C12',
    'C1(=C(C2=NSN=C2C(=C1F)Br)Br)F',
]

if (os.path.isfile(targs.filename)):
    with open(targs.filename, 'r') as f:
        monomers = f.readlines()
    monomers = [x.strip() for x in monomers]
else:
    monomers = ms

data, cols = runScreen(targs.name, monomers, targs.style, targs.length, 'h2o', targs.parallel)
