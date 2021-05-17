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
parser.add_argument('-d', dest='database', default='database.db')
parser.add_argument('--ForceNew', action="store_true")
targs = parser.parse_args()

if (os.path.isfile(targs.filename)):
    with open(targs.filename, 'r') as f:
        monomers = f.readlines()
    monomers = [x.strip() for x in monomers]
else:
    monomers = ms

data, cols = runScreen(targs.name, monomers, targs.style, targs.length, 'h2o', targs.parallel, targs.database)
