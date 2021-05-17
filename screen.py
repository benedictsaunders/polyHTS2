from protocol import *
import os, argparse

# runScreen(name, monomer_list, style, length, solvent):
homedir = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='+', dest='filename', default='monomers.txt')
parser.add_argument('-l', dest='length', default=4)
parser.add_argument('-r', dest='style', default='AB')
parser.add_argument('-s', dest='solvent', default='h2o')
parser.add_argument('-p', dest='parallel', default=4)
parser.add_argument('-n', dest='name', default='screen')
parser.add_argument('-d', dest='database', default='database.db')
parser.add_argument('--exhaustive', action="store_true")
parser.add_argument('--ForceNew', action="store_true")
targs = parser.parse_args()

monomer_list = []

for fn in targs.filename:
    if (os.path.isfile(fn)):
        with open(fn, 'r') as f:
            monomers = f.readlines()
        monomers = [x.strip() for x in monomers]
        monomer_list.append(monomers)
    else:
        raise(f'File {fn} does not exist.')

data, cols = runScreen(targs.name, monomer_list, targs.style, targs.length, targs.solvent, targs.parallel, targs.database, targs.exhaustive, targs.ForceNew)
