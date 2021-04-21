#!/usr/bin/env python
import math, os, shutil, itertools, stk, time, concurrent.futures, string, os.path
import subprocess as sp
import pandas as pd
import numpy as np
import argparse

from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops, SanitizeMol, rdDistGeom, AllChem
from time import sleep
from utils import *

class polyscreen:
    def __init__(self, Id, monomers, style, length, solvent_params, threads):
        self.pwd = os.getcwd()
        self.Id = Id
        self.monomers = monomers
        self.style = style
        self.length = length
        self.solvent_params = solvent_params
        self.threads = threads

    def getPolymerWithConformer(self):
        """
        Using STK, the polymer is created from a combination of n monomer smiles strings.
        At the moment, bromine groups are used as to indicate where the monomers link, but
        boronic acid group may be more appropriate as the former means we cannot include
        monomers that contain genuine bromo groups.
        """
        bbs = []
        
        # Making the building blocks and hence the polymer from the monomer smiles strings.
        for m in self.monomers:
            bbs.append(stk.BuildingBlock(smiles=m, functional_groups=[stk.BromoFactory()]))
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=bbs,
                repeating_unit=self.style,
                num_repeating_units=self.length,
            ),
        )
        
        # Replacing the unused terminal bromines with H
        H = Chem.MolFromSmarts("[H]")
        Br = Chem.MolFromSmarts("[Br]")
        molObj = polymer.with_canonical_atom_ordering().to_rdkit_mol()
        molObj = Chem.rdmolops.ReplaceSubstructs(molObj, Br, H, replaceAll = True)[0]
        molObj = Chem.rdmolops.ReplaceSubstructs(molObj, Br, H)[0]
        molObj = rdmolops.RemoveHs(molObj)

        rdmolops.SanitizeMol(molObj)
        self.smiles = Chem.MolToSmiles(molObj)
        molObj = Chem.AddHs(molObj)
        
        # Setting parameters for ETKDG for the conformer search.
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        params.randomSeed = -1
        params.numThreads = self.threads
        params.useRandomCoords = True

        confNum = 500
        cids = Chem.rdDistGeom.EmbedMultipleConfs(molObj, confNum, params=params)

        # Minimizing the conformer set, and hence saving the geometry of the lowest energy conformer
        res = AllChem.MMFFOptimizeMoleculeConfs(molObj, numThreads = self.threads)
        mmffenergies = []
        for r in res:
            mmffenergies.append(r[1])
        lowest = mmffenergies.index(min(mmffenergies))
        Chem.rdmolfiles.MolToXYZFile(molObj, "{}.xyz".format(str(self.Id)), confId=cids[lowest])

    def runCalculations(self):
        """
        In order for VIP/VEA calculations to run, we have to use xTB 5.6.4SE and as such, must also
        use FakeTime as old versions of xTB appear to be timelocked. The properties of each calculation
        search for a specific string in the output stream, and finding the output values relative to
        the lookup string.
        """
        spawn_xtb = ["faketime", "2018-1-1 00:00:00", "xtb"]
        xyzfile = f"{str(self.Id)}.xyz"

        # xTB OPTIMIZATION

        command = spawn_xtb + [xyzfile, "-opt", self.solvent_params]
        output = run(command)

        with open(f"opt.out", "x+") as f:
            f.write(output)

        lines = output.splitlines()
        E_lines = []
        solv = []
        for idx, line in enumerate(lines):
            if "total E" in line:
                E_lines.append(line)
            if "Etot ala COSMO" in line:
                solv.append(line)
        E_line = E_lines[-1].split()
        solv_line = solv[-1].split()
        if len(self.solvent_params) > 0:
            E_xtb, E_solv = float(E_line[-1]), float(solv_line[-1])
        else:
            E_xtb, E_solv = float(E_line[-1]), None

        shutil.copy("xtbopt.xyz", f"{str(self.Id)}-opt.xyz")
        xyzfile = f"{self.Id}-opt.xyz"

        # xTB VIP CALCULATION
        command = spawn_xtb + [xyzfile, "-vip", self.solvent_params]
        output = run(command)
        with open("vip.out", "x+") as f:
            f.write(output)
        vip = output[output.find("delta SCC IP") :].split()[4]

        # xTB VEA CALCULATION
        command = spawn_xtb + [xyzfile, "-vea", self.solvent_params]
        output = run(command)
        with open("vea.out", "x+") as f:
            f.write(output)
        vea = output[output.find("delta SCC EA") :].split()[4]

        # xTB WAVEFUNCTION CALCULATION (xtb4stda)
        command = ["xtb4stda", xyzfile, self.solvent_params]
        output = run(command)

        # xTB EXCITATIONS CALCULATIONS
        command = ["stda", "-xtb", "-e", "8"]
        output = run(command)
        with open("stda-calc.out", "w") as f:
            f.write(output)
        lines = output.splitlines()
        for idx, line in enumerate(lines):
            if "excitation energies, transition moments" in line:
                gap = lines[idx + 2].split()[1]
                fL = lines[idx + 2].split()[3]

        # Processing results
        results = [self.Id, self.smiles, E_xtb, E_solv, vip, vea, gap, fL]
        cols = ['Id', 'smi', 'E_xtb', 'E_solv', 'vip', 'vea', 'gap', 'fL']
        removeJunk()
        
        # Adding A/B/C etc titles to columns so we can have monomer smiles in an individual column
        for idx, monomer in enumerate(self.monomers):
            name = list(string.ascii_uppercase)[idx]
            results.insert(idx+1, monomer)
            cols.insert(idx+1, name)

        df = pd.DataFrame(results)
        df = df.T
        df.columns = cols
        df.to_csv("output.csv", index=False)
        return results, cols

    ### END OF CLASS ###

def getProps(Id, monomers, style, length, solvent_params, threads):
    """
    This function allows for the parallelisation of the screening protocol using the
    screening class, defined above. Each polymer is given its own directory wherein the
    ETKDG and xTB geometries are saved, as well as the logs of the xTB calculations.
    """
    pwd = os.getcwd()
    os.mkdir(str(Id))
    os.chdir(str(Id))

    ps = polyscreen(Id, monomers, style, length, solvent_params, threads)
    ps.getPolymerWithConformer()
    results, cols = ps.runCalculations()
    os.chdir(pwd)
    return [results, cols]

def runScreen(name, monomer_list, style, length, solvent):
    """
    The screening is initialised, using concurrent.futures from within this function. Presently,
    the number of concurrent jobs can be changed with 'in_parallel', but the total cores and
    OMP_NUM_THREADS must be changed accordingly. OMP_NUM_THREADS is used as that is how xTB
    reads how many threads are to be used.
    """
    # Check solvent validity
    if solvent != None:
        if solvent not in valid_solvents:
            raise Exception(
                "Invalid solvent choice. Valid solvents:",
                [i for i in valid_solvents],
            )
        else:
            solvent_params = f"-gbsa {solvent}"
    else:
        solvent_params = ""

    # Generate all possible combinations of monomers from the provided list.
    n = len(set([char for char in style]))
    combinations = list(itertools.combinations_with_replacement(monomer_list, n))
    threads = int(os.environ['OMP_NUM_THREADS'])

    # Combining the arguments for getProps into iterables to be cycled with concurrent.futures.
    args = []
    for idx, combination in enumerate(combinations):
        l = [idx, combination, style, length, solvent_params, threads]
        args.append(l)
    in_parallel = 2
    chunksize = int(math.ceil(len(combinations) / in_parallel))
    pwd = os.getcwd()
    os.mkdir(name)
    os.chdir(name)

    # Launching the parallelised screening protocol
    with concurrent.futures.ProcessPoolExecutor(max_workers=in_parallel) as executor:
        x = executor.map(getProps, *zip(*args), chunksize=chunksize)
    lst = list(x)
    transpose = list(map(list, zip(*lst)))
    data = transpose[0]
    cols = transpose[1][0]

    props_df = pd.DataFrame(data)
    print(props_df)
    props_df.columns = cols
    props_df.to_csv("output.csv", index=False)
    os.chdir(pwd)
#############################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='filename', default='monomers.txt')
targs = parser.parse_args()

ms = [
    'BrC1=CC2=CC3=C(C=C(Br)N3)C=C2N1',
    'BrC1=CC2=CC3=C(C=C(Br)O3)C=C2O1',
    'BrC1=CC=C(Br)C2=NON=C12',
    'Brc3ccc(c2ccc(c1ccc(Br)cc1)cc2)cc3',
]

if (os.path.isfile(targs.filename)):
    with open(targs.filename, 'r') as f:
        monomers = f.readlines()
    monomers = [x.strip() for x in monomers]
else:
    monomers = ms

runScreen('6n', ms, 'AB', 6, 'h2o')

