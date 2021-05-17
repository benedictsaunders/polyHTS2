#!/usr/bin/env python
import math, os, shutil, itertools, stk, time, concurrent.futures, string, os.path
import subprocess as sp
import pandas as pd
import numpy as np
import argparse

from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops, SanitizeMol, rdDistGeom, AllChem
from time import sleep, time
from utils import *
import sqlite3
import sql
import constants

class polyscreen:
    def __init__(self, Id, monomers, style, length, solvent_params, threads):
        self.pwd = os.getcwd()
        self.Id = Id
        self.monomers = monomers
        self.style = style
        self.length = int(length)
        self.solvent_params = solvent_params
        self.threads = int(threads)

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
        repeat = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=bbs,
                repeating_unit=self.style,
                num_repeating_units=1,
            ),
        )

        # Replacing the unused terminal bromines with H
        H = Chem.MolFromSmarts("[H]")
        Br = Chem.MolFromSmarts("[Br]")

        repeat = repeat.with_canonical_atom_ordering().to_rdkit_mol()
        molObj = polymer.with_canonical_atom_ordering().to_rdkit_mol()

        molObj = Chem.rdmolops.ReplaceSubstructs(molObj, Br, H, replaceAll = True)[0]
        repeat = Chem.rdmolops.ReplaceSubstructs(repeat, Br, H, replaceAll = True)[0]

        molObj = rdmolops.RemoveHs(molObj)
        repeat = rdmolops.RemoveHs(repeat)

        rdmolops.SanitizeMol(molObj)
        rdmolops.SanitizeMol(repeat)

        self.smiles = Chem.MolToSmiles(molObj)
        self.repeat_smiles = Chem.MolToSmiles(repeat)

        molObj = Chem.AddHs(molObj)

        # Setting parameters for ETKDG for the conformer search.
        params = AllChem.ETKDGv3()

        params.useSmallRingTorsions = False
        params.ignoreSmoothingFailures = True
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
        sid = str(self.Id)
        spawn_xtb = ["faketime", "2018-1-1 00:00:00", "xtb"]
        xyzfile = f"{str(self.Id)}.xyz"

        # xTB OPTIMIZATION

        command = spawn_xtb + [xyzfile, "-opt", self.solvent_params]
        output = run(command)

        with open(f"opt.out", "w+") as f:
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
            E_xtb, E_solv = float(E_line[-1]), 0

        shutil.copy("xtbopt.xyz", f"{str(self.Id)}-opt.xyz")
        xyzfile = f"{self.Id}-opt.xyz"

        # xTB VIP CALCULATION
        command = spawn_xtb + [xyzfile, "-vip", self.solvent_params]
        output = run(command)
        with open("vip.out", "w+") as f:
            f.write(output)
        vip = float(output[output.find("delta SCC IP") :].split()[4])+constants.ELECTRODE

        # xTB VEA CALCULATION
        command = spawn_xtb + [xyzfile, "-vea", self.solvent_params]
        output = run(command)
        with open("vea.out", "w+") as f:
            f.write(output)
        vea = float(output[output.find("delta SCC EA") :].split()[4])+constants.ELECTRODE

        # xTB WAVEFUNCTION CALCULATION (xtb4stda)
        command = ["xtb4stda", xyzfile, self.solvent_params]
        output = run(command)

        # xTB EXCITATIONS CALCULATIONS
        command = ["stda", "-xtb", "-e", "8"]
        output = run(command)
        with open("stda-calc.out", "w+") as f:
            f.write(output)
        lines = output.splitlines()
        for idx, line in enumerate(lines):
            if "excitation energies, transition moments" in line:
                ogap = float(lines[idx + 2].split()[1])
                fL = lines[idx + 2].split()[3]

        # Processing results
        vip_dft = convert(vip, constants.DFT_FACTOR["IP"]["High_epsilon"])
        vea_dft = convert(vea, constants.DFT_FACTOR["EA"]["High_epsilon"])
        ogap_dft = convert(ogap, constants.DFT_FACTOR["Gap"]["High_epsilon"])

        egap = vip - vea
        egap_dft = vip_dft - vea_dft

        results = [quotes(self.repeat_smiles), E_xtb, E_solv, vip, vea, egap, ogap, vip_dft, vea_dft, egap_dft, ogap_dft, fL]
        removeJunk()

        return results

def getProps(Id, monomers, style, length, solvent_params, threads, name, database, tableName):
    """
    This function allows for the parallelisation of the screening protocol using the
    screening class, defined above. Each polymer is given its own directory wherein the
    ETKDG and xTB geometries are saved, as well as the logs of the xTB calculations.
    """
    start = int(time.time())
    sid = str(Id)
    log(f"Starting polymer {sid}")

    if not os.path.exists(sid):
        os.mkdir(sid)
    os.chdir(sid)

    ps = polyscreen(Id, monomers, style, length, solvent_params, threads)
    ps.getPolymerWithConformer()
    results = ps.runCalculations()
    os.chdir('..')
    end = int(time.time())
    duration = (end-start)/60

    results.append(duration)
    connection = sql.newConnection(database)
    sql.updateData(connection, tableName, results, Id)
    del connection

    log(f"Done polymer {sid}")
    return results

def getCombinations(monomers, n):
    combinations = list(itertools.combinations_with_replacement(monomers, n))
    enumerated_combinations = enumerate(combinations)
    return combinations

def runScreen(name, monomer_list_raw, style, length, solvent, parallel, database, forceNew = False):
    """
    The screening is initialised, using concurrent.futures from within this function. Presently,
    the number of concurrent jobs can be changed with 'in_parallel', but the total cores and
    OMP_NUM_THREADS must be changed accordingly. OMP_NUM_THREADS is used as that is how xTB
    reads how many threads are to be used.
    """
    # Check solvent validity
    if solvent != None:
        if solvent not in constants.SOLVENTS:
            raise Exception(
                "Invalid solvent choice. Valid solvents:",
                [i for i in constants.SOLVENTS],
            )
        else:
            solvent_params = f"-gbsa {solvent}"
    else:
        solvent_params = ""

    # Conversion of monomers to canonical smiles
    canonical = []
    for m in monomer_list_raw:
        molObj = Chem.MolFromSmiles(m)
        canonical.append(Chem.MolToSmiles(molObj, canonical = True))
    monomer_list = canonical


    # Calculation and database setup/re-setup
    tableName = name
    connection = sql.newConnection(database)
    tableExists = sql.SQLExistenceCheck(connection, tableName)

    if (forceNew and tableExists):
        log(f'forceNew set to True, deleting previous table {tableName} and previous screening directory.')
        sql.dropTable(connectcion, tableName)
        os.rmdir(name)
        
    log(f'Writing to {database}')
    log(f'Table name: {tableName}')
    n = len(set([char for char in style]))
    threads = int(os.environ['OMP_NUM_THREADS'])
    args = []


    if (tableExists == False):
        sql.newTable(connection, tableName, style)

        # Generate all possible combinations of monomers from the provided list.
        combinations = getCombinations(monomer_list, n)
        ids = np.arange(0, len(combinations), 1)

        # If partially populated table already exists, find empty records, read combinations and get IDs
    else:
        df = pd.read_sql_query(f"select * from {tableName} where DurationMins = 0.0", connection)
        columns = list(df.columns)
        monomercols = [x for x in columns if x in list(string.ascii_uppercase)]
        monomers = []
        for col in monomercols:
            monomers.append(df[col].tolist())
        combinations = np.array(monomers).T.tolist()
        ids = df["id"].to_list()

    if len(combinations) == 0:
        log('Quitting. Nothing to do.')
        exit()
    # Combining the arguments for getProps into iterables to be cycled with concurrent.futures.
    chunksize = int(math.ceil(len(combinations) / (parallel)))
    for idx, combination in zip(ids, combinations):
        l = [idx, combination, style, length, solvent_params, threads, name, database, tableName]
        if (tableExists == False):
            empty_data = []
            for c in combination:
                empty_data.append(quotes(c  ))
            empty_data.insert(0, idx)
            empty_data.append([length, "''"] + list(np.zeros(len(constants.DATACOLS)-1)))
            sql.insertData(connection, flatten_list(empty_data), tableName, style)
        args.append(l)

    pwd = os.getcwd()
    if not os.path.exists(name):
        os.mkdir(name)
    os.chdir(name)

    # Launching the parallelised screening protocol
    log("Beginning parallelisation")
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as executor:
        x = executor.map(getProps, *zip(*args), chunksize=chunksize)
    lst = list(x)
    transpose = list(map(list, zip(*lst)))
    os.chdir(pwd)
    log('Done.')
