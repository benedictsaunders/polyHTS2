import os
import subprocess as sp
import numpy as np
from rdkit import Chem
import pandas as pd
import rdkit
from rdkit.Chem import rdmolops
from time import sleep

valid_solvents = [
    "acetone",
    "acetonitrile",
    "benzene",
    "chcl3",
    "cs2",
    "dmso",
    "ether",
    "h2o",
    "methanol",
    "thf",
    "toluene",
]

def log(s, homedir):
    pwd = os.getcwd()
    os.chdir(homedir)
    with open('output.log', 'a') as logfile:
        logfile.write(s)
        logfile.write('\n')
    os.chdir(pwd)

def box(msg):
    row = len(msg) + 2
    h = "".join(["+"] + ["-" * row] + ["+"])
    result = h + "\n" "| " + msg + " |" "\n" + h
    return result

def run(command):
    sleep(1)
    p = sp.Popen(calc_params, stdout=sp.PIPE, encoding='utf8')
    output, _ = p.communicate()
    return output

def remove_junk():
    os.system('rm charges charges3 energy tda.dat wfn.xtb wbo xtbopt* xtbrestart .xtboptok')
