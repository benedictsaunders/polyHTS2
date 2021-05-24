import os
import subprocess as sp
import numpy as np
from rdkit import Chem
import pandas as pd
import rdkit
from rdkit.Chem import rdmolops
import time
from time import sleep


def listelementsToString(lst):
    return [str(i) for i in lst]


def quotes(s):
    return f"'{s}'"


def log(s):
    with open("output.log", "a") as logfile:
        t = time.strftime("%H:%M:%S")
        logfile.write(f"[{t}]: {s}")
        logfile.write("\n")


def box(msg):
    row = len(msg) + 2
    h = "".join(["+"] + ["-" * row] + ["+"])
    result = h + "\n" "| " + msg + " |" "\n" + h
    return result


def run(command):
    sleep(1)
    p = sp.Popen(command, stdout=sp.PIPE, encoding="utf8")
    output, _ = p.communicate()
    return output


def removeJunk():
    os.system(
        "rm charges charges3 energy tda.dat wfn.xtb wbo xtbopt* xtbrestart .xtboptok"
    )


def flatten_list(nd_list):
    flat_list = []
    for i in nd_list:
        if type(i) is list:
            for item in i:
                flat_list.append(item)
        else:
            flat_list.append(i)
    return flat_list


def convert(x, ab):
    a = ab[0]
    b = ab[1]
    return (float(x) * a) + b
