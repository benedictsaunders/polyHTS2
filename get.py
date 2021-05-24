import os, pandas as pd, xarray as xr, numpy as np, matplotlib.pyplot as plt
import itertools, stk
from rdkit import Chem
from rdkit.Chem import Draw, rdmolops

dirs = ["2n", "3n", "4n", "5n", "6n"]
frames = []
homedir = os.getcwd()
for d in dirs:
    os.chdir(d)
    frames.append(pd.read_csv("output.csv"))
    os.chdir(homedir)
concattenated = pd.concat(frames)

x = xr.concat([df.to_xarray() for df in frames], dim="len")

ms = [
    "BrC1=CC2=CC3=C(C=C(Br)N3)C=C2N1",
    "BrC1=CC2=CC3=C(C=C(Br)O3)C=C2O1",
    "BrC1=CC=C(Br)C2=NON=C12",
    "C1(=C(C2=NSN=C2C(=C1F)Br)Br)F",
]
combinations = list(itertools.combinations_with_replacement(ms, 2))

vips = []
veas = []
gaps = []
polymers = []

for l, df in enumerate(frames, 2):
    vips.append(df["vip"].to_list())
    veas.append(df["vea"].to_list())
    gaps.append(df["gap"].to_list())

np_vips = np.array(vips).T
np_veas = np.array(veas).T
np_gaps = np.array(gaps).T

lengths = np.array([2, 3, 4, 5, 6])

fig = plt.figure()
for idx, c in enumerate(combinations):
    ax = fig.add_subplot(2, 5, idx + 1)
    ax.scatter(lengths, np_vips[idx], label="VIP")
    ax.scatter(lengths, np_veas[idx], label="VEA")
    ax.scatter(lengths, np_gaps[idx], label="OG")
    ax.set_title(f"Polymer {str(idx)}")
    ax.set_xlabel("$n$")
    ax.set_ylabel("Energy vs. SHE(eV)")

    bbs = []
    for m in c:
        bbs.append(stk.BuildingBlock(smiles=m, functional_groups=[stk.BromoFactory()]))
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=bbs,
            repeating_unit="AB",
            num_repeating_units=1,
        ),
    )
    H = Chem.MolFromSmarts("[H]")
    Br = Chem.MolFromSmarts("[Br]")
    molObj = polymer.with_canonical_atom_ordering().to_rdkit_mol()
    molObj = rdmolops.RemoveHs(molObj)
    rdmolops.SanitizeMol(molObj)
    molObj = Chem.MolFromSmiles(Chem.MolToSmiles(molObj))
    polymers.append(molObj)


img = Draw.MolsToGridImage(
    polymers,
    molsPerRow=5,
    legends=[
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
    ],
    subImgSize=(250, 250),
)
img.save("key.png")
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper right")
fig.set_size_inches(16, 12)
plt.xticks(lengths)

plt.tight_layout()
plt.savefig("plot.png")
with open("test", "w+") as f:
    f.write("test")
