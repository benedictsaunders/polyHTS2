DATACOLS = [
    "smiles",
    "E_xtb",
    "E_solv",
    "VIPvSHE_xtb",
    "VEAvSHE_xtb",
    "EnergyGap_xtb",
    "OpticalGap_xtb",
    "VIPvSHE_aDFT",
    "VEAvSHE_aDFT",
    "EnergyGap_aDFT",
    "OpticalGap_aDFT",
    "fL",
    "DurationMins",
]

SOLVENTS = [
    ("acetone", 20.7),
    ("acetonitrile", 37.5),
    ("benzene", 2.27),
    ("chcl3", 4.81),
    ("cs2", 2.6),
    ("dmso", 46.68),
    ("ether", 4.33),
    ("h2o", 80.1),
    ("methanol", 32.7),
    ("thf", 7.58),
    ("toluene", 2.38),
]

ELECTRODE = -4.44

DFT_FACTOR = {
    "IP": {"High_epsilon": (0.88, -0.24), "Low_epsilon": (0.91, 0.15)},
    "EA": {"High_epsilon": (0.84, -0.74), "Low_epsilon": (0.92, -1.05)},
    "IP/EA": {"High_epsilon": (1.02, -0.41), "Low_epsilon": (1.27, -0.32)},
    "Gap": {"High_epsilon": (0.85, -0.25), "Low_epsilon": (0.83, -0.21)},
}
