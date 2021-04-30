DATACOLS = ["smiles", "E_xtb", "E_solv", "VIPvSHE_xtb", "VIPvSHE_aDFT", "VEAvSHE_xtb", "VEAvSHE_aDFT", "Gap_xtb", "Gap_aDFT", "fL", "DurationMins"]

SOLVENTS = ["acetone", "acetonitrile", "benzene", "chcl3", "cs2", "dmso", "ether", "h2o", "methanol", "thf", "toluene"]

ELECTRODE = -4.44

DFT_FACTOR = {
                "IP":{"High_epsilon":(0.88, -0.24),
                        "Low_epsilon":(0.91, 0.15)},
                "EA":{"High_epsilon":(0.84, -0.74),
                        "Low_epsilon":(0.92, -1.05)},
                "IP/EA":{"High_epsilon":(1.02, -0.41),
                        "Low_epsilon":(1.27, -0.32)},
                "Gap":{"High_epsilon":(0.85, -0.25),
                        "Low_epsilon":(0.83, -0.21)},
            }
