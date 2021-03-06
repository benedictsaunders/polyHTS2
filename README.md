# polyHTS2
Building on [polyHTS by Liam Wilbraham](https://github.com/LiamWilbraham/polyhts), polyHTS2 is an updated version containing more rigourous methods of retreiving data from output files, as well as some improvements to molecular handling and parallelisation. It was developed during my time as a researcher in the Zwijnenburg Group at UCL.

PolyHTS2 has the following prerequisites:
* RDkit
* stk
* xTB 5.6.4SE
* xTB4sTDA
* sTDA
* sqlite3

## What does it do?

PolyHTS2 builds a series of polymers for computatiol screening. Polymers can be created exhaustively (all possbile combinations of the monomers inputted), or have a set pattern, e.g. give one list for monomer A and another for monomer B. Each polymer undergoes a conformer search; every conformer is relaxed with MMFF94. The lowest energy conformer is then treated with xTB (SE geometry optimisation, VIP, VEA) and also sTDA (optical gap, oscillator strength).

## Usage

The screen can be started from the command line/terminal with the following arguments:
```bash
$ python screen.py -f <list_1.txt> <list_2.txt> # Lists of SMILES strings, default is monomers.txt
                 -l <polymer length>            # Number of repeat units, not separate monomers, default 4
                 -r <repeat_style>              # E.g. A for a homopolymer, AB for a binary copolymer, default AB
                 -p <parallel_workers>          # Should be $NSLOTS/OMP_NUM_THREADS
                 -n <environment_name>          # A new directory where the calculation subdirectories are created, default 'screen'
                 -s <solvent>                   # Default 'h2o'
                 -d <database_path>             # If a database does not exist, one will be created here, default 'database.db'
                 --exhaustive                   # If not set, more than one input file is required.
                 --ForceNew                     # Deletes any exisitng directory and table with the same environment name, if they exist
```

Or, in your own python script, running a screen is done simply with
```python
  runScreen(
    <monomer_lists>, 
    <repeat_unit_style>,
    <number_of_repeat_units>,
    <solvent>,                # as defined in the xtb documentation, and near the top of constants.py
    <parallel_workers>,
    <database_path>,
    <exhaustive>,             # Bool
    <forceNew>                # Bool
    )
```
You must set `OMP_NUM_THREADS` otherwise polyHTS will use as many cores as are available. This will likely results in massive inefficiencies, particularly if you specify more that one parallel worker.

### Converion Script

So that previous screens' tables can be compared properly, one can use the `canonical_conversion.py` script. It for every column headered with A, B, C etc, is converted to canonical SMILES, which is added as another column `Canonical_[...]`

```bash
$ python canonical_conversion.py -d <database, default data.db> -t <table>
```
Species undergoing screening will automatically have their SMILES converted.
## Output

The output SQLite database table has the following fields (columns):
* ID
* Monomer smiles (each having a separate column for the combination chosen, e.g. for an ABCB copolymer, there will be three columns, labelled A, B and C)
* xTB single point energy of optimised geometry
* xTB SCC energy, inclusive of solvation
* xTB VIP, VEA,  energy gap and sTDA optical gap
* Approximated DFT VIP, VEA, fundamental gap and optical gap
* sTDA oscialltor strength
* Duration of each polymer calculation

It is known that the conformer search is by far the rate limiting step ???? of the process, sometimes taking several hours for generation of 500 conformers of an AB copolymer with length 8 on 8 threads. This issue is especially prevelant for larger monomers, and apperas to be driven by RDKit struggling to embed the monomers properly from SMILES. I shall look in to why, and if there is a workaround.

## To-do
- [x] Upload to GitHub
- [x] Complete and log initial tests
- [x] Include commandline capabilities
- [x] Read monomers from user-defined file
- [ ] Allow user to change number of conformers
- [ ] Allow user to retrieve selected geometries as XYZ/MOL files (Probably a good idea to align and minimise RMS!)
- [ ] ~~PostreSQL Integration! CSVs are testing my patience, and pandas is only so good...~~
- [x] SQLite Integration!
- [x] Add conversion/scaling for SE <-> DFT comparisons
- [x] Look into why some combinations lead to slow conformer searching.
- [x] Add xTB and aDFT bandgap columns, rename Gap to OpticalGap
- [x] One-to-many combinations

### Reason for longer conformer searches
Running ETKDG with experiemental torsion angles enabled significanly slowed the screenin; an 8 x fluorene polymer on 8 cores took over 2000 minutes!
