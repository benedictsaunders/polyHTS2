# polyHTS2
Building on [polyHTS by Liam Wilbraham](https://github.com/LiamWilbraham/polyhts), polyHTS2 is an updated version containing more rigourous methods of retreiving data from output files, as well as some improvements to molecular handling and parallelisation. It was developed during my time as a researcher in the Zwijnenburg Group at UCL.

PolyHTS2 has the following prerequisites:
* RDkit
* stk
* xTB 5.6.4SE
* xTB4sTDA
* sTDA
* pandas
* sqlite3

## What does it do?

Given a list of monomers in the form of SMILES strings, polyHTS2 finds all the possible combinations of these monomers, which are then used to define a polymer. Each polymer undergoes a conformer search; every conformer is relaxed with MMFF94. The lowest energy conformer is then treated with xTB (SE geometry optimisation, VIP, VEA) and also sTDA (optical gap, oscillator strength).

## Usage

The screen can be started from the command line/terminal with the following arguments:
```
python screen.py -f <monomers_list>     # List of SMILES strings, default is monomers.txt
                 -l <polymer length>    # Number of repeat units, not separate monomers, default 4
                 -r <repeat_style>      # E.g. A for a homopolymer, AB for a binary copolymer, default AB
                 -p <parallel_workers>  # Should be $NSLOTS/OMP_NUM_THREADS
                 -n <environment_name>  # A new directory where everything is written and saved, default 'screen'
                 -s <solvent>           # Deafult 'h2o'
```

Or, in your own python script, running a screen is done simply with
```python
  runScreen(
    <monomer_list>, 
    <repeat_unit_style>,
    <number_of_repeat_units>,
    <solvent>,                # as defined in the xtb documentation, and near the top of constants.py
    )
```
You must set `OMP_NUM_THREADS` otherwise polyHTS will use as many cores as are available. This will likely results in massive inefficiencies, particularly if you specify more that one parallel worker.

The default values are set by `argparse` when handing command line arguments, and not when executing `runScreen`.

Currently, the output CSV has the form:
ID | Monomer SMILES | Polymer SMILES | E_xtb | E_solv | VIP | VEA | Opitical gap | Osc. strength|Duration|
-- | -------------- | -------------- | ----- | ------ | --- | --- | ------------ | -------------|--------|
⋮|⋮|⋮|⋮|⋮|⋮|⋮|⋮|⋮|⋮

PolyHTS will also write to a SQL database during execution, using `SQLite3`, stored in the `<name>.db` file.


In fact, there will be a separate monomer SMILES column for each of the unique monomers used to create each polymer.

It is known that the conformer search is by far the rate limiting step 😎 of the process, sometimes taking several hours for generation of 500 conformers of an AB copolymer with length 8 on 8 threads. This issue is especially prevelant for larger monomers, and apperas to be driven by RDKit struggling to embed the monomers properly from SMILES. I shall look in to why, and if there is a workaround.

## To-do
- [x] Upload to GitHub
- [x] Complete and log initial tests
- [x] Include commandline capabilities
- [x] Read monomers from user-defined file
- [ ] Allow user to change number of conformers
- [ ] Allow user to retrieve selected geometries as XYZ/MOL files (Probably a good idea to align and minimise RMS!)
- [ ] ~~PostreSQL Integration! CSVs are testing my patience, and pandas is only so good...~~
- [x] SQLite Integration!
- [ ] Provide bencmarks regarding polymer lengths and core counts etc.
- [x] Add conversion/scaling for SE <-> DFT comparisons
- [ ] Look into why some combinations lead to slow conformer searching.
