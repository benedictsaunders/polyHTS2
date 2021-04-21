# polyHTS2
Building on [polyHTS by Liam Wilbraham](https://github.com/LiamWilbraham/polyhts), polyHTS2 is an updated version containing more rigourous methods of retreiving data from output files, as well as some improvements to molecular handling and parallelisation. It was developed during my time as a researcher in the Zwijnenburg Group at UCL.

PolyHTS2 has the same prerequisites as its predecessor:
* RDkit
* stk
* xTB 5.6.4SE
* xTB4sTDA
* sTDA

and also
* pandas


## Usage
The screening protocol can be launched with 
```python
  runScreen(
    <monomer_list>, 
    <repeat_unit_style>,      #as defined in stk, such as 'AB' for an AB copolymer
    <number_of_repeat_units>,
    <solvent>,                # as defined in the xtb documentation, and near the top of utils.py
    )
```
and outputs a csv file:
ID | Monomer SMILES | Polymer SMILES | E_xtb | E_solv | VIP | VEA | Opitical gap | Osc. strength
-- | -------------- | -------------- | ----- | ------ | --- | --- | ------------ | -------------
⋮|⋮|⋮|⋮|⋮|⋮|⋮|⋮|⋮


In fact, there will be a separate monomer SMILES column for each of the unique monomers used to create each polymer.

It is known that the conformer search is by far the rate limiting step 😎 of the process, sometimes taking several hours for generation of 500 conformers of an AB copolymer with length 8 on 8 threads. However, this appears to be limited to repeat units with linkages that are no 180° to eachother, which is to be somewhat expected.

## To-do
- [x] Upload to GitHub
- [x] Complete and log initial tests
- [ ] Include commandline capabilities
- [x] Read monomers from user-defined file
- [ ] Allow user to change number of conformers
- [ ] Allow user to retrieve selected geometries as XYZ/MOL files (Probably a good idea to align and minimise RMS!)
- [ ] SQL Integration! CSVs are testing my patience, and pandas is only so good...
- [ ] Provide bencmarks regarding polymer lengths and core counts etc.
