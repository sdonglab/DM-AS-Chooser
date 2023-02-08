# active-space-chooser
Active Space Chooser (ASC) implements the recommended active space selection algorithm
as described in this [paper](TODO).

ASC has two methods of selection, GDM-AS and EDM-AS see the paper for a more thorough
description of these methods.

## CLI
### Input
The input to ASC will depend on using GDM-AS or EDM-AS, though both share similar interfaces.
All multi-reference calculations are assumed to be have a single file per active space. The name
of the file is encouraged to have the syntax `<active_electrons>-<active_orbitalS>` in the name
though this is not required. For example some file names could be
- `mr-14-14.log`
- `mr_benzene_10-10.log`
- `4-4.csv`
- `benzene.csv`

The multi-reference calculations files can be provided either as log file or as csv files. If they
are log files they must be `Molcas`/`OpenMolcas` log files. ASC will use [molextract](https://github.com/sdonglab/molextract)
to automatically parser these log files and extract the relevant dipole moments.

If your multi-reference calculations were not done using `Molcas`/`OpenMolcas` then the multi-reference
calculation files must instead be CSV files. Refer to the [GDM-AS](#gdm-as) and [EDM-AS](#edm-as)
sections for what the expected syntax of these CSV files.

You are allowed to mix-and-match both log and csv files.

Additionally both GDM-AS and EDM-AS both have an input reference value for the reference dipole moment.
Refer to the [GDM-AS](#gdm-as) and [EDM-AS](#edm-as) for specifics as to how to provide these reference
values.


### Output
Based on the method used, there will be logging messages indicated the values being extract / analyzed
printed to `stderr`. If a multi-reference calculation file does not have any data it will be removed
from the analysis and an appropriate log message will be printed.

The final selected active space will be printed to `stdout` as a JSON with the following format
```
{"num_electrons": <n_elec>, "num_orbitals": <n_orb>, "path": <path to file>}
```

The `num_electrons` and `num_orbitals` will simply be the inferred values from the name of the file
as described in the beginning of the Input section. If they cannot be inferred they will be set to
`nul`. The `path` field will contain the path to the multi-reference calculation that was selected,
again refer to the paper for a reference as to how this selection is made.


## GDM-AS
```
usage: active_space_chooser.py gdm-as [-h] -m MR_FILES [MR_FILES ...] -r REF_DIPOLE

optional arguments:
  -h, --help            show this help message and exit
  -m MR_FILES [MR_FILES ...], --mr-files MR_FILES [MR_FILES ...]
                        The path(s) to the multi-reference calculation files
  -r REF_DIPOLE, --ref-dipole REF_DIPOLE
                        The reference dipole moment given as a floating point number or path to a TD-DFT calculation file
```
To use GDM-AS provide the paths to the multi-reference calculation files as well as provide a reference
ground-state dipole moment. This reference value can be provided in three ways
- A single floating point number (eg. `1.2345`)
- A path to the Gaussian ground-state dipole calculation log file (must end in `.log`)
- A path to a CSV file containing the ground-state dipole (see below; must end in `.csv`)

### CSV Syntax
The CSV file must contain a header line, it does not matter what this header is or how many columns there are.
There can be any number of data lines following the header, but the value that will be read as the total ground
state dipole will be the value in the first row and first column, below is an example CSV file
```csv
Ground_State_Dipole
1.2345
```

## EDM-AS
```
usage: active_space_chooser.py edm-as [-h] -m MR_FILES [MR_FILES ...] [-S S [S ...]] -t TDDFT_FILES [TDDFT_FILES ...]

optional arguments:
  -h, --help            show this help message and exit
  -m MR_FILES [MR_FILES ...], --mr-files MR_FILES [MR_FILES ...]
                        The path(s) to the multi-reference calculation files
  -S S [S ...]          which ground / excited states to use
  -t TDDFT_FILES [TDDFT_FILES ...], --tddft-files TDDFT_FILES [TDDFT_FILES ...]
                        The reference excited state calculation files
```
To use EDM-AS provide the paths to the multi-reference calculation files as well as provide reference
excited-state dipole moment calculation files.

By default EDM-AS will run the algorithm using the first three excited state dipole moments (S1, S2, S3). The
specific excited states to use in the algorithm can be specified with the `-S` flag. Provide the index for which
states to run the analysis on with `-S`. These correspond as
- `0` -> ground state
- `1` -> S1
- `2` -> S2
- etc.

For example specifying `-S 3 4` will run the EDM-AS algorithm for S3 and S4. Using `-S 0` will run the EDM-AS
algorithm using just S0 or the ground state. The numbers passed into `-S` must be unique, non-negative and may
be provided in any order. By default the script will use `-S 1 2 3`.

The values passed into `-S` are key to how the `MR_FILES` and `TDDFT_FILES` are parsed. For the `MR_FILES` it is expected
that each file will contain excited states up to and include the largest excited state specified in `-S`. So if `-S 3` is
provided and a OpenMolcas log files is provided in `MR_FILES` it is expected that the log file contains dipole moments
for S0, S1, S2 and S3.

The `TDDFT_FILES` can be the paths to the corresponding excited state dipole moment calculations. These files should be in
the same order as specified in `-S`. See the section on the CSV [syntax](#td-dft-csv-syntax).

To re-cap the order of `TDDFT_FILES` is expected to match the order of `-S` so if `-S 2 1 4` is provided then the `TDDFT_FILES`
must be the excited state dipole moments for S2, S1 and S4

## Multi-Reference CSV Syntax
The `MR_FILES` syntax for CSV files is also closely tied to the values provided in `-S`. The CSV file must have a header,
though again the contents of the header do not matter. It is then expected that the CSV file will have dipole moments for
the ground state up to the largest excited state specified in `-S`. Each row following the header line will contain the
total dipole moment for that corresponding ground / excited state. So if `-S 1 2 3` is set an example CSV is
```csv
Excited_State_Dipole
1.000
1.500
1.750
2.000
```

If `-S 3` is provided you will still require four rows (S0, S1, S2, S3) though the values in S0, S1, and S2 do not matter
and may be dummy data for example
```
Excited_State_Dipole
-1
-1
-1
2.000
```

## TD-DFT CSV Syntax
These CSV files have similar syntax to the ground state dipole moment CSV syntax in GDM-AS. Again there must be a header
with any name. The excited state dipole moment will then be the value located in the first row and first column. An
example is below
```csv
Excited_State_Dipole
1.2345
```


## Installation
If you only want to use CSV files for running ASC there are no dependencies to run ASC. If you would like to run ASC
with the convenience of parsing OpenMolcas and Gaussian log files you must install [molextract](https://github.com/sdonglab/molextract)
