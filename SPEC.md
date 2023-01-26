# CLI / Help Messages
```
$ python mles-choser.py -h
usage: mles-chooser <command> [args]

    gdm-as     Run the ground-state dipole moment active-space selection algorithim
    emd-as     Run the excited-state dipole moment active-space selection algorithim

See `mles-choose.py gdm-as -h` and `mles-chooser.py edm-as -h` for more information
```

## GDM-AS

```
$ python mles-choser.py gdm-as -h
useage: mles-choser.py gdm-as [--csv] [--data-dir DATA_DIR] --ref-dipole REF_DIPOLE

    -c, --csv              Read the data in DATA_DIR as CSV files instead of calculation log files
    -d, --data-dir         The the directory containing the OpenMolcas multi-reference calculations (default $PWD)
    -r, --ref-dipole       The reference dipole moment given as a floating point number

The DATA_DIR must contain sub-directories whose names are in the form <n_electrons>-<n_orbitals>. Within
each of the <n_electrons>-<n_orbitals> directories, there must be a file that ends in .log that is the
OpenMolcas log file of the corresponding CASSCF Calculation. For example say you provide --data-dir
/foo/bar the contents of the /foo/bar may look like

/foo/bar
    11-10/
        11-10.log
    5-5/
        5-5.log
    7-12/
        benzene_7_12.log

If the --csv flag is provided the same directory structure applies but the log files must instead be files
that end in .csv. The file must have a single column header of any name, and the ground state dipole
moment followed on the next line. For example it may look like

Ground_State_Dipole_Moment
0.06
```
Example usages
- `python mles-choose.py gdm-as --data-dir /scratch/s.dong/test-set --ref-dipole 0.05`
- `python mles-choose.py gdm-as --ref-dipole 0.05`
- `python mles-choose.py gdm-as --csv --data-dir /scratch/s.dong/test-set --ref-dipole 0.05`
- `python mles-choose.py gdm-as --csv --ref-dipole 0.05`

## EDM-AS
```
$ python mles-choser.py edm-as -h
useage: mles-choser.py gdm-as [--csv] [--data-dir DATA_DIR] [--ref-tddft REF_TDDFT]

    -c, --csv              Read the data in DATA_DIR and REF_TDDFT as CSV files instead of calculation log files
    -d, --data-dir         The the directory containing the OpenMolcas multi-reference calculations (default $PWD)
    -r, --ref-tddft        The reference Gaussian TD-DFT calculation log file (default to the .log file in DATA_DIR)

The DATA_DIR must contain sub-directories whose names are in the form <n_electrons>-<n_orbitals>. Within
each of the <n_electrons>-<n_orbitals> directories, there must be a file that ends in .log that is the
OpenMolcas log file of the corresponding CASSCF Calculation. Unless an explicit path is provided for REF_TDDFT
it will be assumed to be the file ending in .log located in DATA_DIR For example say you provide --data-dir
/foo/bar the contents of the /foo/bar may look like

/foo/bar
    11-10/
        11-10.log
    5-5/
        5-5.log
    7-12/
        benzene_7_12.log
    tddft.log

If the --csv flag is provided the same directory structure applies but the log files must instead be files
that end in .csv. The CSV file must have a single column header of any name. On the next three lines must
follow the S1, S2 and S3 excited state dipole moments. For example it may look like

Excited_State_Dipole_Moment
0.01
0.02
0.03
```
Example usages
- `python mles-choose.py edm-as --data-dir /scratch/s.dong/test-set --ref-tddft /work/s.dong/tddft.log`
- `python mles-choose.py edm-as --data-dir /scratch/s.dong/test-set`
- `python mles-choose.py edm-as --ref-tddft /work/s.dong/tddft.log`
- `python mles-choose.py edm-as`
- `python mles-choose.py edm-as --csv --data-dir /scratch/s.dong/test-set --ref-tddft /work/s.dong/tddft.log`
- `python mles-choose.py edm-as --csv --data-dir /scratch/s.dong/test-set`
- `python mles-choose.py edm-as --csv --ref-tddft /work/s.dong/tddft.log`
- `python mles-choose.py edm-as --csv`

# Output Format
The output of either script will be single JSON object of the following format
```
{
    "num_electrons": 12,
    "num_orbitals": 10
}
```