import argparse
import re
import os
import glob
# TODO: try / except import molextract

MR_DIRNAME_RE = re.compile(r'\d+-\d+')
GDM_AS = 'gdm-as'
EDM_AS = 'edm-as'

class GDMSelector:
    pass

class EDMSelector:
    pass


def get_parsers():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        help='The available methods to choose', dest='method', required=True)

    gdm_parser = subparsers.add_parser(
        GDM_AS,
        help=
        'Run the ground-state dipole moment active-space selection algorithim')
    gdm_parser.add_argument(
        '-c',
        '--csv',
        action='store_true',
        help=
        'Read the data in DATA_DIR as CSV files instead of calculation log files'
    )
    gdm_parser.add_argument(
        '-d',
        '--data-dir',
        type=str,
        default=None,
        help=
        'The the directory containing the OpenMolcas multi-reference calculations (default $PWD)'
    )
    gdm_parser.add_argument(
        '-r',
        '--ref-dipole',
        type=float,
        required=True,
        help='The reference dipole moment given as a floating point number')

    edm_parser = subparsers.add_parser(
        EDM_AS,
        help=
        'Run the excited-state dipole moment active-space selection algorithim'
    )
    edm_parser.add_argument(
        '-c',
        '--csv',
        action='store_true',
        help=
        'Read the data in DATA_DIR and REF_TDDFT as CSV files instead of calculation log files'
    )
    edm_parser.add_argument(
        '-d',
        '--data-dir',
        type=str,
        default=None,
        help=
        'The the directory containing the OpenMolcas multi-reference calculations (default $PWD)'
    )
    edm_parser.add_argument(
        '-r',
        '--ref-tddft',
        type=str,
        default=None,
        help=
        'The reference Gaussian TD-DFT calculation log file (default to the .log file in DATA_DIR)'
    )

    return parser, gdm_parser, edm_parser


def process_opts(gdm_parser, edm_parser, opts):
    parser = gdm_parser if opts.method == GDM_AS else edm_parser
    ext = 'csv' if opts.csv else 'log'
    if opts.data_dir is None:
        opts.data_dir = os.getcwd()
    
    if not os.path.exists(opts.data_dir):
        parser.error(f'{opts.data_dir} does not exist')
    
    mr_file_glob = os.path.join(opts.data_dir, '*', f'*.{ext}')
    mr_files = []
    for path in glob.glob(mr_file_glob):
        if os.path.isfile(path):
            parent_dir = os.path.basename(os.path.dirname(path))
            if MR_DIRNAME_RE.fullmatch(parent_dir):
                mr_files.append(path)
    
    opts.mr_files = mr_files
    if len(mr_files) == 0:
        parser.error(f'did not find any multi-reference calculation files in {opts.data_dir}')
    
    if opts.method == EDM_AS and opts.ref_tddt is None:
        tddft_file_glob = os.path.join(opts.data_dir, f'*.{ext}')
        for path in glob.glob(tddft_file_glob):
            if os.path.isfile(path):
                opts.ref_tddft = path
                break

        parser.error(f'did not find td-dft calculation file in {opts.data_dir}')

    if opts.method == EDM_AS and not os.path.exists(opts.ref_tddft):
        parser.error(f'{opts.ref_tddft} does not exist')
    
    return opts
    

def main(args=None):
    parser, gdm_parser, edm_parser = get_parsers()
    opts = parser.parse_args(args)
    opts = process_opts(gdm_parser, edm_parser, opts)


if __name__ == '__main__':
    main()
