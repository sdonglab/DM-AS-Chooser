import argparse
import re
import os
import glob
import csv
import logging
from typing import List, Tuple, Optional
from dataclasses import dataclass

MR_DIRNAME_RE = re.compile(r'\d+-\d+')
GDM_AS = 'gdm-as'
EDM_AS = 'edm-as'
logger = logging.getLogger('active_space_chooser')


@dataclass
class MultiRefCalc:
    num_electrons: int
    num_orbitals: int
    path: str


_RASSCF_MOL_PROPS_RULE = None
_TDDFT_DIPOLE_MOMENT_RULE = None


def get_mr_parser():
    """
    Lazily import molextract and return a new parser to parse RASSCF Mol Props
    """
    global _RASSCF_MOL_PROPS_RULE
    try:
        from molextract.rules.molcas import log, general
        from molextract.parser import Parser
    except ModuleNotFoundError:
        raise ValueError(
            'you must have molextract installed to parse log files')

    if _RASSCF_MOL_PROPS_RULE is None:

        class RASSCFMolProps(log.ModuleRule):
            def __init__(self):
                super().__init__('rasscf', [general.MolProps()])

            def reset(self):
                return self.rules[0].reset()

        _RASSCF_MOL_PROPS_RULE = RASSCFMolProps

    return Parser(_RASSCF_MOL_PROPS_RULE())


def get_tddft_parser():
    """
    Lazily import molextract and return a new parser to parse Gaussian dipole
    moments
    """
    global _TDDFT_DIPOLE_MOMENT_RULE
    try:
        from molextract.rules.gaussian import general
        from molextract.parser import Parser
    except ModuleNotFoundError:
        raise ValueError(
            'you must have molextract installed to parse log files')

    if _TDDFT_DIPOLE_MOMENT_RULE is None:
        _TDDFT_DIPOLE_MOMENT_RULE = general.DipoleMoment

    return Parser(_TDDFT_DIPOLE_MOMENT_RULE)


class GDMSelector:
    def __init__(self, mr_calcs: List[MultiRefCalc], ref_dipole: float):
        self.mr_calcs = mr_calcs
        self.ref_dipole = ref_dipole

    def select(self) -> MultiRefCalc:
        dipole_moments = [
            self.get_ground_state_dipole(mr_calc.path)
            for mr_calc in self.mr_calcs
        ]

        dipole_errors = [abs(dm - self.ref_dipole) for dm in dipole_moments]
        i = dipole_errors.index(min(dipole_errors))

        return self.mr_calcs[i]

    @staticmethod
    def get_ground_state_dipole(path: str) -> float:
        if path.endswith('.log'):
            return GDMSelector._parse_mr_log(path)
        elif path.endswith('.csv'):
            return GDMSelector._parse_csv(path)
        else:
            raise ValueError(f'unrecognized file extension {path}')

    @staticmethod
    def _parse_mr_log(path: str) -> float:
        parser = get_mr_parser()
        with open(path, 'r') as f:
            raw_log = f.read()

        data = parser.feed(raw_log)
        ground_state_idx = 0
        return data[ground_state_idx]['dipole']['total']

    @staticmethod
    def _parse_csv(path: str) -> float:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            dipole_field = reader.fieldnames[0]
            data_line = next(reader)
            return float(data_line[dipole_field])


class EDMSelector:
    DEFAULT_MAX_EXCITED_STATE = 3

    def __init__(self,
                 mr_calcs: List[MultiRefCalc],
                 ref_tddft: str,
                 max_es: int = DEFAULT_MAX_EXCITED_STATE):
        self.mr_calcs = mr_calcs
        self.ref_tddft = ref_tddft
        self.max_es = max_es

    def select(self) -> MultiRefCalc:
        tddft_dipoles = self.get_tddft_es_dipoles(self.ref_tddft)
        all_mr_dipoles = [
            self.get_mr_es_dipoles(mr_calc.path) for mr_calc in self.mr_calcs
        ]

        all_mr_errors = []
        for mr_dipoles in all_mr_dipoles:
            mr_errors = [
                abs(mr_dm - tddft_dm)
                for mr_dm, tddft_dm in zip(mr_dipoles, tddft_dipoles)
            ]
            all_mr_errors.append(mr_errors)

        max_mr_errors = [max(mr_errors) for mr_errors in all_mr_errors]
        i = max_mr_errors.index(min(max_mr_errors))

        return self.mr_calcs[i]

    def get_mr_es_dipoles(self, path: str) -> List[float]:
        if path.endswith('.log'):
            dipoles = self._parse_mr_log(path)
        elif path.endswith('.csv'):
            dipoles = self._parse_csv(path)
        else:
            raise ValueError(f'unrecognized file extension {path}')

        if len(dipoles) != self.max_es:
            raise ValueError(
                f'did not find {self.max_es} excited state dipole moments in {path}'
            )

        return dipoles

    def get_tddft_es_dipoles(self, path: str) -> List[float]:
        if path.endswith('.log'):
            dipoles = self._parse_tddft_log(path)
        elif path.endswith('.csv'):
            dipoles = self._parse_csv(path)
        else:
            raise ValueError(f'unrecognized file extension {path}')

        if len(dipoles) != self.max_es:
            raise ValueError(
                f'did not find {self.max_es} excited state dipole moments in {path}'
            )

        return dipoles

    @staticmethod
    def _parse_mr_log(path: str) -> List[float]:
        parser = get_mr_parser()
        with open(path, 'r') as f:
            raw_log = f.read()

        data = parser.feed(raw_log)
        first_es_idx = 1
        return [es['dipole']['total'] for es in data[first_es_idx:]]

    @staticmethod
    def _parse_tddft_log(path: str) -> List[float]:
        raise NotImplementedError

    @staticmethod
    def _parse_csv(path: str) -> List[float]:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            dipole_field = reader.fieldnames[0]
            data_rows = list(reader)

        return [float(row[dipole_field]) for row in data_rows]


def get_parsers() -> Tuple[argparse.ArgumentParser, argparse.ArgumentParser,
                           argparse.ArgumentParser]:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='The available methods to choose',
                                       dest='method',
                                       required=True)

    gdm_parser = subparsers.add_parser(
        GDM_AS,
        help=
        'Run the ground-state dipole moment active-space selection algorithim')
    gdm_parser.add_argument(
        '-m',
        '--mr-files',
        type=str,
        required=True,
        nargs='+',
        help='The path(s) to the multi-reference calculation files')
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
        '-m',
        '--mr-files',
        type=str,
        required=True,
        nargs='+',
        help='The path(s) to the multi-reference calculation files')
    edm_parser.add_argument('-t',
                            '--tddft-files',
                            type=str,
                            required=True,
                            nargs='+',
                            help='The reference TD-DFT log files ')

    return parser, gdm_parser, edm_parser


def process_opts(gdm_parser: argparse.ArgumentParser,
                 edm_parser: argparse.ArgumentParser,
                 opts: argparse.Namespace):
    parser = gdm_parser if opts.method == GDM_AS else edm_parser
    valid_exts = ('.log', '.csv')
    for mr_file in opts.mr_files:
        if not os.path.exists(mr_file):
            parser.error(f'{mr_file} does not exist')

        _, ext = os.path.splitext(mr_file)
        if ext not in valid_exts:
            parser.error(f'{mr_file} must be a .log or .csv file')

    if opts.method == EDM_AS:
        if len(opts.tddft_files) == 1:
            td_file = opts.tddft_files[0]
            if not td_file.endswith('.csv'):
                parser.error(
                    'if a single TD-DFT file is given it must be a .csv file')
        else:
            for td_file in opts.tddft_files:
                if not td_file.endswith('.log'):
                    parser.error(
                        'if more than one TD-DFT file is given they must all be .log files'
                    )

        for td_file in opts.tddft_files:
            if not os.path.exists(td_file):
                parser.error(f'{td_file} does not exist')


def main(args: Optional[List[str]] = None):
    parser, gdm_parser, edm_parser = get_parsers()
    opts = parser.parse_args(args)
    process_opts(gdm_parser, edm_parser, opts)
    print(opts)


if __name__ == '__main__':
    main()
