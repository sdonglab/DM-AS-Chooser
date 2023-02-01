import argparse
import re
import os
import json
import csv
import logging
from typing import List, Tuple, Optional, Union
import dataclasses
from dataclasses import dataclass

ACTIVE_SPACE_RE = re.compile(r'(\d+)-(\d+)')
GDM_AS = 'gdm-as'
EDM_AS = 'edm-as'
logger = logging.getLogger('active_space_chooser')


class DipoleNotFoundError(Exception):
    pass


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

    return Parser(_TDDFT_DIPOLE_MOMENT_RULE())


class GDMSelector:
    def __init__(self, mr_calcs: List[MultiRefCalc], ref_dipole: Union[float,
                                                                       str]):
        self.mr_calcs = mr_calcs
        try:
            ref_dipole = float(ref_dipole)
        except ValueError:
            if ref_dipole.endswith('.csv'):
                ref_dipole = self._parse_csv(ref_dipole)
            elif ref_dipole.endswith('.log'):
                ref_dipole = self._parse_tddft_log(ref_dipole)
            else:
                raise ValueError(f'unrecognized file extension {ref_dipole}')

        self.ref_dipole = ref_dipole

    def select(self) -> MultiRefCalc:
        dipole_moments = []
        valid_mr_calcs = []
        dipole_errors = []
        name = 'ref_dipole'
        print(f"{name:15s} -> dipole={self.ref_dipole:.6f} err={0:.6f}")
        for mr_calc in self.mr_calcs:
            basename = os.path.basename(mr_calc.path)
            try:
                dm = self.get_ground_state_dipole(mr_calc.path)
                dipole_moments.append(dm)
                valid_mr_calcs.append(mr_calc)
            except DipoleNotFoundError:
                print(
                    f'did not find dipole moment in {mr_calc.path}; removing from analysis...'
                )
                continue

            err = abs(dm - self.ref_dipole)
            print(f"{basename:15s} -> dipole={dm:.6f} err={err:.6f}")
            dipole_errors.append(err)

        i = dipole_errors.index(min(dipole_errors))
        return valid_mr_calcs[i]

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
        if data is None:
            raise DipoleNotFoundError(f'no dipole for {path}')

        ground_state_idx = 0
        return data[ground_state_idx]['dipole']['total']

    @staticmethod
    def _parse_tddft_log(path: str) -> float:
        parser = get_tddft_parser()
        with open(path, 'r') as f:
            raw_log = f.read()

        dipole = parser.feed(raw_log)
        if dipole is None:
            raise ValueError(f'did not find dipole moment for {path}')

        return dipole['total']

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
                 tddft_calcs: List[str],
                 max_es: int = DEFAULT_MAX_EXCITED_STATE):
        if len(tddft_calcs) < max_es:
            raise ValueError(f'there must be at least {max_es} tddft calcs')
        self.mr_calcs = mr_calcs
        self.tddft_calcs = tddft_calcs
        self.max_es = max_es

    def select(self) -> MultiRefCalc:
        tddft_dipoles = self.get_tddft_es_dipoles(self.tddft_calcs)
        all_mr_dipoles = []
        valid_mr_calcs = []
        all_mr_errors = []
        max_es = self.max_es
        name = 'ref_dipole'
        fmt_dipoles = '  '.join([f'{dm:.6f}' for dm in tddft_dipoles])
        err = 0
        fmt_errors = '  '.join([f'{err:.6f}'] * max_es)
        print(f"{name:15s} -> dipoles=({fmt_dipoles})   err=({fmt_errors})   max_err={err:.6f}")
        for mr_calc in self.mr_calcs:
            basename = os.path.basename(mr_calc.path)
            try:
                dipoles = self.get_mr_es_dipoles(mr_calc.path)
                all_mr_dipoles.append(dipoles)
                valid_mr_calcs.append(mr_calc)
            except DipoleNotFoundError:
                print(f'did not find dipole moments in {mr_calc.path}; removing from analysis...')

            zipped = zip(dipoles[:max_es], tddft_dipoles[:max_es])
            mr_errors = [abs(mr_dm - tddft_dm) for mr_dm, tddft_dm in zipped]
            fmt_dipoles = '  '.join([f'{dm:.6f}' for dm in dipoles])
            fmt_errors = '  '.join([f'{err:.6f}' for err in mr_errors])
            print(f"{basename:15s} -> dipoles=({fmt_dipoles})   err=({fmt_errors})   max_err={max(mr_errors):.6f}")
            all_mr_errors.append(mr_errors)

        max_mr_errors = [max(mr_errors) for mr_errors in all_mr_errors]
        i = max_mr_errors.index(min(max_mr_errors))
        return valid_mr_calcs[i]

    def get_mr_es_dipoles(self, path: str) -> List[float]:
        if path.endswith('.log'):
            dipoles = self._parse_mr_log(path)
        elif path.endswith('.csv'):
            dipoles = self._parse_mr_csv(path)
        else:
            raise ValueError(f'unrecognized file extension {path}')

        if len(dipoles) < self.max_es:
            raise ValueError(
                f'did not find at least {self.max_es} excited state dipole moments in {path}'
            )

        return dipoles

    def get_tddft_es_dipoles(self, paths: List[str]) -> List[float]:
        dipoles = []
        for path in paths:
            if path.endswith('.log'):
                dipole = self._parse_tddft_log(path)
            elif path.endswith('.csv'):
                dipole = self._parse_tddft_csv(path)
            else:
                raise ValueError(f'unrecognized file extension {path}')
            dipoles.append(dipole)

        if len(dipoles) < self.max_es:
            raise ValueError(
                f'did not find at least {self.max_es} excited state dipole moments in {path}'
            )

        return dipoles

    @staticmethod
    def _parse_mr_log(path: str) -> List[float]:
        parser = get_mr_parser()
        with open(path, 'r') as f:
            raw_log = f.read()

        data = parser.feed(raw_log)
        if data is None:
            raise DipoleNotFoundError(
                f'did not find dipole moments for {path}')

        first_es_idx = 1
        return [es['dipole']['total'] for es in data[first_es_idx:]]

    @staticmethod
    def _parse_mr_csv(path: str) -> List[float]:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            dipole_field = reader.fieldnames[0]
            data_rows = list(reader)

        return [float(row[dipole_field]) for row in data_rows]

    @staticmethod
    def _parse_tddft_log(path: str) -> float:
        parser = get_tddft_parser()
        with open(path, 'r') as f:
            raw_log = f.read()

        dipole = parser.feed(raw_log)
        if dipole is None:
            raise ValueError(f'did not find dipole moment for {path}')

        return dipole['total']

    @staticmethod
    def _parse_tddft_csv(path: str) -> float:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            dipole_field = reader.fieldnames[0]
            data_line = next(reader)
            return float(data_line[dipole_field])


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
        required=True,
        help=
        'The reference dipole moment given as a floating point number or path to a TD-DFT calculation file'
    )

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
    files = opts.mr_files.copy()
    if opts.method == EDM_AS:
        files.extend(opts.tddft_files)

    for file in files:
        if not os.path.exists(file):
            parser.error(f'{file} does not exist')

        _, ext = os.path.splitext(file)
        if ext not in valid_exts:
            parser.error(f'{file} must be a .log or .csv file')


def infer_mr_calc(path: str) -> MultiRefCalc:
    num_electrons = None
    num_orbitals = None
    m = ACTIVE_SPACE_RE.search(path)
    if m:
        num_electrons = int(m.group(1))
        num_orbitals = int(m.group(2))
    else:
        logging.warning(f'unable to infer active space for {path}')

    return MultiRefCalc(num_electrons, num_orbitals, path)


def main(args: Optional[List[str]] = None):
    parser, gdm_parser, edm_parser = get_parsers()
    opts = parser.parse_args(args)
    process_opts(gdm_parser, edm_parser, opts)

    mr_calcs = [infer_mr_calc(path) for path in opts.mr_files]
    if opts.method == GDM_AS:
        selector = GDMSelector(mr_calcs, opts.ref_dipole)
    else:
        selector = EDMSelector(mr_calcs, opts.tddft_files)

    selection = dataclasses.asdict(selector.select())
    print(json.dumps(selection))


if __name__ == '__main__':
    main()
