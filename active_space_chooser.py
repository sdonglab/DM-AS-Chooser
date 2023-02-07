import argparse
import re
import os
import json
import csv
import logging
from typing import List, Tuple, Optional, Union, Dict
import dataclasses
from dataclasses import dataclass

ACTIVE_SPACE_RE = re.compile(r'(\d+)-(\d+)')
EXCITED_STATES_RE = re.compile(r'^([1-9][0-9]*,)*([1-9][0-9]*){1}$')
GDM_AS = 'gdm-as'
EDM_AS = 'edm-as'
CSV_EXT = '.csv'
LOG_EXT = '.log'
VALID_EXTS = (CSV_EXT, LOG_EXT)
NO_MOLEXTRACT_ERROR = 'you must have molextract installed to parse log files'
logger = logging.getLogger('active_space_chooser')
handler = logging.StreamHandler()
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


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
        raise ValueError(NO_MOLEXTRACT_ERROR)

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
        raise ValueError(NO_MOLEXTRACT_ERROR)

    if _TDDFT_DIPOLE_MOMENT_RULE is None:
        _TDDFT_DIPOLE_MOMENT_RULE = general.DipoleMoment

    return Parser(_TDDFT_DIPOLE_MOMENT_RULE())


def parse_mr_log(path: str) -> List[Dict[str, float]]:
    parser = get_mr_parser()
    with open(path, 'r') as f:
        raw_log = f.read()

    try:
        data = parser.feed(raw_log)
    except ValueError:
        # Issue when attempting to parse
        data = None

    if data is None or len(data) == 0:
        raise DipoleNotFoundError(f'no dipole(s) for {path}')

    return data


def parse_tddft_log(path: str) -> Dict[str, float]:
    parser = get_tddft_parser()
    with open(path, 'r') as f:
        raw_log = f.read()

    try:
        dipole = parser.feed(raw_log)
    except ValueError:
        # Issue when attempting to parse
        dipole = None

    if dipole is None:
        raise ValueError(f'did not find dipole moment for {path}')

    return dipole


class GDMSelector:
    def __init__(self, mr_calcs: List[MultiRefCalc], ref_dipole: Union[float,
                                                                       str]):
        self.mr_calcs = mr_calcs
        try:
            ref_dipole = float(ref_dipole)
        except ValueError:
            if ref_dipole.endswith(CSV_EXT):
                ref_dipole = self._parse_csv(ref_dipole)
            elif ref_dipole.endswith(LOG_EXT):
                ref_dipole = parse_tddft_log(ref_dipole)
                ref_dipole = ref_dipole['total']
            else:
                raise ValueError(f'unrecognized file extension {ref_dipole}')

        self.ref_dipole = ref_dipole

    def select(self) -> MultiRefCalc:
        dipole_moments = []
        valid_mr_calcs = []
        dipole_errors = []
        name = 'ref_dipole'
        logger.debug(f"{name:15s} -> dipole={self.ref_dipole:.6f} err={0:.6f}")
        for mr_calc in self.mr_calcs:
            basename = os.path.basename(mr_calc.path)
            try:
                dm = self.get_ground_state_dipole(mr_calc.path)
                dipole_moments.append(dm)
                valid_mr_calcs.append(mr_calc)
            except DipoleNotFoundError:
                logger.warning(
                    f'did not find dipole moment in {mr_calc.path}; removing from analysis...'
                )
                continue

            err = abs(dm - self.ref_dipole)
            logger.debug(f"{basename:15s} -> dipole={dm:.6f} err={err:.6f}")
            dipole_errors.append(err)

        if len(valid_mr_calcs) == 0:
            raise ValueError(
                'no valid multi reference calcs found; cannot do analysis')

        i = dipole_errors.index(min(dipole_errors))
        return valid_mr_calcs[i]

    @staticmethod
    def get_ground_state_dipole(path: str) -> float:
        if path.endswith('.log'):
            dipoles = parse_mr_log(path)
            ground_state_index = 0
            return dipoles[ground_state_index]['dipole']['total']
        elif path.endswith('.csv'):
            return GDMSelector._parse_csv(path)
        else:
            raise ValueError(f'unrecognized file extension {path}')

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
                 es_spec: List[int]):
        if len(tddft_calcs) != len(es_spec):
            raise ValueError(f'mismatch between es_spec and tddft_calcs')
        self.mr_calcs = mr_calcs
        self.tddft_calcs = tddft_calcs
        self.es_spec = es_spec

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
        logger.debug(
            f"{name:15s} -> dipoles=({fmt_dipoles})   err=({fmt_errors})   max_err={err:.6f}"
        )
        for mr_calc in self.mr_calcs:
            basename = os.path.basename(mr_calc.path)
            try:
                dipoles = self.get_mr_es_dipoles(mr_calc.path)
                all_mr_dipoles.append(dipoles)
                valid_mr_calcs.append(mr_calc)
            except DipoleNotFoundError:
                logger.debug(
                    f'did not find dipole moments in {mr_calc.path}; removing from analysis...'
                )
                continue

            zipped = zip(dipoles[:max_es], tddft_dipoles[:max_es])
            mr_errors = [abs(mr_dm - tddft_dm) for mr_dm, tddft_dm in zipped]
            fmt_dipoles = '  '.join([f'{dm:.6f}' for dm in dipoles])
            fmt_errors = '  '.join([f'{err:.6f}' for err in mr_errors])
            logger.debug(
                f"{basename:15s} -> dipoles=({fmt_dipoles})   err=({fmt_errors})   max_err={max(mr_errors):.6f}"
            )
            all_mr_errors.append(mr_errors)

        if len(valid_mr_calcs) == 0:
            raise ValueError(
                'no valid multi reference calcs found; cannot do analysis')

        max_mr_errors = [max(mr_errors) for mr_errors in all_mr_errors]
        i = max_mr_errors.index(min(max_mr_errors))
        return valid_mr_calcs[i]

    def get_mr_es_dipoles(self, path: str) -> List[float]:
        if path.endswith(LOG_EXT):
            dipoles = parse_mr_log(path)
            first_es_idx = 1
            dipoles = [es['dipole']['total'] for es in dipoles[first_es_idx:]]
        elif path.endswith(CSV_EXT):
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
            if path.endswith(LOG_EXT):
                dipole = parse_tddft_log(path)
                dipole = dipole['total']
            elif path.endswith(CSV_EXT):
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
    def _parse_mr_csv(path: str) -> List[float]:
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            dipole_field = reader.fieldnames[0]
            data_rows = list(reader)

        return [float(row[dipole_field]) for row in data_rows]

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
    edm_parser.add_argument(
        '-S',
        type=str,
        default='1,2,3',
        help='A comma seperated string of the provided excited states')
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
    files = opts.mr_files.copy()
    if opts.method == EDM_AS:
        files.extend(opts.tddft_files)
        if not EXCITED_STATES_RE.match(opts.S):
            parser.error(
                "excited state specification must be a comma seperated list of integer"
            )

        opts.S = [int(n) for n in opts.S.split(',')]
        if len(opts.S) != len(opts.tddft_files):
            parser.error(
                f"number of excited states in -S ({len(opts.S)}) do not match the number of tddft reference vals ({len(opts.tddft_files)})"
            )

    for file in files:
        if not os.path.exists(file):
            parser.error(f'{file} does not exist')

        _, ext = os.path.splitext(file)
        if ext not in VALID_EXTS:
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
        selector = EDMSelector(mr_calcs, opts.tddft_files, opts.S)

    selection = dataclasses.asdict(selector.select())
    print(json.dumps(selection))


if __name__ == '__main__':
    main()
