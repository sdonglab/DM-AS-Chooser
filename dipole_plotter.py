"""
Plots a series of multireference calculation dipole moments
"""
import argparse
from string import digits
import os
import logging
import pathlib

from molextract.rules.molcas import log, general, rasscf
from molextract.parser import Parser

import matplotlib.pyplot as plt

LOG_EXT = '.log'
DEFAULT_COLOR = 'pink'
SIMPLE_CPK_COLORS = {
    'H': 'white',
    'C': 'black',
    'N': 'blue',
    'O': 'red',
    'F': 'green',
    'CL': 'green',
    'P': 'orange',
    'S': 'yellow',
}
ATOM_MARKER_SIZE = 250

logger = logging.getLogger('active_space_chooser')
handler = logging.StreamHandler()
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


class RASSCFModule(log.ModuleRule):
    def __init__(self):
        rules = [rasscf.RASSCFCartesianCoords(), general.MolProps()]
        super().__init__('rasscf', rules)


class DipolePlotter:
    def __init__(self, mr_files):
        self.mr_files = mr_files
        self.fig = plt.figure()
        # self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')

    def plot(self):
        # Move left y-axis and bottom x-axis to centre, passing through (0,0)
        # self.ax.spines['left'].set_position('zero')
        # self.ax.spines['bottom'].set_position('zero')

        # Hide upper and right axes
        # self.ax.spines['right'].set_visible(False)
        # self.ax.spines['top'].set_visible(False)
        # self.ax.set_ylim(bottom=-.4e-6, top=.4e-6)
        # self.ax.set_xlim(left=-.4e-6, right=.4e-6)


        # Show ticks in the left and lower axes only
        # self.ax.xaxis.set_ticks_position('bottom')
        # self.ax.yaxis.set_ticks_position('left')

        plot_coords = True
        for path in self.mr_files:
            stem = pathlib.Path(path).stem
            try:
                coords, mol_props = self._parse_mr_file(path)
            except ValueError:
                continue

            if plot_coords:
                self._plot_coords(coords)
                plot_coords = False

            ground_dipole = mol_props[0]['dipole']
            xs = (0, ground_dipole['x'])
            ys = (0, ground_dipole['y'])
            zs = (0, ground_dipole['z'])
            self.ax.plot(xs, ys, zs, 'o-', label=stem)

        self.ax.set_xlabel('X Dipole (Ang)')
        self.ax.set_ylabel('Y Dipole (Ang)')
        self.ax.set_zlabel('Z Dipole (Ang)')
        self.ax.legend()

    def _parse_mr_file(self, path):
        parser = Parser(RASSCFModule())
        print(path)
        try:
            with open(path, 'r') as f:
                raw_data = f.read()
            parsed_data = parser.feed(raw_data)
        except ValueError:
            logger.warning(f'unabled to parse {path}, removing from plot')
            raise

        if parsed_data is None or len(parsed_data) == 0:
            logger.warning(f'no data parsed from {path}, removing from plot')
            raise ValueError('no data parsed')

        coords, mol_props = parsed_data
        # TODO: verify coords are consistent across files
        return coords, mol_props 
    
    def _plot_coords(self, coords):
        xs, ys, zs = [], [], []
        colors = []
        for element, x, y, z in coords:
            colors.append(self._get_color(element))
            xs.append(x)
            ys.append(y)
            zs.append(z)
        
        self.ax.scatter(xs, ys, zs, s=ATOM_MARKER_SIZE, c=colors)
    
    def _get_color(self, element):
        element = element.rstrip(digits)
        return SIMPLE_CPK_COLORS.get(element, DEFAULT_COLOR)


    def show(self):
        plt.show()

    def save(self, path):
        self.fig.savefig(path)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'mr_files',
        type=str,
        nargs='+',
        help='The path(s) to the multi-reference calculation log files')

    return parser


def process_opts(parser, opts):
    for file in opts.mr_files:
        if not os.path.exists(file):
            parser.error(f'{file} does not exit')
        if not file.endswith(LOG_EXT):
            parser.error(f'{file} must be .log file')


def main(args=None):
    parser = get_parser()
    opts = parser.parse_args(args)
    print(args, opts)
    process_opts(parser, opts)

    plotter = DipolePlotter(opts.mr_files)
    plotter.plot()
    plotter.show()
    plotter.save('foo.png')


if __name__ == '__main__':
    main()
