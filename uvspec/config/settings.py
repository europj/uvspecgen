# Generate UV-Vis spectra from Gaussian09 TDHF/TDDFT log files. 
# Copyright (C) 2013  Gaussian Toolkit
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Command line user interface for the UVSpecGen script.

An instance of the ArgumentParser class from the argparse module is created
for providing program usage information and for parsing command-line input.

"""
from uvspec import get_version

# The argparse module is not included with Python versions less than 2.7.
# A copy of the module is included with the uvspec package.
try:
    import argparse as arg
except ImportError:
    from uvspec.lib import argparse as arg


defaults = dict(
    grid  = 0.02,
    range = 2.50,
    sigma = 0.12,
    shift = 0.00)


class CommandLineInput():
    """Process command-line input and generate help documentation.

    An instance of this class parses command-line input and provides program
    help documentation.  All needed user-input is handled with this class,
    including input/output file names, fit parameters, and output formatting.

    """
    def __init__(self, defaults=defaults):
        self.version = get_version()
        self.header = '%(prog)s' + ' ' + self.version

        self.clinput = self._parse_command_line_input(defaults)

        self.logfile = self.clinput.logfile
        self.outfile = self.clinput.outfile
        self.join = self.clinput.join
        self.plot = self.clinput.plot
        self.output = self.clinput.output
        self.nometa = self.clinput.nometa
        self.parameters = {'grid': self.clinput.grid,
                           'range': self.clinput.range,
                           'sigma': self.clinput.sigma,
                           'shift': self.clinput.shift}

    def _parse_command_line_input(self, defaults):
        """Create an ArgumentParser object and define program options.

        The defaults dictionary contains the default values for several
        parameters required to perform the Gaussian fit.  The dictionary
        is defined at the top of this file so that users may easily set
        the default values themself.  Additionally, the dictionary allows
        a single object to be passed into the CommandLineInput class and
        the herein described function.

        """
        parser = arg.ArgumentParser(
            formatter_class = arg.ArgumentDefaultsHelpFormatter,
            description = 'Generate UV-Vis spectrum from TDHF/TDDFT data',
            epilog = 'Report bugs to <gaussiantoolkit@gmail.com>')
        
        # Required positional argument, the input (.log) file
        parser.add_argument(
            'logfile',
            nargs = '+',
            help = 'Gaussian logfile')
        
        # Optional positional argument, the output file name; option to
        # modify output data
        parser.add_argument(
            'outfile',
            nargs = '?',
            default = None,
            help = ('output filename, .spec.txt will be appended, '
                    '<logfile>.spec.txt is used as default'))
        parser.add_argument(
            '-o',
            '--output',
            default = 'both',
            choices = ['both', 'curve', 'sticks'],
            help = 'specify results to print in output file')
        parser.add_argument(
            '--nometa',
            default = 'False',
            action = 'store_true',
            help = 'do not print spectrum metadata to output file')
        
        # Optional arguments, to modify the parameters of the Gaussian fit
        parser.add_argument(
            '-g',
            '--grid',
            default = defaults['grid'],
            type = float,
            help = 'set grid spacing value for the Gaussian curve')
        parser.add_argument(
            '-r',
            '--range',
            default = defaults['range'],
            type = float,
            help = 'set the spacing below and above the smallest and largest\
                    excited state energy for plotting')
        parser.add_argument(
            '-s',
            '--sigma',
            default = defaults['sigma'],
            type = float,
            help = 'set value for Gaussian broadening constant sigma')
        parser.add_argument(
            '--shift',
            default = defaults['shift'],
            type = float,
            help = 'set shift for the starting point of the Gaussian curve')

        # Option to join multiple log files into one spectrum
        parser.add_argument(
            '-j',
            '--join',
            action = 'store_true',
            help = 'join multiple logfile spectra into one spectrum')
        
        # Optional flag to plot the spectrum using matplotlib, if available
        parser.add_argument(
            '-p',
            '--plot',
            action = 'store_true',
            help = 'display a plot of the absorbance spectrum if matplotlib\
                    is available')
        
        # Print program version
        parser.add_argument(
            '--version',
            action = 'version',
            version = self.header)
        return parser.parse_args()
