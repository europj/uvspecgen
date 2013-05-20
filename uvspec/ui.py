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

"""UVSpecGen

This is the main user interface module for UVSpecGen.  An instance of the
ArgumentParser class is created for providing program usage and parsing
command-line input.

"""

import sys
if sys.version_info > (2, 7):
    import argparse as arg
else:
    try:
        from lib import argparse as arg
    except ImportError:
        print " [ERROR] The Python module ArgParse is required"
        sys.exit(1)

class CommandLineInput():
    """Process command-line input and generate help documentation.

    An instance of this class parses command-line input and provides program
    help documentation.  All needed user-input is handled with this class,
    including input/output file names, fit parameters, and output formatting.

    """
    def __init__(self, version, default):
        self.header = '%(prog)s' + ' ' + version
        self.input = self.parse_command_line_input(self.header, default)
        self.logfile = self.input.logfile
        self.outfile = self.input.outfile
        self.plot = self.input.plot
        self.output = self.input.output
        self.nometa = self.input.nometa
        self.parameters = dict(
            grid = self.input.grid,
            range = self.input.range,
            sigma = self.input.sigma,
            shift = self.input.shift)


    def parse_command_line_input(self, header, default):
        """Create an ArgumentParser object and define program options.

        The default dictionary contains the default values for several
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
            help = 'Gaussian logfile')
        
        # Optional positional argument, the output file name; option to modify
        # output data
        parser.add_argument(
            'outfile',
            nargs = '?',
            default = '<logfile>.spec.txt',
            help = 'specify output filename, .spec.txt will be appended')
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
            default = default['grid'],
            type = float,
            help = 'set grid spacing value for the Gaussian curve')
        parser.add_argument(
            '-r',
            '--range',
            default = default['range'],
            type = float,
            help = 'set the spacing below and above the smallest and largest\
                    excited state energy for plotting')
        parser.add_argument(
            '-s',
            '--sigma',
            default = default['sigma'],
            type = float,
            help = 'set value for Gaussian broadening constant sigma')
        parser.add_argument(
            '--shift',
            default = default['shift'],
            type = float,
            help = 'set shift for the starting point of the Gaussian curve')

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
            version = header)
        args = parser.parse_args()
        return args
