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
import argparse
import configparser
import os.path

from uvspec.version import get_version


class ConfigFile():
    """Process configuration file."""
    def __init__(self):
        self.path = os.path.expanduser('~/.uvspecgen')
        self.filename = os.path.join(self.path, 'uvspecgen.conf')
        self.config = configparser.ConfigParser()
        self.default = {'grid': '0.02',
                        'range': '2.50',
                        'sigma': '0.12',
                        'shift': '0.00'}
        if not os.path.exists(self.path):
            self.create()

    def create(self):
        try:
            os.mkdir(self.path)
        except OSError:
            pass
        self.config['FitParameters'] = self.default
        self._write_config()

    def read(self):
        self.config.read(self.filename) 
        params = self.config['FitParameters']
        # Convert the parameter values from strings to floats
        params = dict([k, float(v)] for k, v in params.iteritems())
        return params 

    def reset(self):
        self.create()

    def update(self, parameter, value):
        self.config.read(self.filename)
        self.config['FitParameters'][parameter] = str(value)
        self._write_config()

    def _write_config(self):
        with open(self.filename, 'w') as configfile:
            self.config.write(configfile)


class Settings():
    """Process command-line input and generate help documentation.

    An instance of this class parses command-line input and provides program
    help documentation.  All needed user-input is handled with this class,
    including input/output file names, fit parameters, and output formatting.

    """
    def __init__(self):
        self.version = get_version()
        self.header = '%(prog)s' + ' ' + self.version

        defaults = ConfigFile().read()
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
        self.save = self.clinput.save
        self.reset = self.clinput.reset

    def _parse_command_line_input(self, defaults):
        """Create an ArgumentParser object and define program options.

        The defaults dictionary contains the default values for several
        parameters required to perform the Gaussian fit.  The dictionary
        is defined at the top of this file so that users may easily set
        the default values themself.  Additionally, the dictionary allows
        a single object to be passed into the CommandLineInput class and
        the herein described function.

        """
        parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Generate UV-Vis spectrum from TDHF/TDDFT data',
            epilog = 'Report bugs to <gaussiantoolkit@gmail.com>')
        
        # Required positional argument, the input (.log) file
        parser.add_argument(
            'logfile',
            nargs = '*',
            help = 'input TDHF/TDDFT logfile')
        
        # Option to specify the output file name
        parser.add_argument(
            '-O',
            '--outfile',
            default = None,
            help = ('output filename, .spec.txt will be appended, '
                    '<logfile>.spec.txt is used as default'))

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
        
        # Options to modify the output data
        parser.add_argument(
            '--nometa',
            #default = 'False',
            action = 'store_true',
            help = 'do not print spectrum metadata to output file')
        parser.add_argument(
            '-o',
            '--output',
            default = 'both',
            choices = ['both', 'curve', 'sticks'],
            help = 'specify results to print in output file')
        
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
            help = ('set the spacing below and above the smallest and '
                    'largest excited state energy for plotting'))
        parser.add_argument(
            '-s',
            '--sigma',
            default = defaults['sigma'],
            type = float,
            help = 'set value for Gaussian broadening constant sigma')
        parser.add_argument(
            '-S',
            '--shift',
            default = defaults['shift'],
            type = float,
            help = 'set shift for the starting point of the Gaussian curve')
       
        # Flags for configuring program-wide default fit parameters
        config = parser.add_mutually_exclusive_group()
        config.add_argument(
            '--save',
            action = 'store_true',
            help = 'save the specified fit parameters as the default values')
        config.add_argument(
            '--reset',
            action = 'store_true',
            help = 'reset the fit parameters to the original default values')

        # Print program version
        parser.add_argument(
            '--version',
            action = 'version',
            version = self.header)
        return parser.parse_args()
