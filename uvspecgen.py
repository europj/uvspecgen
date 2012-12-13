#!/usr/bin/python

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
#                                   UVSPECGEN
#
# Copyright (C) 2012; Li Research Group; University of Washington; Seattle, WA
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

'''
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*                               FILE DESCRIPTION
*
* Filename      : uvspecgen.py
* Version       : 1.0.0
* Programmer(s) : JWM
* Updated on    : 12.12.12
* Created on    : 07.11.12
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
'''

# SET PROGRAM DEFAULTS HERE
defaults = dict(
    grid  = 0.02,
    range = 2.50,
    sigma = 0.12,
    shift = 0.00)


__version__ = '1.0.0'

#
# MODULES
#

import argparse as arp
import sys

from math import *
try:
    from matplotlib.pyplot import *
    mkplt = True
except ImportError:
    print 'matplotlib.pylot is required to plot the spectrum'
    mkplt = False

#
# CLASSES
#

class CommandLineInput():
    def __init__(self, version=__version__, default=defaults):
        self.header = '%(prog)s' + ' ' + version
        self.input = self._parse_command_line_input(self.header, default)
        self.logfile = self.input.logfile
        self.outfile = self.input.outfile
        self.fit_parameters = dict(
            grid = self.input.grid,
            range = self.input.range,
            sigma = self.input.sigma,
            shift = self.input.shift)


    def _parse_command_line_input(self, header, default):
        """_parse_command_line_input(header : string, default : dictionary) ->
        [list]

        Create an ArgumentParser object and define program options.  The
        default dictionary contains the default values for several parameters
        required to perform the Gaussian fit.  The dictionary is defined at the
        top of this file so that users may easily set the default values
        themself.  Additionally, the dictionary allows a single object to be
        passed into the CommandLineInput class and the herein described
        function.

        """
        parser = arp.ArgumentParser(
            formatter_class = arp.ArgumentDefaultsHelpFormatter,
            description = 'Generate UV-vis spectrom from TDHF/TDDFT data',
            epilog = 'Report bugs to <gaussiantoolkit@gmail.com>')
        
        # Required positional argument, the input (.log) file
        parser.add_argument(
            'logfile',
            type = arp.FileType('r'))
        
        # Optional positional argument, the output file name
        parser.add_argument(
            'outfile',
            nargs = '?',
            default = gen_outfile,
            type = arp.FileType('w'))
        
        # Optional arguments, to control the parameters of the Gaussian fit
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
        
        # Print program version
        parser.add_argument(
            '--version',
            action = 'version',
            version = header)
        args = parser.parse_args()
        return args

#
# LOCAL FUNCTIONS
#

def gen_outfile(logfile):
    """gen_outfile(logfile : string) -> [string]
    
    Removes the log extension from the input file (of the form filename.log)
    and appends the spec extension for the output file

    """
    outfile = logfile[:-3]+'spec'
    return outfile


def genAbsSpec(logfile):
    ''' This function takes in the direction to a .log file and creates a .spec text file
            in the same location containing the excited states from the logfile and their
            oscillator strengths, as well as the coordinates for points of a spectrum
            (osc. strength vs. energy). It also generates this graph. 
    '''
    Emin=0.0                                        # Energy shift
    sigma=0.12                                      # Broadening parameter
    dEstep=0.02                                     # Step between energy points calculated for graph
    dErange=2.5                                     # Range of graph/output is between dErange less than the smallest excited state
                                                    # energy and dErange greater than the largest excited state energy (see *)
    log=open(logfile)                               # open file
    ExcitedState=[]                                 # create lists for excited states, oscilator strength, and energy and absorbanvce vectors for graph
    Intensity=[]
    Evec=[] 
    A=[]                                    
    
    for line in log:                                # read logfile and extract the excited states and oscilator strengths, storing
        if line.startswith(' Excited State '):      # them in the appropriate lists
            words=line.split(' ')
            words=delspaces(words)
            ExcitedState+=[float(words[4])]
            Intensity+=[float(words[8][2:])]
            
    maxE=dErange+max(ExcitedState)                  # *
    minE=min(ExcitedState)-dErange
    n=minE
    while n<=maxE:                                  # create energy points for graph
        Evec+=[n]
        n+=dEstep

    for En in Evec:                                 # create absorbance spectra points for each energy point for graph
        a=0
        for index in range(len(ExcitedState)):
            a+=Intensity[index]*e**(-0.5*((En+Emin-ExcitedState[index])**2)/(sigma**2))
        A+=[a]
        
#    newfilename=logfile[:-3]+'spec'                 # write new file with information gleaned
    newfile=open(newfilename,'w')
    sn=1
    newfile.write(repr('Excited States:')[1:-1].rjust(18))
    newfile.write(repr('Energy Level(eV)')[1:-1].rjust(22))
    newfile.write(repr('Intensity')[1:-1].rjust(20))
    newfile.write('\n')
    for state in range(len(ExcitedState)):
        newfile.write(repr('Excited State')[1:-1].rjust(15))
        newfile.write(repr(sn).rjust(4))
        newfile.write(repr(ExcitedState[state]).rjust(20))
        newfile.write(repr(Intensity[state]).rjust(20))
        newfile.write('\n')
        sn+=1
    newfile.write('\n')
    newfile.write(repr('Gaussian Fit Absorbance Energy(eV)')[1:-1].rjust(35))
    newfile.write(repr('Intensity')[1:-1].rjust(30))
    newfile.write('\n')
    for level in range(len(Evec)):
        newfile.write(repr(Evec[level]).rjust(25))
        newfile.write(repr(A[level]).rjust(40))
        newfile.write('\n')
        
    if mkplt:
        xlabel("Energy(eV)")                            # creat absorbance graph
        ylabel("Oscillator Strength")
        plot(Evec,A,'k')
        show()

    
  
def delspaces(L):
    ''' This function takes a list and recursively removes any empty strings from the list.
            (returning the list without empty strings)
    '''
    if len(L)==0:
        return L
    elif L[0]=='':
        return delspaces(L[1:])
    elif L[0]!='':
        return [L[0]]+delspaces(L[1:])

#
# MAIN PROGRAM
#

options = CommandLineInput()
print options.logfile
print options.outfile
print options.fit_parameters['grid']
print options.fit_parameters['range']
print options.fit_parameters['sigma']
print options.fit_parameters['shift']
outfile = gen_outfile(options.input.logfile.name)
print outfile

#logfilename=raw_input('Enter logfile: ')            # run this function as a script-prompting for input

#genAbsSpec(logfilename)
