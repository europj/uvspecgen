#!/usr/bin/python

# uvspecgen 
# Generate a UV-vis absorption spectrum from a Gaussian TDHF/TDDFT calculation

from math import *
try:
    from matplotlib.pyplot import *
    mkplt = True
except ImportError:
    print 'matplotlib.pyplot is required to plot the spectrum'
    mkplt = False


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
        
    newfilename=logfile[:-3]+'spec'                 # write new file with information gleaned
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

logfilename=raw_input('Enter logfile: ')            # run this function as a script-prompting for input

genAbsSpec(logfilename)

    
