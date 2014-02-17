#!/usr/local/env python

from uvspec.spectrum import AbsorptionSpectrum


AS = AbsorptionSpectrum()

#AS.excited_state_energy = [1,2,3]
#AS.oscillator_strength = [0.1,0.2,0.3]

AS.extract('tests/ethene.log')

AS.generate()
AS.write('test-api')
AS.plot()
