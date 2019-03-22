#!/usr/bin/env python
# -*- coding: utf-8 -*-

linear = False

bonds = {'N-N': 1, 'H-N': 4}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 2

energy = {'cbs-qb3': Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/N2H4/composite_a1700/output.out')}

geometry = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/N2H4/freq_a1701/output.out')

frequencies = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/N2H4/freq_a1701/output.out')



rotors = [HinderedRotor(scanLog=Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/N2H4/scan_a1693/output.out'), pivots=[1, 2], top=[1, 3, 4], symmetry=1, fit='fourier')]

