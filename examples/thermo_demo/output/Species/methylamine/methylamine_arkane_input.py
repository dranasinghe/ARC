#!/usr/bin/env python
# -*- coding: utf-8 -*-

linear = False

bonds = {'C-N': 1, 'H-N': 2, 'C-H': 3}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {'cbs-qb3': Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/methylamine/composite_a1685/output.out')}

geometry = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/methylamine/freq_a1689/output.out')

frequencies = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/methylamine/freq_a1689/output.out')



rotors = [HinderedRotor(scanLog=Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/methylamine/scan_a1690/output.out'), pivots=[1, 2], top=[1, 7, 6], symmetry=3, fit='fourier')]

