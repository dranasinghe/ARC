#!/usr/bin/env python
# -*- coding: utf-8 -*-

linear = False

bonds = {'C=C': 1, 'C-C': 1, 'C-H': 6}

externalSymmetry = 1

spinMultiplicity = 1

opticalIsomers = 1

energy = {'cbs-qb3': Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/propene/composite_a1691/output.out')}

geometry = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/propene/freq_a1694/output.out')

frequencies = Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/propene/freq_a1694/output.out')



rotors = [HinderedRotor(scanLog=Log('/home/dranasinghe/Software/ARC/examples/thermo_demo/calcs/Species/propene/scan_a1695/output.out'), pivots=[1, 2], top=[1, 4, 5, 6], symmetry=3, fit='fourier')]

