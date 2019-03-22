#!/usr/bin/env python
# encoding: utf-8

name = "minimal"
shortDesc = u""
longDesc = u"""
ARC v1.0.0
ARC project minimal

Levels of theory used:

Conformers:       b97d/6-31g
TS guesses:       b3lyp/6-31+g(d,p)
Optimization:     b3lyp/cbsb7 (using a fine grid)
Frequencies:      b3lyp/cbsb7
Single point:     b3lyp/6-311+g(3df,2p)
Rotor scans:      b3lyp/cbsb7
Using bond additivity corrections for thermo

Using the following settings: {'molpro': 'pharos', u'ssh': True, 'qchem': 'pharos', 'gaussian': 'c3ddb'}

Considered the following species and TSs:
Species H2 (run time: 0:01:15)

Overall time since project initiation: 00:02:17
"""
entry(
    index = 0,
    label = "H2",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.49732,0.000137375,-6.62237e-07,9.20434e-10,-3.07911e-13,-1462.6,-4.2753], Tmin=(10,'K'), Tmax=(1121.84,'K')),
            NASAPolynomial(coeffs=[3.57841,-0.000557471,8.09322e-07,-2.76439e-10,3.06534e-14,-1455.26,-4.56204], Tmin=(1121.84,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (-12.1613,'kJ/mol'),
        Cp0 = (29.1007,'J/(mol*K)'),
        CpInf = (37.4151,'J/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Bond corrections: {'H-H': 1}
Treated as a linear species

External symmetry: 2, optical isomers: 1

Geometry:
H       1.19359800   -0.03722400    0.06557400
H       1.93777100   -0.03722400    0.06557400
""",
)

