#!/usr/bin/env python
# encoding: utf-8

name = "minimal"
shortDesc = u""
longDesc = u"""

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
            NASAPolynomial(coeffs=[3.49863,8.15038e-05,-2.56227e-07,5.54135e-12,3.17271e-13,-1061.95,-4.27706], Tmin=(10,'K'), Tmax=(606.905,'K')),
            NASAPolynomial(coeffs=[3.68632,-0.000783318,9.61266e-07,-3.21285e-10,3.56255e-14,-1091.59,-5.14531], Tmin=(606.905,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (-8.83004,'kJ/mol'),
        Cp0 = (29.1007,'J/(mol*K)'),
        CpInf = (37.4151,'J/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Bond corrections: {'H-H': 1}
Treated as a linear species


Geometry:
H       1.19804800    0.01874300    0.06873000
H       1.94230700    0.01874300    0.06873000
""",
)

