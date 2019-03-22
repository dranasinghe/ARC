# Coordinates for vinoxy in Input Orientation (angstroms):
#   O    0.9864   -0.1541    0.2605
#   C    2.8154    1.2670    0.5785
#   C    2.1641    0.1261    0.0279
#   H    3.8523    1.4873    0.3555
#   H    2.2569    1.9247    1.2333
#   H    2.7729   -0.5177   -0.6355
conformer(
    label = 'vinoxy',
    E0 = (2.10971, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(43.0184, 'amu')),
        NonlinearRotor(
            inertia = ([7.50143, 44.1211, 51.6225], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([444.03, 499.256, 746.906, 964.384, 969.909, 1139.28, 1381.96, 1454.54, 1536.87, 2899.38, 3110.88, 3226.05], 'cm^-1'),
        ),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

# Thermodynamics for vinoxy:
#   Enthalpy of formation (298 K)   =     3.319 kcal/mol
#   Entropy of formation (298 K)    =    61.895 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      12.728       3.344      61.980     -15.250
#            400      15.277       4.747      65.998     -21.652
#            500      17.432       6.386      69.646     -28.437
#            600      19.236       8.222      72.988     -35.571
#            800      22.129      12.373      78.939     -50.778
#           1000      24.252      17.022      84.117     -67.095
#           1500      27.356      30.025      94.619    -111.904
#           2000      28.911      44.126     102.722    -161.318
#           2400      29.764      55.869     108.072    -203.504
#    =========== =========== =========== =========== ===========
thermo(
    label = 'vinoxy',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.06896, -0.00648651, 8.39881e-05, -1.46829e-07, 8.48452e-11, 254.379, 7.29715],
                Tmin = (10, 'K'),
                Tmax = (527.075, 'K'),
            ),
            NASAPolynomial(
                coeffs = [1.8099, 0.018851, -1.14378e-05, 3.36298e-09, -3.82274e-13, 378.706, 15.6694],
                Tmin = (527.075, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (2.11095, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (133.032, 'J/(mol*K)'),
    ),
)

