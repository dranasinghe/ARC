# Coordinates for N2H4 in Input Orientation (angstroms):
#   N   -0.6529    0.0171   -0.3019
#   N    0.6171    0.0550    0.3659
#   H   -1.1314   -0.8216    0.0054
#   H   -1.2389    0.8126   -0.0573
#   H    0.5688    0.5496    1.2542
#   H    1.2586    0.5622   -0.2325
conformer(
    label = 'N2H4',
    E0 = (73.7842, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(32.0375, 'amu')),
        NonlinearRotor(
            inertia = ([3.50009, 20.8275, 20.8456], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(
            frequencies = ([808.247, 988.666, 1108.97, 1290.99, 1321.08, 1663.69, 1680.08, 3401.51, 3413.3, 3512.61, 3519.96], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (0.877501, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-9.99514, 12.2005, -4.0032, -0.0541511, -0.418044],
                    [0.0211491, -0.0515621, 0.0253709, 0.000567373, 0.00432144],
                ],
                'kJ/mol',
            ),
        ),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 2,
)

# Thermodynamics for N2H4:
#   Enthalpy of formation (298 K)   =    18.174 kcal/mol
#   Entropy of formation (298 K)    =    58.089 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      12.193      18.198      58.171       0.747
#            400      13.877      19.502      61.910      -5.262
#            500      15.499      20.972      65.182     -11.620
#            600      17.004      22.598      68.143     -18.288
#            800      19.587      26.266      73.402     -32.455
#           1000      21.654      30.398      78.003     -47.605
#           1500      25.142      42.182      87.511     -89.085
#           2000      27.098      55.287      95.038    -134.789
#           2400      28.080      66.334     100.071    -173.836
#    =========== =========== =========== =========== ===========
thermo(
    label = 'N2H4',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.92977, 0.0050715, 1.19915e-05, -1.6603e-08, 6.58568e-12, 7673.05, 4.93308],
                Tmin = (10, 'K'),
                Tmax = (671.133, 'K'),
            ),
            NASAPolynomial(
                coeffs = [2.56681, 0.0131948, -6.1643e-06, 1.43187e-09, -1.32382e-13, 7856, 10.9651],
                Tmin = (671.133, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (63.7852, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (128.874, 'J/(mol*K)'),
    ),
)

