# Coordinates for OH in Input Orientation (angstroms):
#   O    0.0000    0.0000   -0.1221
#   H    0.0000    0.0000    0.8531
conformer(
    label = 'OH',
    E0 = (29.1135, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(17.0027, 'amu')),
        LinearRotor(inertia=(0.901549, 'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([3668.48], 'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

# Thermodynamics for OH:
#   Enthalpy of formation (298 K)   =     9.031 kcal/mol
#   Entropy of formation (298 K)    =    42.578 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300       6.954       9.045      42.624      -3.742
#            400       6.946       9.740      44.624      -8.110
#            500       6.951      10.435      46.174     -12.652
#            600       6.974      11.131      47.443     -17.335
#            800       7.081      12.535      49.462     -27.035
#           1000       7.251      13.967      51.059     -37.092
#           1500       7.720      17.711      54.089     -63.422
#           2000       8.105      21.672      56.365     -91.058
#           2400       8.310      24.958      57.863    -113.912
#    =========== =========== =========== =========== ===========
thermo(
    label = 'OH',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.49682, 0.000188979, -1.03683e-06, 1.65099e-09, -6.50931e-13, 3500.3, 1.48065],
                Tmin = (10, 'K'),
                Tmax = (972.456, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.43865, -0.000262833, 7.26068e-07, -2.88353e-10, 3.55075e-14, 3544.29, 1.92767],
                Tmin = (972.456, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (29.1021, 'kJ/mol'),
        Cp0 = (29.1007, 'J/(mol*K)'),
        CpInf = (37.4151, 'J/(mol*K)'),
    ),
)

