# Coordinates for H2 in Input Orientation (angstroms):
#   H    1.1980    0.0187    0.0687
#   H    1.9423    0.0187    0.0687
conformer(
    label = 'H2',
    E0 = (-8.82033, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(2.01565, 'amu')),
        LinearRotor(inertia=(0.279128, 'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([4373.9], 'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

# Thermodynamics for H2:
#   Enthalpy of formation (298 K)   =    -0.035 kcal/mol
#   Entropy of formation (298 K)    =    31.137 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300       6.961      -0.022      31.183      -9.376
#            400       6.953       0.674      33.185     -12.600
#            500       6.947       1.369      34.735     -15.999
#            600       6.950       2.064      36.002     -19.537
#            800       7.005       3.458      38.007     -26.948
#           1000       7.111       4.869      39.581     -34.712
#           1500       7.492       8.516      42.533     -55.284
#           2000       7.878      12.362      44.743     -77.125
#           2400       8.115      15.563      46.202     -95.321
#    =========== =========== =========== =========== ===========
thermo(
    label = 'H2',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.49863, 8.15038e-05, -2.56227e-07, 5.54135e-12, 3.17271e-13, -1061.95, -4.27706],
                Tmin = (10, 'K'),
                Tmax = (606.905, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.68632, -0.000783318, 9.61266e-07, -3.21285e-10, 3.56255e-14, -1091.59, -5.14531],
                Tmin = (606.905, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (-8.83004, 'kJ/mol'),
        Cp0 = (29.1007, 'J/(mol*K)'),
        CpInf = (37.4151, 'J/(mol*K)'),
    ),
)

