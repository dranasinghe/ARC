# Coordinates for H2 in Input Orientation (angstroms):
#   H    1.1936   -0.0372    0.0656
#   H    1.9378   -0.0372    0.0656
conformer(
    label = 'H2',
    E0 = (-12.1513, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(2.01565, 'amu')),
        LinearRotor(inertia=(0.279063, 'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([4273.39], 'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

# Thermodynamics for H2:
#   Enthalpy of formation (298 K)   =    -0.832 kcal/mol
#   Entropy of formation (298 K)    =    31.136 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300       6.958      -0.818      31.183     -10.172
#            400       6.950      -0.122      33.183     -13.396
#            500       6.948       0.573      34.734     -16.794
#            600       6.956       1.268      36.001     -20.333
#            800       7.012       2.663      38.008     -27.743
#           1000       7.124       4.076      39.584     -35.508
#           1500       7.522       7.735      42.546     -56.083
#           2000       7.908      11.596      44.764     -77.933
#           2400       8.143      14.809      46.228     -96.139
#    =========== =========== =========== =========== ===========
thermo(
    label = 'H2',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.49732, 0.000137375, -6.62237e-07, 9.20434e-10, -3.07911e-13, -1462.6, -4.2753],
                Tmin = (10, 'K'),
                Tmax = (1121.84, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.57841, -0.000557471, 8.09322e-07, -2.76439e-10, 3.06534e-14, -1455.26, -4.56204],
                Tmin = (1121.84, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (-12.1613, 'kJ/mol'),
        Cp0 = (29.1007, 'J/(mol*K)'),
        CpInf = (37.4151, 'J/(mol*K)'),
    ),
)

