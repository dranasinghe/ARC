# Coordinates for methylamine in Input Orientation (angstroms):
#   N    0.8243    0.2062   -0.3489
#   C   -0.5698   -0.0228    0.0434
#   H   -1.0848   -0.5630   -0.7554
#   H   -0.7138   -0.5880    0.9781
#   H   -1.0755    0.9400    0.1548
#   H    1.3055   -0.6810   -0.4569
#   H    1.3141    0.7086    0.3847
conformer(
    label = 'methylamine',
    E0 = (-36.4512, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(31.0422, 'amu')),
        NonlinearRotor(
            inertia = ([4.88657, 22.2735, 23.1584], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([844.739, 969.64, 1045.46, 1160.05, 1339.64, 1447, 1482.2, 1503.91, 1652.12, 2918.29, 3021.17, 3056.11, 3453.02, 3526.63], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (1.12175, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [-0.0635472, 0.0718386, -4.47214, -0.0681254, 0.0447407],
                    [0.101423, 0.124489, -0.471867, 0.0931588, 0.0965851],
                ],
                'kJ/mol',
            ),
        ),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

# Thermodynamics for methylamine:
#   Enthalpy of formation (298 K)   =    -5.917 kcal/mol
#   Entropy of formation (298 K)    =    57.803 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      11.958      -5.893      57.882     -23.258
#            400      14.199      -4.586      61.627     -29.237
#            500      16.440      -3.053      65.038     -35.572
#            600      18.535      -1.303      68.223     -42.237
#            800      22.110       2.775      74.064     -56.476
#           1000      24.946       7.492      79.315     -71.823
#           1500      29.632      21.257      90.414    -114.363
#           2000      32.153      36.766      99.320    -161.873
#           2400      33.365      49.886     105.296    -202.825
#    =========== =========== =========== =========== ===========
thermo(
    label = 'methylamine',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.00472, -0.000514021, 3.57069e-05, -4.40694e-08, 1.7711e-11, -4384.43, 5.19361],
                Tmin = (10, 'K'),
                Tmax = (644.485, 'K'),
            ),
            NASAPolynomial(
                coeffs = [0.915559, 0.0186585, -8.91525e-06, 2.08772e-09, -1.93353e-13, -3986.24, 18.74],
                Tmin = (644.485, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (-36.4551, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (153.818, 'J/(mol*K)'),
    ),
)

