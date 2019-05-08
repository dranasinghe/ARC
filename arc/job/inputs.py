#!/usr/bin/env python
# encoding: utf-8


##################################################################


"""
parameters for input files:

memory (in MB for gaussian, MW for molpro)
method
basis set
slash is '', unless this is gaussian NOT running a composite method, in which case it is '/'
charge
multiplicity/spin
xyz

gaussian:
    job_type_1: '' for sp, irc, or composite methods, 'opt=calcfc', 'opt=(calcfc,ts,noeigen)',
    job_type_2: '' or 'freq iop(7/33=1)' (cannot be combined with CBS-QB3)
                'scf=(tight,direct) int=finegrid irc=(rcfc,forward,maxpoints=100,stepsize=10) geom=check' for irc f
                'scf=(tight,direct) int=finegrid irc=(rcfc,reverse,maxpoints=100,stepsize=10) geom=check' for irc r
    scan: '\nD 3 1 5 8 S 36 10.000000' (with the line break)
    restricted: '' or 'u' for restricted / unrestricted
    `iop(2/9=2000)` makes Gaussian print the geometry nn eee input orientation even for molecules with more
      than 50 atoms (important so it matches the hessian, and so that Arkane can parse the geometry)

qchem:
    job_type_1: 'opt', 'ts', 'sp'
    job_type_2: 'freq'.
    fine: '\n   GEOM_OPT_TOL_GRADIENT 15\n   GEOM_OPT_TOL_DISPLACEMENT 60\n   GEOM_OPT_TOL_ENERGY 5\n   XC_GRID SG-3'
    restricted: 'false' or 'true' for restricted / unrestricted
"""

input_files = {
    'gaussian': """%chk=check.chk
%mem={memory}mb
%nproc=8

#P {job_type_1} {restricted}{method}{slash}{basis} {job_type_2} {fine} {trsh} iop(2/9=2000)

name

{charge} {multiplicity}
{xyz}
{scan}{scan_trsh}


""",

    'qchem': """$molecule
{charge} {multiplicity}
{xyz}
$end

$rem
   JOBTYPE       {job_type_1}
   METHOD        {method}
   UNRESTRICTED  {restricted}
   BASIS         {basis}{fine}{trsh}
$end

""",

    'molpro': """***,name
memory,{memory},m;
geometry={{angstrom;
{xyz}}}

basis={basis}

int;

{{hf;{shift}
maxit,1000;
wf,spin={spin},charge={charge};}}

{restricted}{method};
{job_type_1}
{job_type_2}
---;

""",

    'mrci': """***,name
memory,{memory},m;
geometry={{angstrom;
{xyz}}}

gprint,orbitals;

basis={basis}

{{hf;shift,-1.0,-0.5;
maxit,1000;
wf,spin={spin},charge={charge};}}

{{multi;
{occ}noextra,failsafe,config,csf;
wf,spin={spin},charge={charge};
natorb,print,ci;}}

{{mrci;
{occ}wf,spin={spin},charge={charge};}}

E_mrci=energy;
E_mrci_Davidson=energd;

table,E_mrci,E_mrci_Davidson;

---;

""",

    'arkane_species': """#!/usr/bin/env python
# -*- coding: utf-8 -*-

linear = {linear}{bonds}

externalSymmetry = {symmetry}

spinMultiplicity = {multiplicity}

opticalIsomers = {optical}

energy = {{'{model_chemistry}': Log('{sp_path}')}}

geometry = Log('{opt_path}')

frequencies = Log('{freq_path}')

{rotors}

""",
    'arkane_hindered_rotor':
        """HinderedRotor(scanLog=Log('{rotor_path}'), pivots={pivots}, top={top}, symmetry={symmetry}, fit='fourier')""",

    'arkane_free_rotor':
        """FreeRotor(pivots={pivots}, top={top}, symmetry={symmetry})""",
    'kinbot':"""{"username": "duminda",
"methodclass": "dft",
"slurm_feature": "",
"verbose": 0,
"rotor_scan": 1,
"single_point_ppn": 1,
"scan_step": 30,
"zf": 4,
"irc_maxpoints": 30,
"homolytic_scissions": 0,
"ppn": 1,
"single_point_qc": "molpro",
"ModelEnergyLimit": 400,
"queuing": "slurm",
"high_level_method": "B3LYP",
"TemperatureList": [500,
2000],
"max_dihed": 5,
"me_code": "mess",
"EnergyRelaxationFactor": 200,
"plot_hir_profiles": 1,
"rotation_restart": 3,
"mess_command": "mess",
"title": "ARC",
"break_bonds": [],
"queue_name": "defq",
"ExcessEnergyOverTemperature": 30,
"pes": 0,
"EnergyRelaxationPower": 0.85,
"charge": 0,
"mult": 2,
"scratch": "/scratch/users/duminda/kinbot",
"EnergyStepOverTemperature": 0.2,
"method": "b3lyp",
"Sigmas": [2.576,
6.0],
"CalculationMethod": "direct",
"high_level_basis": "CBSB7",
"smiles": "CCCCO[O]",
"families": "HO2_Elimination_from_PeroxyRadical",
"simultaneous_kinbot": 5,
"integral": "",
"PressureList": [760],
"irc_stepsize": 20,
"specific_reaction": 0,
"high_level": 0,
"mesmer_command": "mesmer",
"ChemicalEigenvalueMax": 0.2,
"form_bonds": [],
"barrier_threshold": 100.0,
"structure": "",
"me": 1,
"random_conf": 500,
"qc_command": "g09",
"single_point_template": "",
"nrotation": 45,
"conformer_search": 0,
"basis": "6-31G",
"reaction_search": 1,
"queue_template": "",
"Epsilons": [7.08,
310.387],
"qc": "gauss",
"Masses": [4.0,
87.0],
"single_point_key": "MYENERGY",
"EnergyRelaxationExponentCutoff": 15,
"dimer": 0}
"""
}


#