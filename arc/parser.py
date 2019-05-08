#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function, unicode_literals)
import logging
import numpy as np
import os

from arkane.statmech import Log

from arc.species.converter import get_xyz_string
from arc.arc_exceptions import InputError

"""
Various ESS parsing tools
"""

##################################################################


def parse_frequencies(path, software):
    if not os.path.isfile(path):
        raise InputError('Could not find file {0}'.format(path))
    freqs = np.array([], np.float64)
    if software.lower() == 'qchem':
        with open(path, 'rb') as f:
            for line in f:
                if ' Frequency:' in line:
                    items = line.split()
                    for i, item in enumerate(items):
                        if i:
                            freqs = np.append(freqs, [(float(item))])
    elif software.lower() == 'gaussian':
        with open(path, 'rb') as f:
            line = f.readline()
            while line != '':
                if 'Frequencies --' in line:
                    freqs = np.append(freqs, [float(frq) for frq in line.split()[2:]])
                line = f.readline()
    else:
        raise ValueError('parse_frequencies() can curtrently only parse QChem and gaussian files,'
                         ' got {0}'.format(software))
    logging.debug('Using parser.parse_frequencies. Determined frequencies are: {0}'.format(freqs))
    return freqs


def parse_t1(path):
    """
    Parse the T1 parameter from a Molpro coupled cluster calculation
    """
    if not os.path.isfile(path):
        raise InputError('Could not find file {0}'.format(path))
    t1 = None
    with open(path, 'rb') as f:
        for line in f:
            if 'T1 diagnostic:' in line:
                t1 = float(line.split()[-1])
    return t1


def parse_e0(path):
    """
    Parse the zero K energy, E0, from an sp job
    """
    if not os.path.isfile(path):
        raise InputError('Could not find file {0}'.format(path))
    log = Log(path='')
    log.determine_qm_software(fullpath=path)
    try:
        e0 = log.loadEnergy(frequencyScaleFactor=1.) * 0.001  # convert to kJ/mol
    except Exception:
        e0 = None
    return e0


def parse_xyz_from_file(path):
    """
    Parse xyz coordinated from:
    .xyz - XYZ file
    .gjf - Gaussian input file
    .out or .log - ESS output file (Gaussian, QChem, Molpro)
    other - Molpro or QChem input file
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    _, file_extension = os.path.splitext(path)

    xyz = None
    relevant_lines = list()

    if file_extension == '.xyz':
        relevant_lines = lines[2:]
    elif file_extension == '.gjf':
        for line in lines[5:]:
            if line and line != '\n' and line != '\r\n':
                relevant_lines.append(line)
            else:
                break
    elif 'out' in file_extension or 'log' in file_extension:
        log = Log(path='')
        log.determine_qm_software(fullpath=path)
        coord, number, mass = log.software_log.loadGeometry()
        xyz = get_xyz_string(xyz=coord, number=number)
    else:
        record = False
        for line in lines:
            if '$end' in line or '}' in line:
                break
            if record and len(line.split()) == 4:
                relevant_lines.append(line)
            elif '$molecule' in line:
                record = True
            elif 'geometry={' in line:
                record = True
        if not relevant_lines:
            raise InputError('Could not parse xyz coordinates from file {0}'.format(path))
    if xyz is None and relevant_lines:
        xyz = ''.join([line for line in relevant_lines if line])
    return xyz

def read_rotors(inputfile):
    """read rotor block of mess input and return rotor dictionary
    'Pivots', 'PES', 'Symmetry', 'Top'
    """
    rotor={}
    rotorline=inputfile.readline()
    while rotorline != '':
        if 'End' in rotorline:
            break
        elif 'Group' in rotorline:
            Top= [int(x) for x in rotorline.split()[1:]]
            rotor['Top']=Top
        elif 'Axis' in rotorline:
            Pivots=[int(x) for x in rotorline.split()[1:]]
            rotor['Pivots']=Pivots
        elif 'Symmetry' in rotorline:
            Symmetry=int(rotorline.split()[1])
            rotor['Symmetry']=Symmetry
        elif 'Potential' in rotorline:
            rotorline=inputfile.readline()
            Points=[float(x) for x in rotorline.split()]
            rotor['PES']=Points
        rotorline=inputfile.readline()
    return rotor

def read_mess(path):
    """read Mess inputfile and return a dictionary containg
     'Freq',
     'ImaginaryFrequency',
     'SymmetryFactor',
     'geom',
     'name',
     'rotors is a array of dictionary with Pivots', 'PotentialPoints', 'Symmetry', 'Top' for each rotor
     """

    data = {}
    data['ImaginaryFrequency'] = []
    rotors = []
    Freq = []
    inp = open(path, 'r')
    line = inp.readline()
    while line != '':

        if 'end barrier' in line:
            break

        elif 'Barrier' in line:
            data['name'] = str(line.split()[5])

        elif 'Geometry' in line:
            geom = []
            line = inp.readline()
            while line != '':
                if 'Core' in line:
                    break
                else:

                    coor = line.split()
                    if len(coor) > 2:
                        geom.append(coor)
                line = inp.readline()
            data['geom'] = geom

        elif 'SymmetryFactor' in line:
            data['SymmetryFactor'] = float(line.split()[1])

        elif 'Frequencies' in line:
            line = inp.readline()
            while line != '':
                if len(line.strip()) == 0:
                    break
                else:

                    vib = [float(x) for x in line.split()]
                    if len(vib) > 0:
                        Freq.append(vib)
                line = inp.readline()
            data['Freq'] = Freq

        elif 'Rotor' in line:

            rotors.append(read_rotors(inp))
            data['rotors'] = rotors
        elif 'ImaginaryFrequency' in line:
            data['ImaginaryFrequency'] = float(line.split()[1])
        line = inp.readline()
        # Close file when finished
    inp.close()
    # data.temperature=temp
    # data.pressure=pressure
    return data