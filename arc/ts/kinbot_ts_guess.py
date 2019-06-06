#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import logging
from arc.arc_exceptions import TSError, ReactionError, InputError
import arc.job.inputs as inputs
from arc.settings import arc_path
from arc.ts.atst import get_reaction_label
# from arc.species.species import _get_possible_conformers_rdkit,get_min_energy_conformer
import json
from kinbot.stationary_pt import StationaryPoint
from kinbot.reaction_generator import ReactionGenerator
from kinbot.qc import QuantumChemistry
from kinbot.reaction_finder import ReactionFinder
# from kinbot.parameters import Parameters
from kinbot.modify_geom import modify_coordinates
from rdkit import Chem
from rdkit.Chem import AllChem
from rmgpy.molecule.molecule import Molecule
from rmgpy.data.rmg import getDB, RMGDatabase


##loading RMG Database, database have to load only once. getrmg_db() will reload the database form memory
# rmg_db = RMGDatabase()
# base_path = '/home/dranasinghe/Software/rmg/RMG-database/input'
# rmg_db.load(path=base_path,kineticsFamilies='all')
# loading sample json. we will chage parameters later.
class KinbotPara:
    # this is to make sure we ave all the parameters kinbot need. if in future kinbot update parameters
    # we need to update the input_file dictionary
    def __init__(self):
        KinbotJson = inputs.input_files['kinbot']
        self.par = json.loads(KinbotJson)


par = KinbotPara()

# reaction supported by kibot
kinbot_rxn_types = ["intra_H_migration", "intra_H_migration_suprafacial", "intra_R_migration", "intra_OH_migration",
                    "cpd_H_migration", "Intra_RH_Add_Endocyclic_F", "Intra_RH_Add_Endocyclic_R",
                    "Cyclic_Ether_Formation", "Intra_RH_Add_Exocyclic_F", "Intra_RH_Add_Exocyclic_R", "Retro_Ene",
                    "Intra_R_Add_Endocyclic_F", "Intra_R_Add_ExoTetCyclic_F", "Intra_R_Add_Exocyclic_F", "Korcek_step2",
                    "r22_cycloaddition", "r12_cycloaddition", "r12_insertion_R", "r13_insertion_CO2",
                    "r13_insertion_ROR", "Diels_alder_addition", "Intra_Diels_alder_R", "ketoenol",
                    "HO2_Elimination_from_PeroxyRadical", "R_Addition_COm3_R", "R_Addition_MultipleBond",
                    "12_shift_S_F", "12_shift_S_R", "R_Addition_CSm_R", "r13_insertion_RSR"]


# TODO rewrite to accep RMG reactions
def ts_guess(reaction_label=None, rmg_reaction=None, reaction_family=None):
    """provide TS guess as xyz
    :param reaction_label: 
    :param rmg_reaction: 
    :param reaction_family: 
    :return: 
    """
    CreateKinbotInput = False
    xyz = None
    rxnReacIndexfamily = reaction_family

    if rmg_reaction is not None and reaction_label is None:
        reaction_label = get_reaction_label(rmg_reaction)

    elif reaction_label is None:
        raise TSError('Must get either reaction_label or rmg_reaction')
    # rxnReacIndexfamily, rReacIndex, pReacIndex = reactant_and_family(data, ReacIndex)
    reactants = rmg_reaction.reactants
    products = rmg_reaction.products
    if len([reactants]) > 1:
        raise TSError("Only unimoleculre reactions as reactant is possible")
    rReacIndex = reactants.molecule[0]
    pReacIndex = [product.molecule[0] for product in products]
    print('{0} <-> {1} {2}'.format(rReacIndex, pReacIndex, rxnReacIndexfamily))

    TSinfo = from_kinbot_get_TS_guess(rReacIndex, pReacIndex, rxnReacIndexfamily, reaction_label, par,
                                      CreateKinbotInput)
    xyz = TSinfo[2]
    return xyz


def from_kinbot_get_TS_guess(rReacIndex=None, pReacIndex=None, rxnReacIndexfamily=None, reaction_label=None, par=par,
                             CreateKinbotInput=False):
    """for given rection index find the TS guess
      rReacIndex                reactant RMG species
      pReacIndex                product RMG species
      rxnReacIndexfamily        reaction family
      par                       kinbot parameters as an object
      CreateKinbotInput         Flag for creating kinbot
      :rtype: object
     """
    if rReacIndex is None:
        raise TSError('reactant is mising')
    elif len([rReacIndex]) > 1:
        raise TSError('only unimolecular reactions are possible')

    if pReacIndex is None:
        raise TSError('producnts are mising')
    if not isinstance(pReacIndex, list):
        raise TSError('provide producs as an array')

    if rxnReacIndexfamily is None:
        raise TSError('need a reaction lable')

    if rxnReacIndexfamily not in kinbot_rxn_types:
        raise TSError('reaction type not supported by Kinbot')
    if reaction_label is None:
        raise TSError('need a reaction lable')

    kb_objReacIndex = setup_kinbot(rReacIndex, rxnReacIndexfamily, par, CreateKinbotInput, reaction_label)

    try:
        rReacIndexres, pReacIndexres = generate_resonance_structure(r=rReacIndex, p=pReacIndex)
        print(' Reaction : {0} kinbot instants list: {1}'.format(rxnReacIndexfamily, kb_objReacIndex.species.reac_inst))

        num = go_through_all_resonance(rReacIndexres, pReacIndexres, rxnReacIndexfamily, kb_objReacIndex)
        return [num,kb_objReacIndex.species.reac_inst[num],kb_objReacIndex.species.reac_name[num],get_ts_guess(kb_objReacIndex, num, check=True)]
    except:
        raise TSError('Failed to find a TS guess using Kinbot')
        pass


def FromJsonlistGetTSGuess(data, ReacIndex, CreateKinbotInput):
    """for given rection index find the TS guess
    data from json
    ReacIndex reaction index
    CreateKinbotInput boolean generate Kinbot input
    """
    rxnReacIndexfamily, rReacIndex, pReacIndex = reactant_and_family(data, ReacIndex)
    print('{0} <-> {1} {2}'.format(rReacIndex, pReacIndex, rxnReacIndexfamily))
    kb_objReacIndex = setup_kinbot(rReacIndex, rxnReacIndexfamily, par, CreateKinbotInput, ReacIndex)

    try:
        rReacIndexres, pReacIndexres = generate_resonance_structure(rReacIndex, pReacIndex)
        print(' Reaction : {0} kinbot instants list: {1}'.format(rxnReacIndexfamily, kb_objReacIndex.species.reac_inst))

        num = go_through_all_resonance(rReacIndexres, pReacIndexres, rxnReacIndexfamily, kb_objReacIndex)
        print_structure(get_ts_guess(kb_objReacIndex, num, check=True))
        return [num, kb_objReacIndex.species.reac_name[num], get_ts_guess(kb_objReacIndex, num, check=True)]
    except:
        raise TSError('Failed to find a TS guess using Kinbot')
        pass


def print_structure(sp):
    """Print the geometry as strings"""
    for line in sp:
        print('{0} {1} {2} {3}'.format(str(line[0]), str(line[1]), str(line[2]), str(line[3])))


def combine_atom_geom(atoms, geoms):
    """comnine atom labels with cartiasin coordinates"""
    struc = []
    for i in range(len(atoms)):
        struc.append([atoms[i], geoms[i][0], geoms[i][1], geoms[i][2]])
    return struc


def getPossibleConformersRDKitUFF(mol):
    """
     use only for molecules.atosm < 3 use FF optimize geomatry.
    """
    # print "using UFFOptimizeMolecule"
    rdmol, rdInds = mol.toRDKitMol(removeHs=False, returnMapping=True)
    rdIndMap = dict()
    for k, atm in enumerate(mol.atoms):
        ind = rdInds[atm]
        rdIndMap[ind] = k
    if len(mol.atoms) < 3:
        AllChem.EmbedMultipleConfs(rdmol, numConfs=1, randomSeed=1, )
    else:
        AllChem.EmbedMultipleConfs(rdmol, numConfs=(len(mol.atoms) - 2) * 6, randomSeed=1, )
    # AllChem.EmbedMultipleConfs(rdmol,numConfs=(len(mol.atoms)-3)*60,randomSeed=1,)
    energies = []
    xyzs = []
    for i in xrange(rdmol.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(rdmol, confId=i)
        E = AllChem.UFFGetMoleculeForceField(rdmol, confId=i).CalcEnergy()

        energies.append(E)
        cf = rdmol.GetConformer(i)
        xyz = []
        for j in xrange(cf.GetNumAtoms()):
            pt = cf.GetAtomPosition(j)
            xyz.append([pt.x, pt.y, pt.z])
        xyz = [xyz[rdIndMap[i]] for i in xrange(len(xyz))]
        xyzs.append(xyz)

    return xyzs, energies


def getMinEnergyConformer(xyzs, energies):
    """get the minimum conformer"""
    minval = min(energies)
    minind = energies.index(minval)
    return xyzs[minind]


def get_strucuture(spc):
    """get minimum energy confomformer structure"""
    xyzs, E = getPossibleConformersRDKitUFF(spc)
    xyzmin = getMinEnergyConformer(xyzs, E)
    geom = []
    atomList = [x.element.symbol for x in spc.atoms]
    for i in range(len(xyzmin)):
        geom.append(str(atomList[i]))
        geom.append(str(xyzmin[i][0]))
        geom.append(str(xyzmin[i][1]))
        geom.append(str(xyzmin[i][2]))
    return geom


def setup_kinbot(spc=None, rxnfamily=None, Par=None, CreateKinbotInput=False, ReacIndex="ARC"):
    """spc           rmg molecule object
        rxnfamily    string one of reaction families supported by kinbot
        Par          parameters for json creat by kinbot tomake sure I did miss any parameters
        CreateKinbotInput Flag to create Kibot input
        generate kinbot reaction object for specific reaction.rg object conatin all the information need.eg.
    # rg.species.reac_type have reaction type
    # rg.species.reac_name reaction name and file names for kinbot
    # rg.species.reac_inst atom index list
    :type CreateKinbotInput: bool
            """
    if spc is None:
        TSError('need a reactant object')
    if not isinstance(spc, Molecule):
        TSError('need a reactant as a RMG species object')
    if Par is None:
        TSError('Kinbot parameter object is not available')
    if rxnfamily is None:
        TSError('require reaction family')

    Par.par["smiles"] = spc.toSMILES()
    Par.par['charge'] = 0  # RMG only support neutral species
    Par.par['mult'] = spc.multiplicity
    Par.par['families'] = rxnfamily
    Par.par['structure'] = get_strucuture(spc)

    # kinbot generate Well0 and it properties
    well0 = StationaryPoint('well0',
                            Par.par['charge'],
                            Par.par['mult'],
                            # smiles=Par.Par['smiles'],
                            structure=Par.par['structure'],
                            )

    well0.calc_chemid()
    well0.bond_mx()
    well0.find_cycle()
    well0.find_atom_eqv()
    well0.find_conf_dihedral()

    # nessary for determine in level of theroy. have to have this to find reactions
    qc = QuantumChemistry(Par)
    rf = ReactionFinder(well0, Par, qc)
    rf.find_reactions()
    # generate reactions
    rg = ReactionGenerator(well0, Par, qc)  # rg is reaction object
    # rg.species.reac_type have reaction type
    # rg.species.reac_name reaction name and file names for kinbot
    # rg.species.reac_inst atom index list

    # generating a kinbot inputfile
    if CreateKinbotInput:
        f = open(str(arc_path) + str(ReacIndex) + "_" + str(rxnfamily) + ".json", 'w')
        json.dump(Par.par, f, separators=(',\n', ": "))
        f.close()

    return rg


def get_ts_guess(rxn=None, index=None, check=True):
    """rxn   Kinbot reaction object
       index  which index
       check   flag to print succcess Kinbot determine if TS guess is reasonable """
    if rxn is None:
        TSError('No kinbot reaction generator object was given')

    if index is None or not isinstance(index, int):
        TSError('No kinbot reaction index was given, please provide a integer')

    print('kinbot index guess close to RMG:{0}'.format(index))
    obj = rxn.species.reac_obj[index]
    instance_name = obj.instance_name
    geom_well0 = obj.species.geom

    step, fix, change, release = obj.get_constraints(12, geom_well0)

    change_starting_zero = []
    for c in change:
        c_new = [ci - 1 for ci in c[:-1]]
        c_new.append(c[-1])
        change_starting_zero.append(c_new)
    if len(change_starting_zero) > 0:
        success, geom_ts = modify_coordinates(obj.species, obj.instance_name, geom_well0, change_starting_zero,
                                              obj.species.bond)
        for c in change:
            fix.append(c[:-1])
        change = []

    if check and success == 1:
        print('succesfully found a TS guess')
        return [obj.species.atom, geom_ts]
    if not check:
        print('TS guess from Kinbot failed internal testing')
        return [obj.species.atom, geom_ts]


def match_with_rmg(react=None, prod=None, rxnfamily=None, kinbot_output=None, rmg_db=None):
    """This will match which Kinbot guess match with RMG reaction"""
    if react is None:
        raise TSError('reactant is mising')
    elif len([react]) > 1:
        raise TSError('only unimolecular reactions are possible')

    if prod is None:
        raise TSError('producnts are mising')
    if not isinstance(prod, list):
        raise TSError('provide producs as an array')

    if rxnfamily is None:
        raise TSError('need a reaction lable')
    if rmg_db is None:
        rmg_db = getDB()

    print('kinbot atom lists:{0}'.format(kinbot_output))

    family = rmg_db.kinetics.families[rxnfamily]
    reactant_structures, product_structures = family.getLabeledReactantsAndProducts(reactants=react, products=prod)
    reactant_structure = reactant_structures[0]

    atom_index_dict = dict()
    for i, atom in enumerate(reactant_structure.atoms):
        atom_index_dict[atom] = i  # starting from 0
    # print atom_index_dict

    labeled_atoms_dict = reactant_structure.getLabeledAtoms()

    ordered_list_of_labeled_atoms = list()
    for i in range(len(labeled_atoms_dict)):
        atom = labeled_atoms_dict['*{0}'.format(i + 1)]
        ordered_list_of_labeled_atoms.append(atom_index_dict[atom])

    print('rmg atom list:{0}'.format(ordered_list_of_labeled_atoms))

    # print 'in the kinbot output this will be TS number...:'
    # special case for Cyclic_Ether_Formation is needed RMG template is longer than kinbot
    if str(rxnfamily) == "Cyclic_Ether_Formation":
        sorted_kinbot_output = [sorted(x) for x in kinbot_output]
        sorted_atom_list = sorted(ordered_list_of_labeled_atoms)
        for A in sorted_kinbot_output:
            C = list(set(A) & set(sorted_atom_list))
            if len(C) == len(A):
                return sorted_kinbot_output.index(C)
    elif ordered_list_of_labeled_atoms in kinbot_output:

        print('exact match found')
        kinbot_output.index(ordered_list_of_labeled_atoms)
        return kinbot_output.index(ordered_list_of_labeled_atoms)

    else:

        print("no exact match found. Sorting arrays to find a match.....")
        sorted_kinbot_output = [sorted(x) for x in kinbot_output]
        sorted_atom_list = sorted(ordered_list_of_labeled_atoms)
        if sorted_atom_list in sorted_kinbot_output:
            print("found a match after sorting")

            sorted_kinbot_output.index(sorted_atom_list)
            return sorted_kinbot_output.index(sorted_atom_list)


# TODO remove this section
def reactant_and_family(RXNS, index):
    """obtain reaction info relataed to index from RXNS json file"""
    f = str(RXNS[index]['Family'])
    r = Molecule(SMILES=str(RXNS[index]['Reactants'][0]["SMILES"][0]))
    p = []
    for i in range(len(RXNS[index]['Products'])):
        p.append(Molecule(SMILES=str(RXNS[index]['Products'][i]["SMILES"][0])))
    return f, r, p


# generare resonance structures
def generate_resonance_structure(r=None, p=None):
    """ generate resonace for reactants and products
    :param r: 
    :param p: 
    :return: 
    """
    rRes = []
    if r is None:
        raise TSError("need reactant as species")
    if p is None:
        raise TSError("need product as species")
    if not isinstance(p, list):
        raise TSError('provide producs as an array')

    rtemp = r.generate_resonance_structures()
    for i in rtemp:
        rRes.append([i])
    print("reactant resonace structures {0}".format(rRes))
    pRes = []
    iRes = []
    for i in p:
        iRes.append(i.generate_resonance_structures())
    if len(p) == 1:
        # print "one product"
        for i in iRes[0]:
            pRes.append([i])
    if len(p) == 2:
        # print "two product"
        for i in iRes[0]:
            for j in iRes[1]:
                pRes.append([i, j])
    print("reactant resonace structures {0}".format(pRes))

    return rRes, pRes


def go_through_all_resonance(rRes, pRes, rxnfamily, kb_obj):
    """go through all the resonce stucures and reactions and find whicj one is a match
    :return: 
    :param rRes: 
    :param pRes: 
    :param rxnfamily: 
    :param kb_obj: 
    :return: 
    """
    Indexes = []
    for r in rRes:
        for p in pRes:
            print('reactant: {0} and products:{1}'.format(r, p))
            try:
                # Matched_index=match_with_rmg(r,p,rxnfamily,kb_obj.species.reac_inst)
                Indexes.append(
                    match_with_rmg(react=r, prod=p, rxnfamily=rxnfamily, kinbot_output=kb_obj.species.reac_inst))
            except:
                print("going through other resonace combinations")
                pass
    # return Matched_index
    return Indexes[0]
