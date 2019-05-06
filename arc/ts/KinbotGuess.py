#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import logging
from arc.arc_exceptions import TSError, ReactionError, InputError
import arc.job.inputs as inputs
from arc.settings import arc_path
from arc.ts.atst import get_reaction_label
#from arc.species.species import _get_possible_conformers_rdkit,get_min_energy_conformer
import json
from kinbot.stationary_pt import StationaryPoint
from kinbot.reaction_generator import ReactionGenerator
from kinbot.qc import QuantumChemistry
from kinbot.reaction_finder import ReactionFinder
#from kinbot.parameters import Parameters
from kinbot.modify_geom import modify_coordinates
from rdkit import Chem
from rdkit.Chem import AllChem
from rmgpy.molecule.molecule import Molecule
from rmgpy.data.rmg import getDB, RMGDatabase

##loading RMG Database, database have to load only once. getDB() will reload the database form memory
#db = RMGDatabase()
#base_path = '/home/dranasinghe/Software/rmg/RMG-database/input'
#db.load(path=base_path,kineticsFamilies='all')
#loaing sample json. we will chage parameters later.
class KinbotPara:
    def __init__(self):
        KinbotJson=inputs.input_files['kinbot']
        self.par=json.loads(KinbotJson)

par=KinbotPara()

#par = Parameters("sample.json")
#CreateKinbotInput=False #Flag fo r creating Kinbot input

#reaction supported by kibot
kinbot_rxn_types=["intra_H_migration","intra_H_migration_suprafacial","intra_R_migration","intra_OH_migration","cpd_H_migration","Intra_RH_Add_Endocyclic_F","Intra_RH_Add_Endocyclic_R","Cyclic_Ether_Formation","Intra_RH_Add_Exocyclic_F","Intra_RH_Add_Exocyclic_R","Retro_Ene","Intra_R_Add_Endocyclic_F","Intra_R_Add_ExoTetCyclic_F","Intra_R_Add_Exocyclic_F","Korcek_step2","r22_cycloaddition","r12_cycloaddition","r12_insertion_R","r13_insertion_CO2","r13_insertion_ROR","Diels_alder_addition","Intra_Diels_alder_R","ketoenol","HO2_Elimination_from_PeroxyRadical","R_Addition_COm3_R","R_Addition_MultipleBond","12_shift_S_F","12_shift_S_R","R_Addition_CSm_R","r13_insertion_RSR"]
#TODO rewrite to accep RMG reactions
def TSGuess(reaction_label=None, rmg_reaction=None, reaction_family=None):

    xyz = None
    rxnReacIndexfamily=reaction_family
    ReacIndex=reaction_label
    if rmg_reaction is not None and reaction_label is None:
        reaction_label = get_reaction_label(rmg_reaction)
        ReacIndex = reaction_label
    elif reaction_label is None:
        raise TSError('Must get either reaction_label or rmg_reaction')
    #rxnReacIndexfamily, rReacIndex, pReacIndex = ReactantandFamily(data, ReacIndex)
    reactants = rmg_reaction.reactants
    products = rmg_reaction.products
    if len([reactants])>1:
        raise TSError("Only unimoleculre reactions as reactant is possible")
    rReacIndex=reactants.molecule[0]
    pReacIndex= [product.molecule[0] for product in products]
    logging.info('{0} <-> {1} {2}'.format(rReacIndex, pReacIndex, rxnReacIndexfamily))
    TSinfo=FromKinbotGetTSGuess(rReacIndex, pReacIndex, rxnReacIndexfamily,reaction_label, par, CreateKinbotInput)
    xyz=TSinfo[2]
    return xyz


def  FromKinbotGetTSGuess(rReacIndex=None, pReacIndex=None, rxnReacIndexfamily=None, reaction_label= None, par=par, CreateKinbotInput=False):
    """for given rection index find the TS guess
      rReacIndex                reactant RMG species
      pReacIndex                product RMG species
      rxnReacIndexfamily        reaction family
      par                       data form json
      CreateKinbotInput         Flag for creating kinbot
     """
    if rReacIndex is None:
        raise TSError('reactant is mising')
    elif len([rReacIndex]) > 1:
        raise TSError('only unimolecular reactions are possible')

    if pReacIndex is None:
        raise TSError('producnts are mising')
    if rxnReacIndexfamily is None:
        raise TSError('need a reaction lable')

    if rxnReacIndexfamily not in kinbot_rxn_types:
        raise TSError('reaction type not supported by Kinbot')
    if reaction_label is None:
        raise TSError('need a reaction lable')


    kb_objReacIndex = SetupKinbot(rReacIndex, rxnReacIndexfamily, par, CreateKinbotInput,reaction_label)

    try:
        rReacIndexres, pReacIndexres = GenerateResonaceStructure(rReacIndex, pReacIndex)
        logging.info(' Reaction : {0} kinbot instants list: {1}'.format(rxnReacIndexfamily, kb_objReacIndex.species.reac_inst))

        num = GothourghAllresonace(rReacIndexres, pReacIndexres, rxnReacIndexfamily, kb_objReacIndex)
        PrintStrucuture(GetTSGuess(kb_objReacIndex, num, check=True))
        return [num, kb_objReacIndex.species.reac_name[num],GetTSGuess(kb_objReacIndex, num, check=True)]
    except:
        raise TSError('Failed to find a TS guess using Kinbot')
        pass

def FromJsonlistGetTSGuess(data, ReacIndex, CreateKinbotInput):
    """for given rection index find the TS guess
    data from json
    ReacIndex reaction index
    CreateKinbotInput boolean generate Kinbot input
    """
    rxnReacIndexfamily, rReacIndex, pReacIndex = ReactantandFamily(data, ReacIndex)
    logging.info('{0} <-> {1} {2}'.format(rReacIndex, pReacIndex, rxnReacIndexfamily))
    kb_objReacIndex = SetupKinbot(rReacIndex, rxnReacIndexfamily, par, CreateKinbotInput, ReacIndex)

    try:
        rReacIndexres, pReacIndexres = GenerateResonaceStructure(rReacIndex, pReacIndex)
        logging.info(' Reaction : {0} kinbot instants list: {1}'.format(rxnReacIndexfamily, kb_objReacIndex.species.reac_inst))

        num = GothourghAllresonace(rReacIndexres, pReacIndexres, rxnReacIndexfamily, kb_objReacIndex)
        PrintStrucuture(GetTSGuess(kb_objReacIndex, num, check=True))
        return [num, kb_objReacIndex.species.reac_name[num],GetTSGuess(kb_objReacIndex, num, check=True)]
    except:
        raise TSError('Failed to find a TS guess using Kinbot')
        pass

def PrintStrucuture(sp):
    """Print the geometry as strings"""
    for line in sp:
         logging.info('{0} {1} {2} {3}'.format(str(line[0]),str(line[1]),str(line[2]), str(line[3])))

def combineatomgeom(atoms,geoms):
    """comnine atom labels with cartiasin coordinates"""
    struc=[]
    for i in range(len(atoms)):
        struc.append([atoms[i],geoms[i][0],geoms[i][1],geoms[i][2]])
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

def getMinEnergyConformer(xyzs,energies):
    """get the minimum conformer"""
    minval = min(energies)
    minind = energies.index(minval)
    return xyzs[minind]

def GetStrucuture(spc):
    """get minimum energy confomformer structure"""
    xyzs, E=getPossibleConformersRDKitUFF(spc)
    xyzmin=getMinEnergyConformer(xyzs,E)
    geom=[]
    atomList = [x.element.symbol for x in spc.atoms]
    for i in range(len(xyzmin)):
            geom.append(str(atomList[i]))
            geom.append(str(xyzmin[i][0]))
            geom.append(str(xyzmin[i][1]))
            geom.append(str(xyzmin[i][2]))
    return geom


def SetupKinbot(spc, rxnfamily, par, CreateKinbotInput=False, ReacIndex="ARC"):
    """spc           rmg molecule object
        rxnfamily    string one of reaction families supported by kinbot
        par          parameters for json creat by kinbot tomake sure I did miss any parameters
        CreateKinbotInput Flag to create Kibot input
        generate kinbot reaction object for specific reaction.rg object conatin all the information need.eg.
    # rg.species.reac_type have reaction type
    # rg.species.reac_name reaction name and file names for kinbot
    # rg.species.reac_inst atom index list
            """

    par.par["smiles"] = spc.toSMILES()
    par.par['charge'] = 0  # RMG only support neutral species
    par.par['mult'] = spc.multiplicity
    par.par['families'] = rxnfamily
    par.par['structure'] = GetStrucuture(spc)

    # kinbot generate Well0 and it properties
    well0 = StationaryPoint('well0',
                            par.par['charge'],
                            par.par['mult'],
                            # smiles=par.par['smiles'],
                            structure=par.par['structure'],
                            )

    well0.calc_chemid()
    well0.bond_mx()
    well0.find_cycle()
    well0.find_atom_eqv()
    well0.find_conf_dihedral()

    # nessary for determine in level of theroy. have to have this to find reactions
    qc = QuantumChemistry(par)
    rf = ReactionFinder(well0, par, qc)
    rf.find_reactions()
    # generate reactions
    rg = ReactionGenerator(well0, par, qc)  # rg is reaction object
    # rg.species.reac_type have reaction type
    # rg.species.reac_name reaction name and file names for kinbot
    # rg.species.reac_inst atom index list

    if CreateKinbotInput:
        f = open(str(ReacIndex) + "_" + str(rxnfamily) + ".json", 'w')
        json.dump(par.par, f, separators=(',\n', ": "))
        f.close()

    return rg


def GetTSGuess(rxn, index, check):
    """rxn   Kinbot reaction object
       index  which index
       check   flag to print succcess Kinbot determine if TS guess is reasonable """

    logging.info('kinbot index guess close to RMG:{0}'.format(index))
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
        logging.info('succesfully found a TS guess')
        return combineatomgeom(obj.species.atom, geom_ts)
    if not check:
        logging.warning('TS guess from Kinbot failed internal testing')
        return combineatomgeom(obj.species.atom, geom_ts)


def MatchWithRMG(react, prod, rxnfamily, kinbot_output):
    """This will match which Kinbot guess match with RMG reaction"""
    db = getDB()

    logging.info('kinbot atom lists:{0}'.format(kinbot_output))

    family = db.kinetics.families[rxnfamily]
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

    logging.info('rmg atom list:{0}'.format(ordered_list_of_labeled_atoms))

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

        logging.info('found exact match')
        kinbot_output.index(ordered_list_of_labeled_atoms)
        return kinbot_output.index(ordered_list_of_labeled_atoms)

    else:

        logging.info("no exact match found. Sorting arrays to find a match.....")
        sorted_kinbot_output = [sorted(x) for x in kinbot_output]
        sorted_atom_list = sorted(ordered_list_of_labeled_atoms)
        if sorted_atom_list in sorted_kinbot_output:

            logging.info("found a match after sorting")

            sorted_kinbot_output.index(sorted_atom_list)
            return sorted_kinbot_output.index(sorted_atom_list)
#TODO remove this section
def ReactantandFamily(RXNS,index):
    """obtain reaction info relataed to index from RXNS json file"""
    f= str(RXNS[index]['Family'])
    r=Molecule(SMILES=str(RXNS[index]['Reactants'][0]["SMILES"][0]))
    p=[]
    for i in range(len(RXNS[index]['Products'])):
        p.append(Molecule(SMILES=str(RXNS[index]['Products'][i]["SMILES"][0])))
    return f,r,p


# generare resonance structures
def GenerateResonaceStructure(r, p):
    """ generate resonace for reactants and products"""
    rRes = []
    rtemp = r.generate_resonance_structures()
    for i in rtemp:
        rRes.append([i])
    # print rRes
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
    # print (pRes)

    return rRes, pRes

def GothourghAllresonace(rRes,pRes,rxnfamily,kb_obj):
    """go through all the resonce stucures and reactions and find whicj one is a match"""
    Indexes=[]
    for r in rRes:
        for p in pRes:
            logging.info('reactant: {0} and products:{1}'.format(r,p))
            try:
                #Matched_index=MatchWithRMG(r,p,rxnfamily,kb_obj.species.reac_inst)
                Indexes.append(MatchWithRMG(r,p,rxnfamily,kb_obj.species.reac_inst))
            except:
                logging.info("going through other resonace combinations")
                pass
    #return Matched_index
    return Indexes[0]




