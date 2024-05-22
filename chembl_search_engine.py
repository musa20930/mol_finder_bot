from loguru import logger
import json
import csv
from io import StringIO, BytesIO

from aiogram.types import input_file

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.utils import utils
from validation_models import Smiles
from pydantic import ValidationError


async def retrieve_by_name(name: str) -> list:
    logger.debug(f'Searching \'{name}\' by name')
    molecule = new_client.molecule
    result = molecule.filter(
        molecule_synonyms__molecule_synonym__iexact=name).only(['molecule_chembl_id', 
                                                                'pref_name', 
                                                                'molecule_structures'])
    logger.debug(f'Search by name {result = }')
    return result


async def retrieve_by_id(chembl_id: int) -> list:
    logger.debug(f'Searching by id \'{chembl_id}\'')
    molecule = new_client.molecule
    search_text = f'{chembl_id}' if chembl_id.casefold().startswith('chembl') else f'CHEMBL{chembl_id}'
    result = molecule.filter(chembl_id=search_text).only(['molecule_chembl_id', 
                                                          'pref_name', 
                                                          'molecule_structures']) 
    logger.debug(f'Search by id {result = }')
    return result


async def retrieve_by_inchi(inchi: str) -> list:
    logger.debug(f'Searching \'{inchi}\' by INCHI key')
    molecule = new_client.molecule
    result = molecule.filter(
        molecule_structures__standard_inchi_key=inchi).only(['molecule_chembl_id', 
                                                            'pref_name', 
                                                            'molecule_structures'])
    logger.debug(f'Search by INCHI key {result = }')
    return result


async def retrieve_by_smiles(smiles: str) -> list:
    logger.debug(f'Validating SMILES: \'{smiles}\'')
    try:
        s = Smiles(smiles=smiles)
        mol_smiles = s.smiles
    except ValidationError as e:
        logger.error(e.errors())
        return None
    except Exception as e:
        logger.error(e)
        return None

    logger.debug(f'Searching \'{smiles}\' by SMILES')
    molecule = new_client.molecule
    result = molecule.filter(smiles=mol_smiles).only(['molecule_chembl_id', 
                                                      'pref_name', 
                                                      'molecule_structures'])
    logger.debug(f'Search by SMILES {result = }')
    return result


async def retrieve_similar_molecules(mol_smiles: str, similarity_percent: int) -> list:
    logger.debug(f"Looking for molecules similar to \'{mol_smiles}\' by {similarity_percent}%")
    similarity = new_client.similarity
    # Similarity score must be between 40 and 100
    result = similarity.filter(smiles=mol_smiles, similarity=similarity_percent).only(['molecule_chembl_id', 
                                                                                       'similarity',
                                                                                       'molecule_structures'])
    return result


async def retrieve_properties(mol_smiles: str) -> str:
    mol = utils.smiles2ctab(mol_smiles)
    descs = json.loads(utils.chemblDescriptors(mol))
    return descs


async def save_properties(descs, mol_smiles: str, chembl_id):
    for desc in descs:
        desc['ChEMBLID'] = chembl_id
        desc['SMILES'] = mol_smiles
    file = StringIO()
    columns = [
        'ChEMBLID', 'SMILES', 'qed', 'MolWt', 'TPSA', 'HeavyAtomCount', 
        'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds', 
        'MolLogP', 'MolecularFormula', 'Ro3Pass', 'NumRo5', 'MonoisotopicMolWt'
        ]
    dict_writer = csv.DictWriter(file, columns, delimiter=';')
    dict_writer.writeheader()
    dict_writer.writerows(descs)
    file.seek(0)
    sdf_file = input_file.BufferedInputFile(BytesIO(file.read().encode('utf-8')).getbuffer(), 
                                            filename=f"{chembl_id}.csv")
    return sdf_file


async def standardize_mol(mol_smiles: str, chembl_id):
    logger.debug(f'Standardizing \'{mol_smiles}\'')
    mol = utils.smiles2ctab(mol_smiles)
    st = json.loads(utils.standardize(mol))
    standard_molblock = st[0]['standard_molblock']
    file = StringIO()
    file.write(standard_molblock)
    file.seek(0)
    sdf_file = input_file.BufferedInputFile(BytesIO(file.read().encode('utf-8')).getbuffer(), 
                                            filename=f"{chembl_id}_standard.sdf")
    return sdf_file


async def retrieve_by_mass(mol_mass: float) -> list:
    logger.debug(f'Searching by mass \'{mol_mass}\'')
    molecule = new_client.molecule
    mols = molecule.filter(molecule_properties__mw_freebase__lte=mol_mass).only(['molecule_chembl_id', 
                                                                                 'molecule_structures',
                                                                                 'pref_name'])
    return mols


async def retrieve_by_(mol_mass: float) -> list:
    logger.debug(f'Searching by mass \'{mol_mass}\'')
    molecule = new_client.molecule
    mols = molecule.filter(molecule_properties__mw_freebase__lte=mol_mass).only(['molecule_chembl_id', 
                                                                                 'molecule_structures',
                                                                                 'pref_name'])
    return mols

