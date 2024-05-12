from loguru import logger

from aiogram import types, Router
from aiogram.fsm.context import FSMContext
from aiogram.filters.command import Command
from aiogram import F

from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw

from pydantic import ValidationError
from validation_models import SmilesSimilarity, Smiles
from states import SearchInfo


search_router = Router()


@search_router.message(Command('search'))
async def search(message: types.Message, state: FSMContext) -> None:
    await state.set_state(SearchInfo.search_type)
    keyboard = [
        [types.KeyboardButton(text='By Name'), 
         types.KeyboardButton(text='By Chembl ID'),
         types.KeyboardButton(text='By Similarity')],
    ]
    reply_markup = types.ReplyKeyboardMarkup(keyboard=keyboard, 
                                             resize_keyboard=True, 
                                             input_field_placeholder="Choose option")
    await message.answer(
        text="Choose search type", 
        reply_markup=reply_markup
    )


@search_router.message(SearchInfo.search_type, F.text.casefold() == 'by name')
@search_router.message(SearchInfo.molecule_next_step, F.text.casefold() == 'search by name')
async def search_by_name(message: types.Message, state: FSMContext) -> None:
    await state.update_data(search_type=message.text)
    await state.set_state(SearchInfo.molecule_name)
    await message.reply(
        'Enter name of molecule', 
        reply_markup=types.ReplyKeyboardRemove()
    )


@search_router.message(SearchInfo.search_type, F.text.casefold() == 'by chembl id')
@search_router.message(SearchInfo.molecule_next_step, F.text.casefold() == 'search by chembl id')
async def search_by_id(message: types.Message, state: FSMContext) -> None:
    if 'search' not in message.text:
        await state.update_data(search_type=message.text)
    else:
        await state.update_data(molecule_next_step=message.text)
    await state.set_state(SearchInfo.chembl_id)
    await message.reply(
        'Enter id', 
        reply_markup=types.ReplyKeyboardRemove()
    )


@search_router.message(SearchInfo.molecule_name)
async def get_mol_by_name(message: types.Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    await state.set_state(SearchInfo.molecule_next_step)
    name = message.text
    logger.debug(f'Searching \'{name}\' by name')
    molecule = new_client.molecule
    result = molecule.filter(
        molecule_synonyms__molecule_synonym__iexact=name).only(['molecule_chembl_id', 
                                                                'pref_name', 
                                                                'molecule_structures'])
    logger.debug(f'Search by name {result = }')
    if result:
        pref_name = result[0]['pref_name']
        await message.answer(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        keyboard = [
            [types.KeyboardButton(text='Search by Similarity❄️❄️'), 
            types.KeyboardButton(text='Get Molecule Structure🧬'),
            types.KeyboardButton(text='Search by ChEMBL ID')],
        ]
        reply_markup = types.ReplyKeyboardMarkup(keyboard=keyboard, 
                                                resize_keyboard=True, 
                                                input_field_placeholder="Choose option")
        await message.answer(
            text="You can now search for new molecule by name "
            "or rerun /search command to search by ChEMBL ID", 
            reply_markup=reply_markup
        )
    else:
        await state.set_state(SearchInfo.molecule_name)
        await message.answer(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.chembl_id)
async def get_mol_by_id(message: types.Message, state: FSMContext) -> None:
    await state.update_data(chembl_id=message.text)
    await state.set_state(SearchInfo.molecule_next_step)
    chembl_id = message.text
    logger.debug(f'Searching by id \'{chembl_id}\'')
    molecule = new_client.molecule
    result = molecule.filter(molecule_chembl_id=f'CHEMBL{chembl_id}').only(['molecule_chembl_id', 
                                                                            'pref_name', 
                                                                            'molecule_structures']) 
    logger.debug(f'Search by id {result = }')
    if result:
        pref_name = result[0]['pref_name']
        await message.reply(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        keyboard = [
            [types.KeyboardButton(text='Search by Similarity❄️❄️'), 
            types.KeyboardButton(text='Get Molecule Structure🧬'), 
            types.KeyboardButton(text='Search by name')], 
        ]
        reply_markup = types.ReplyKeyboardMarkup(keyboard=keyboard, 
                                                resize_keyboard=True, 
                                                input_field_placeholder="Choose option")
        await message.answer(
            text="You can now search for new molecule by ChEMBL ID "
            "or rerun /search command to search by name", 
            reply_markup=reply_markup
        )
    else:
        await state.set_state(SearchInfo.molecule_name)
        await message.reply(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.molecule_next_step, F.text.casefold().startswith('get molecule structure'))
async def get_mol_structure(message: types.Message, state: FSMContext) -> None:
    await state.update_data(molecule_next_step=message.text)
    await state.set_state(SearchInfo.molecule_next_step)
    data = await state.get_data()
    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    logger.debug(f'Drawing structure of \'{mol_smiles}\'')
    
    pl = Chem.MolFromSmiles(mol_smiles)
    pl_image = Draw.MolToImage(pl, size=(400, 400))
    pl_image.save('assets/image.png')

    mol_img = types.FSInputFile('assets/image.png')
    await message.answer_photo(photo=mol_img)


@search_router.message(SearchInfo.search_type, F.text.casefold() == 'by similarity')
@search_router.message(SearchInfo.molecule_next_step, F.text.casefold().startswith('search by similarity'))
async def search_by_similarity(message: types.Message, state: FSMContext) -> None:
    if 'search' in message.text.casefold():
        await state.update_data(molecule_next_step=message.text)
        await state.set_state(SearchInfo.similarity_percent)
        await message.reply(
            'Enter similarity percent (between 40 and 100)', 
            reply_markup=types.ReplyKeyboardRemove()
        )
    else:
        await state.update_data(molecule_next_step=message.text)
        await state.set_state(SearchInfo.similarity_percent)
        await message.reply(
            'Enter SMILES and similarity percent (between 40 and 100), separated by a space', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.similarity_percent)
async def get_mols_by_similarity(message: types.Message, state: FSMContext) -> None:
    await state.update_data(similarity_percent=message.text)
    await state.set_state(SearchInfo.similarity_percent)
    data = await state.get_data()
    logger.debug(data)
    logger.debug(data.get('molecule_next_step'))
    if 'search' in data.get('molecule_next_step').casefold():
        mol_info = data.get('mol_info')
        mol_smiles = mol_info['molecule_structures']['canonical_smiles']
        similarity_percent = data.get('similarity_percent')
    else:
        smiles_similarity = data.get('similarity_percent')
        logger.debug(smiles_similarity)
        try:
            s = SmilesSimilarity(smiles_percent=smiles_similarity)
            smiles, similarity_percent = s.smiles_percent.split(' ', maxsplit=1)
            similarity_percent = int(similarity_percent)
            s = Smiles(smiles=smiles)
            mol_smiles = s.smiles
        except ValidationError as e:
            logger.error(e.errors())        
    
    if 40 <= similarity_percent <= 100:
        logger.debug(f"Looking for molecules similar to \'{mol_smiles}\' by {similarity_percent}%")
        similarity = new_client.similarity
        # Similarity score must be between 40 and 100
        res = similarity.filter(smiles=mol_smiles, similarity=similarity_percent).only(['molecule_chembl_id', 'similarity'])
        for i in res:
            logger.debug(i)
        reply_text = f"Results: \n{res}" if len(res) != 1 else f'No molecules with similarity of {similarity_percent}%'
        await message.reply(
                reply_text
            )
    else:
        await message.reply(
                f"Wrong similarity percent, should be between 40 and 100"
            )
