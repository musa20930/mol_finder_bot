from loguru import logger

from aiogram import types, Router
from aiogram.fsm.context import FSMContext
from aiogram.filters.command import Command
from aiogram import F

from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw

from pydantic import ValidationError
from validation_models import Smiles
from states import SearchInfo


search_router = Router()

# BUTTONS
search_by_id_btn = types.KeyboardButton(text='Search by ChEMBL ID')
search_by_name_btn = types.KeyboardButton(text='Search by Name')
search_by_inchi_btn = types.KeyboardButton(text='Search by INCHI key')
search_by_smiles_btn = types.KeyboardButton(text='Search by SMILES')
get_mol_structure_btn = types.KeyboardButton(text='Get Molecule Structure🧬')
search_by_similarity_btn = types.KeyboardButton(text='Search by Similarity❄️❄️')

keyboard_after_search = [
    [search_by_similarity_btn, 
    get_mol_structure_btn],
    [search_by_name_btn, 
     search_by_id_btn,
    search_by_inchi_btn,
    search_by_smiles_btn],
]
reply_markup_after_search = types.ReplyKeyboardMarkup(keyboard=keyboard_after_search, 
                                                      resize_keyboard=True, 
                                                      input_field_placeholder="Choose option")


@search_router.message(Command('search'))
async def search(message: types.Message, state: FSMContext) -> None:
    await state.set_state(SearchInfo.search_type)
    keyboard = [
        [types.KeyboardButton(text='By Name'), 
         types.KeyboardButton(text='By Chembl ID'),
         types.KeyboardButton(text='By INCHI key'),
         types.KeyboardButton(text='By SMILES')],
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
    if 'search' not in message.text:
        await state.update_data(search_type=message.text)
    else:
        await state.update_data(molecule_next_step=message.text)
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


@search_router.message(SearchInfo.search_type, F.text.casefold() == 'by inchi key')
@search_router.message(SearchInfo.molecule_next_step, F.text.casefold() == 'search by inchi key')
async def search_by_inchi(message: types.Message, state: FSMContext) -> None:
    if 'search' not in message.text:
        await state.update_data(search_type=message.text)
    else:
        await state.update_data(molecule_next_step=message.text)
    await state.set_state(SearchInfo.inchi_key)
    await message.reply(
        "Enter INCHI key", 
        reply_markup=types.ReplyKeyboardRemove()
    )


@search_router.message(SearchInfo.search_type, F.text.casefold() == 'by smiles')
@search_router.message(SearchInfo.molecule_next_step, F.text.casefold() == 'search by smiles')
async def search_by_inchi(message: types.Message, state: FSMContext) -> None:
    if 'search' not in message.text:
        await state.update_data(search_type=message.text)
    else:
        await state.update_data(molecule_next_step=message.text)
    await state.set_state(SearchInfo.smiles)
    await message.reply(
        "Enter SMILES", 
        reply_markup=types.ReplyKeyboardRemove()
    )


@search_router.message(SearchInfo.molecule_name)
async def get_mol_by_name(message: types.Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    name = message.text
    logger.debug(f'Searching \'{name}\' by name')
    molecule = new_client.molecule
    result = molecule.filter(
        molecule_synonyms__molecule_synonym__iexact=name).only(['molecule_chembl_id', 
                                                                'pref_name', 
                                                                'molecule_structures'])
    logger.debug(f'Search by name {result = }')
    if result:
        await state.set_state(SearchInfo.molecule_next_step)
        pref_name = result[0]['pref_name']
        await message.answer(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        
        await message.answer(
            text="You can now search for new molecule by name "
            "or rerun /search command to search by ChEMBL ID", 
            reply_markup=reply_markup_after_search
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
    chembl_id = message.text
    logger.debug(f'Searching by id \'{chembl_id}\'')
    molecule = new_client.molecule
    result = molecule.filter(molecule_chembl_id=f'CHEMBL{chembl_id}').only(['molecule_chembl_id', 
                                                                            'pref_name', 
                                                                            'molecule_structures']) 
    logger.debug(f'Search by id {result = }')
    if result:
        await state.set_state(SearchInfo.molecule_next_step)
        pref_name = result[0]['pref_name']
        await message.reply(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        
        await message.answer(
            text="You can now search for new molecule by ChEMBL ID "
            "or rerun /search command to search by name", 
            reply_markup=reply_markup_after_search
        )
    else:
        await state.set_state(SearchInfo.chembl_id)
        await message.reply(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.inchi_key)
async def get_mol_by_name(message: types.Message, state: FSMContext) -> None:
    await state.update_data(inchi_key=message.text)
    inchi = message.text
    logger.debug(f'Searching \'{inchi}\' by INCHI key')
    molecule = new_client.molecule
    result = molecule.filter(
        molecule_structures__standard_inchi_key=inchi).only(['molecule_chembl_id', 
                                                            'pref_name', 
                                                            'molecule_structures'])
    logger.debug(f'Search by INCHI key {result = }')
    if result:
        await state.set_state(SearchInfo.molecule_next_step)
        pref_name = result[0]['pref_name']
        await message.answer(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        
        await message.answer(
            text="You can now search for new molecule by INCHI key "
            "or rerun /search command to search by ChEMBL ID", 
            reply_markup=reply_markup_after_search
        )
    else:
        await state.set_state(SearchInfo.inchi_key)
        await message.answer(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.smiles)
async def get_mol_by_smiles(message: types.Message, state: FSMContext) -> None:
    await state.update_data(smiles=message.text)
    smiles = message.text

    logger.debug(f'Validating SMILES: \'{smiles}\'')
    try:
        s = Smiles(smiles=smiles)
        mol_smiles = s.smiles
    except ValidationError as e:
        logger.error(e.errors())
        await message.reply(
            f"Wrong format, should be: <smiles> <percent>"
        )
    except Exception as e:
        logger.error(e)
        await message.reply(
            f"Something went wrong"
        )

    logger.debug(f'Searching \'{smiles}\' by SMILES')
    molecule = new_client.molecule
    result = molecule.filter(smiles=mol_smiles).only(['molecule_chembl_id', 
                                                      'pref_name', 
                                                      'molecule_structures'])
    logger.debug(f'Search by SMILES {result = }')
    if result:
        await state.set_state(SearchInfo.molecule_next_step)
        pref_name = result[0]['pref_name']
        await message.answer(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL database'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )

        await state.update_data(mol_info=result[0])
        
        await message.answer(
            text="You can now search for new molecule by SMILES "
            "or rerun /search command to search by ChEMBL ID", 
            reply_markup=reply_markup_after_search
        )
    else:
        await state.set_state(SearchInfo.smiles)
        await message.answer(
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


@search_router.message(SearchInfo.molecule_next_step, F.text.casefold().startswith('search by similarity'))
async def search_by_similarity(message: types.Message, state: FSMContext) -> None:
    if 'search' in message.text.casefold():
        await state.update_data(molecule_next_step=message.text)
        await state.set_state(SearchInfo.similarity_percent)
        await message.reply(
            'Enter similarity percent (between 40 and 100)', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@search_router.message(SearchInfo.similarity_percent)
async def get_mols_by_similarity(message: types.Message, state: FSMContext) -> None:
    await state.update_data(similarity_percent=message.text)
    await state.set_state(SearchInfo.similarity_percent)
    data = await state.get_data()

    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    similarity_percent = data.get('similarity_percent')
    similarity_percent = float(similarity_percent)
    
    if 40 <= similarity_percent <= 100:
        logger.debug(f"Looking for molecules similar to \'{mol_smiles}\' by {similarity_percent}%")
        similarity = new_client.similarity
        # Similarity score must be between 40 and 100
        res = similarity.filter(smiles=mol_smiles, similarity=similarity_percent).only(['molecule_chembl_id', 'similarity'])
        res_to_display = ''
        for i in res:
            logger.debug(i)
            res_to_display += f"'{i['molecule_chembl_id']}', {i['similarity']}%\n"
        
        if len(res) > 1:
            reply_text = f"Molecules with similarity of {similarity_percent}%: \n{res_to_display}"
        else:
            reply_text = f'No molecules found with similarity of {similarity_percent}%'
        await message.reply(
                reply_text
            )
    else:
        await message.reply(
                f"Wrong similarity percent, should be between 40 and 100"
            )
