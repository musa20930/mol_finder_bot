from loguru import logger
import re

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.types import (
    Message,
    CallbackQuery
)

from states import SearchInfo
from keyboards import Keyboard
from chembl_search_engine import (
    retrieve_similar_molecules,
    retrieve_by_name, 
    retrieve_by_id,
    retrieve_by_inchi,
    retrieve_by_smiles
)
from message_processing import display_mol_info_for_sim


router = Router()


@router.callback_query(F.data == 'similarity')
async def search_by_similarity(callback: CallbackQuery, state: FSMContext) -> None:
    data = await state.get_data()
    if 'mol_info' in data:
        await state.set_state(SearchInfo.similarity_percent)
        await callback.message.answer('Enter similarity percent (between 40 and 100)')
    else:
        await callback.message.edit_text(
            "Choose similarity search type", 
            reply_markup=Keyboard().similarity_search
        )


@router.callback_query(F.data == 'name_sim')
async def search_by_name_sim(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.molecule_name_sim)
    await callback.message.edit_text('Enter Name of molecule')
    await callback.answer()


@router.callback_query(F.data == 'chembl_id_sim')
async def search_by_name_sim(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.chembl_id_sim)
    await callback.message.edit_text('Enter ChEMBL ID')
    await callback.answer()


@router.callback_query(F.data == 'inchi_sim')
async def search_by_name_sim(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.inchi_key_sim)
    await callback.message.edit_text('Enter INCHI key')
    await callback.answer()


@router.callback_query(F.data == 'smiles_sim')
async def search_by_name_sim(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.smiles_sim)
    await callback.message.edit_text('Enter SMILES')
    await callback.answer()


@router.message(SearchInfo.molecule_name_sim)
async def get_mol_for_similarity(message: Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    result = await retrieve_by_name(message.text)
    await state.update_data(mol_info=result[0])
    if result:
        await display_mol_info_for_sim(message, state, result)
        await state.set_state(SearchInfo.similarity_percent)
        await message.answer('Enter similarity percent (between 40 and 100)')
    else:
        await state.set_state(SearchInfo.molecule_name_sim)
        await message.answer('No such molecule found by name')
    

@router.message(SearchInfo.chembl_id_sim)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    result = await retrieve_by_id(message.text)
    await state.update_data(mol_info=result[0])
    if result:
        await display_mol_info_for_sim(message, state, result)
        await state.set_state(SearchInfo.similarity_percent)
        await message.answer('Enter similarity percent (between 40 and 100)')
    else:
        await state.set_state(SearchInfo.chembl_id_sim)
        await message.answer('No such molecule found by ChEMBL ID')


@router.message(SearchInfo.inchi_key_sim)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    result = await retrieve_by_inchi(message.text)
    await state.update_data(mol_info=result[0])
    if result:
        await display_mol_info_for_sim(message, state, result)
        await state.set_state(SearchInfo.similarity_percent)
        await message.answer('Enter similarity percent (between 40 and 100)')
    else:
        await state.set_state(SearchInfo.inchi_key_sim)
        await message.answer('No such molecule found by INCHI key')


@router.message(SearchInfo.smiles_sim)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    result = await retrieve_by_smiles(message.text)
    await state.update_data(mol_info=result[0])
    if result:
        await display_mol_info_for_sim(message, state, result)
        await state.set_state(SearchInfo.similarity_percent)
        await message.answer('Enter similarity percent (between 40 and 100)')
    else:
        await state.set_state(SearchInfo.smiles_sim)
        await message.answer('No such molecule found by SMILES')


@router.message(SearchInfo.similarity_percent)
async def get_mols_by_similarity(message: Message, state: FSMContext) -> None:
    await state.update_data(similarity_percent=message.text)
    await state.set_state(SearchInfo.similarity_percent)
    data = await state.get_data()

    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    similarity_percent = data.get('similarity_percent')
    if re.search(f"^[+]?([0-9]+([.][0-9]*)?|[.][0-9]+)$", similarity_percent):
        similarity_percent = float(similarity_percent)
    else:
        await message.reply(f"Wrong format, should be a positive number")
    
    if 40 <= similarity_percent <= 100:
        res = await retrieve_similar_molecules(mol_smiles, similarity_percent)
        res_to_display = ''
        for i in res[:5]:
            res_to_display += f"'{i['molecule_chembl_id']}', {i['similarity']}%\n"
        
        if len(res) > 1:
            if len(res) > 5:
                reply_text = f"Molecules with similarity of {similarity_percent}%: \n{res_to_display}\n...and {len(res) - 5} more"
            else:
                reply_text = f"Molecules with similarity of {similarity_percent}%: \n{res_to_display}"
        else:
            reply_text = f'No molecules found with similarity of {similarity_percent}%'
        await message.reply(reply_text)
        await message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_criteria_search_res
        )
        await state.update_data(similarity_percent=None)
        await state.update_data(mol_series_info=res)
        await state.update_data(search_multiple=True)
    else:
        await message.reply(f"Wrong similarity percent, should be between 40 and 100")

