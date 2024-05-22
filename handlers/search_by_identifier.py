from loguru import logger
import re

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram.filters.command import Command
from aiogram import F
from aiogram.types import (
    Message,
    CallbackQuery
)

from states import SearchInfo
from keyboards import Keyboard
from chembl_search_engine import (
    retrieve_by_name, 
    retrieve_by_id,
    retrieve_by_inchi,
    retrieve_by_smiles
)
from message_processing import display_mol_info


router = Router()


@router.message(Command('search'))
async def search(message: Message, state: FSMContext) -> None:
    await state.clear()
    await message.answer(
        "Choose search type", 
        reply_markup=Keyboard().search
    )


@router.callback_query(F.data == 'name')
async def search_by_name(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.molecule_name)
    await callback.message.edit_text('Enter Name of molecule')
    await callback.answer()


@router.callback_query(F.data == 'chembl_id')
async def search_by_id(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.chembl_id)
    await callback.message.edit_text('Enter ChEMBL ID')
    await callback.answer()


@router.callback_query(F.data == 'inchi')
async def search_by_id(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.inchi_key)
    await callback.message.edit_text('Enter INCHI key')
    await callback.answer()


@router.callback_query(F.data == 'smiles')
async def search_by_id(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.smiles)
    await callback.message.edit_text('Enter SMILES')
    await callback.answer()


@router.message(SearchInfo.molecule_name)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(molecule_name=message.text)
    result = await retrieve_by_name(message.text)
    if result:
        await display_mol_info(message, state, result)
    else:
        await state.set_state(SearchInfo.molecule_name)
        await message.answer('No such molecule found by name')


@router.message(SearchInfo.chembl_id)
async def get_mol_by_id(message: Message, state: FSMContext) -> None:
    await state.update_data(chembl_id=message.text)
    result = await retrieve_by_id(message.text)
    if result:
        await display_mol_info(message, state, result)
    else:
        await state.set_state(SearchInfo.chembl_id)
        await message.answer('No such molecule found by ChEMBL ID')


@router.message(SearchInfo.inchi_key)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(inchi_key=message.text)
    result = await retrieve_by_inchi(message.text)
    if result:
        await display_mol_info(message, state, result)
    else:
        await state.set_state(SearchInfo.inchi_key)
        await message.answer('No such molecule found by INCHI key')


@router.message(SearchInfo.smiles)
async def get_mol_by_smiles(message: Message, state: FSMContext) -> None:
    await state.update_data(smiles=message.text)
    result = await retrieve_by_smiles(message.text)
    if result:
        await display_mol_info(message, state, result)
    else:
        await state.set_state(SearchInfo.smiles)
        await message.answer('No such molecule found by SMILES')

