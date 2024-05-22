from loguru import logger

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.enums import ParseMode
from aiogram.types import (
    CallbackQuery
)

from keyboards import Keyboard
from chembl_search_engine import (
    retrieve_properties,
    save_properties,
    standardize_mol
)
from message_processing import send_structure_png


router = Router()


@router.callback_query(F.data == 'structure')
async def get_mol_structure(callback: CallbackQuery, state: FSMContext) -> None:
    data = await state.get_data()
    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    mol_img = await send_structure_png(mol_smiles, mol_info['molecule_chembl_id'])
    # await callback.message.answer_document(mol_img)
    await callback.message.answer_photo(photo=mol_img)
    await callback.answer()


@router.callback_query(F.data == 'standardize')
async def get_standard_mol(callback: CallbackQuery, state: FSMContext) -> None:
    data = await state.get_data()
    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    sdf_file = await standardize_mol(mol_smiles, mol_info['molecule_chembl_id'])
    await callback.message.answer_document(sdf_file)
    await callback.answer()


@router.callback_query(F.data == 'properties')
async def get_properties(callback: CallbackQuery, state: FSMContext) -> None:
    data = await state.get_data()
    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    descs = await retrieve_properties(mol_smiles)
    if len(descs) == 1:
        reply_text = f"<b>{mol_info['molecule_chembl_id']} Properties</b>:\n\n"
        for desc, value in descs[0].items():
            reply_text += f"<b>{desc}</b>: {value}\n"
        await callback.message.answer(reply_text, parse_mode=ParseMode.HTML)
        await callback.message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().properties
        )
        await state.update_data(mol_descriptors=descs)
    else:
        csv_file = await save_properties(descs, mol_smiles, mol_info['molecule_chembl_id'])
        await callback.message.answer_document(csv_file)
    await callback.answer()


@router.callback_query(F.data == 'properties_csv')
async def get_properties_in_csv(callback: CallbackQuery, state: FSMContext) -> None:
    data = await state.get_data()
    mol_info = data.get('mol_info')
    mol_smiles = mol_info['molecule_structures']['canonical_smiles']
    descs = data.get('mol_descriptors')
    csv_file = await save_properties(descs, mol_smiles, mol_info['molecule_chembl_id'])
    await callback.message.answer_document(csv_file)
    await callback.message.delete()
    await callback.answer()

