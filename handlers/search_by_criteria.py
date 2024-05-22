from loguru import logger

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.enums import ParseMode
from aiogram.types import (
    Message,
    CallbackQuery
)

from states import SearchInfo
from keyboards import Keyboard

from chembl_search_engine import (
    retrieve_by_mass,
    retrieve_by_logp
)


router = Router()


# MOLECULAR MASS SEARCH
@router.callback_query(F.data == 'mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await callback.message.edit_text(
        "Choose mass range type", 
        reply_markup=Keyboard().by_mass_search
    )
    await callback.answer()


@router.callback_query(F.data == 'mass_range')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_mass_comparison='range')
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text('Enter molecular mass range in g/mol with format "min mass, max mass"')
    await callback.answer()


@router.callback_query(F.data == 'gt_mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_mass_comparison='>')
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text("Enter molecular mass in g/mol ('>')")
    await callback.answer()


@router.callback_query(F.data == 'gte_mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_mass_comparison='>=')
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text("Enter molecular mass in g/mol ('>=')")
    await callback.answer()


@router.callback_query(F.data == 'lt_mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_mass_comparison='<')
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text("Enter molecular mass in g/mol ('<')")
    await callback.answer()


@router.callback_query(F.data == 'lte_mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_mass_comparison='<=')
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text("Enter molecular mass in g/mol ('<=')")
    await callback.answer()


@router.message(SearchInfo.mol_mass)
async def get_mol_by_mass(message: Message, state: FSMContext) -> None:
    data = await state.get_data()
    comparison_type = data.get('mol_mass_comparison')
    result = await retrieve_by_mass(message.text, comparison_type)

    if result is None:
        await message.reply(f"Wrong format of mass range")
    elif result:
        res_to_display = ''
        for i in result[:5]:
            res_to_display += f"'{i['molecule_chembl_id']}'\n"
        
        if len(result) > 1:
            if len(result) > 5:
                reply_text = f"{res_to_display[:600]}\n...and {len(result) - 5} more"
            else:
                reply_text = f"{res_to_display}"
        await message.reply(reply_text)
        await message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_similarity_res
        )
        await state.update_data(mol_series_info=result)
        await state.update_data(search_multiple=True)
    else:
        await state.set_state(SearchInfo.mol_mass)
        await message.answer('No such molecule found by mass')


# LOGP SEARCH
@router.callback_query(F.data == 'logp')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await callback.message.edit_text(
        "Choose mass range type", 
        reply_markup=Keyboard().by_logp_search
    )
    await callback.answer()


@router.callback_query(F.data == 'gt_logp')
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_logp_comparison='>')
    await state.set_state(SearchInfo.mol_logp)
    await callback.message.edit_text("Enter molecular LogP ('>')")
    await callback.answer()


@router.callback_query(F.data == 'gte_logp')
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_logp_comparison='>=')
    await state.set_state(SearchInfo.mol_logp)
    await callback.message.edit_text("Enter molecular LogP ('>=')")
    await callback.answer()


@router.callback_query(F.data == 'lt_logp')
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_logp_comparison='<')
    await state.set_state(SearchInfo.mol_logp)
    await callback.message.edit_text("Enter molecular LogP ('<')")
    await callback.answer()


@router.callback_query(F.data == 'lte_logp')
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    await state.update_data(mol_logp_comparison='<=')
    await state.set_state(SearchInfo.mol_logp)
    await callback.message.edit_text("Enter molecular LogP ('<=')")
    await callback.answer()


@router.message(SearchInfo.mol_logp)
async def get_mol_by_logp(message: Message, state: FSMContext) -> None:
    data = await state.get_data()
    comparison_type = data.get('mol_logp_comparison')
    result = await retrieve_by_logp(message.text, comparison_type)

    if result is None:
        await message.reply(f"Wrong format of LogP range")
    elif result:
        res_to_display = ''
        for i in result[:5]:
            res_to_display += f"'{i['molecule_chembl_id']}'\n"
        
        if len(result) > 1:
            if len(result) > 5:
                reply_text = f"{res_to_display[:600]}\n...and {len(result) - 5} more"
            else:
                reply_text = f"{res_to_display}"
        await message.reply(reply_text)
        await message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_similarity_res
        )
        await state.update_data(mol_series_info=result)
        await state.update_data(search_multiple=True)
    else:
        await state.set_state(SearchInfo.mol_logp)
        await message.answer('No such molecule found by LogP')

