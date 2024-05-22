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
    retrieve_by_logp,
    retrieve_by_lipinski_rule, 
    retrieve_drugs
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
@router.callback_query(F.data == 'gt_mass')
@router.callback_query(F.data == 'gte_mass')
@router.callback_query(F.data == 'lt_mass')
@router.callback_query(F.data == 'lte_mass')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    comparison_dict = {'mass_range': 'range', 'gt_mass': '>', 'gte_mass': '>=', 'lt_mass': '<', 'lte_mass': '<='}
    await state.update_data(mol_mass_comparison=comparison_dict[callback.data])
    await state.set_state(SearchInfo.mol_mass)
    if comparison_dict[callback.data] == 'range':
        reply_text = 'Enter molecular mass range in g/mol with format "min mass, max mass"'
    else:
        reply_text = f"Enter molecular mass in g/mol ('{comparison_dict[callback.data]}')"
    await callback.message.edit_text(reply_text)
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
            res_to_display += f"{i['molecule_chembl_id']}\n"
        
        if len(result) > 5:
            reply_text = f"{res_to_display}\n...and {len(result) - 5} more"
        else:
            reply_text = f"{res_to_display}"
        await message.reply(reply_text)
        await message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_criteria_search_res
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
@router.callback_query(F.data == 'gte_logp')
@router.callback_query(F.data == 'lt_logp')
@router.callback_query(F.data == 'lte_logp')
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    comparison_dict = {'gt_logp': '>', 'gte_logp': '>=', 'lt_logp': '<', 'lte_logp': '<='}
    await state.update_data(mol_logp_comparison=comparison_dict[callback.data])
    await state.set_state(SearchInfo.mol_logp)
    await callback.message.edit_text(f"Enter molecular LogP ('{comparison_dict[callback.data]}')")
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
            res_to_display += f"{i['molecule_chembl_id']}\n"
        
        if len(result) > 5:
            reply_text = f"{res_to_display}\n...and {len(result) - 5} more"
        else:
            reply_text = f"{res_to_display}"
        await message.reply(reply_text)
        await message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_criteria_search_res
        )
        await state.update_data(mol_series_info=result)
        await state.update_data(search_multiple=True)
    else:
        await state.set_state(SearchInfo.mol_logp)
        await message.answer('No such molecule found by LogP')


# LIPINSKI SEARCH
@router.callback_query(F.data == 'lipinski')
async def search_by_mass(callback: CallbackQuery, state: FSMContext) -> None:
    await callback.message.edit_text(
        "Choose acceptable number of 'Rule-of-Five' violations", 
        reply_markup=Keyboard().lipinski_search
    )
    await callback.answer()


@router.callback_query(F.data.endswith('_violations'))
async def search_by_logp(callback: CallbackQuery, state: FSMContext) -> None:
    result = await retrieve_by_lipinski_rule(int(callback.data[0]))

    if result:
        res_to_display = ''
        for i in result[:5]:
            res_to_display += f"{i['molecule_chembl_id']}\n"
        
        if len(result) > 5:
            reply_text = f"{res_to_display}\n...and {len(result) - 5} more"
        else:
            reply_text = f"{res_to_display}"
        await callback.message.reply(reply_text)
        await callback.message.answer(
            'Want to save to file?', 
            reply_markup=Keyboard().save_criteria_search_res
        )
        await state.update_data(mol_series_info=result)
        await state.update_data(search_multiple=True)
    else:
        await callback.message.answer('No such molecule found by Lipinski rule')
    await callback.answer()


# DRUG SEARCH
@router.callback_query(F.data == 'drugs')
async def search_drugs_by_year(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.drugs_year)
    await callback.message.edit_text("Enter year")
    await callback.answer()


@router.message(SearchInfo.drugs_year)
async def get_(message: Message, state: FSMContext) -> None:
    await state.set_state(SearchInfo.drugs_amount)
    if not message.text.isdigit():
        await message.reply("Wrong format of year")
    else:
        await state.update_data(drugs_year=message.text)
        await message.answer("Enter amount of approved drugs you would like to get")


@router.message(SearchInfo.drugs_amount)
async def search_drugs(message: Message, state: FSMContext) -> None:
    await state.set_state(SearchInfo.drugs_amount)
    if not message.text.isdigit():
        await message.reply("Wrong format of amount of drugs")
    else:
        data = await state.get_data()
        logger.debug(f"data: {data}, message.text: {message.text}")
        year, amount = int(data.get('drugs_year')), int(message.text)

        result = await retrieve_drugs(year, amount)
        if result:
            res_to_display = ''
            for i in result[:5]:
                res_to_display += f"{i['molecule_chembl_id']}\n"
            
            if len(result) > 5:
                reply_text = f"{res_to_display}\n...and {len(result) - 5} more"
            else:
                reply_text = f"{res_to_display}"
            await message.reply(reply_text)
            await message.answer(
                'Want to save to file?', 
                reply_markup=Keyboard().save_criteria_search_res
            )
            await state.update_data(mol_series_info=result)
            await state.update_data(search_multiple=True)
        else:
            await message.answer('No such drugs found')

