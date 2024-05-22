from loguru import logger
import re

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram.filters.command import Command
from aiogram import F
from aiogram.enums import ParseMode
from aiogram.types import (
    Message,
    CallbackQuery
)

from states import SearchInfo
from keyboards import Keyboard

from chembl_search_engine import (
    retrieve_by_mass
)


router = Router()


@router.callback_query(F.data == 'mass')
async def search_by_name(callback: CallbackQuery, state: FSMContext) -> None:
    await state.set_state(SearchInfo.mol_mass)
    await callback.message.edit_text('Enter molecular mass in g/mol ')
    await callback.answer()


@router.message(SearchInfo.mol_mass)
async def get_mol_by_name(message: Message, state: FSMContext) -> None:
    await state.update_data(mol_mass=message.text)
    if re.search(f"^[+]?([0-9]+([.][0-9]*)?|[.][0-9]+)$", message.text):
        similarity_percent = float(similarity_percent)
    else:
        await message.reply(f"Wrong format, should be a positive number")
    result = await retrieve_by_mass(message.text)
    if result:
        res_to_display = ''
        for i in result[:5]:
            res_to_display += f"'{i['molecule_chembl_id']}', {i['similarity']} g/mol\n"
        
        if len(result) > 1:
            if len(result) > 5:
                reply_text = f"Molecules with mass of {similarity_percent} g/mol: \n{res_to_display[:600]}\n...and {len(result) - 5} more"
            else:
                reply_text = f"Molecules with mass of {similarity_percent} g/mol: \n{res_to_display}"
        else:
            reply_text = f'No molecules found with mass of {similarity_percent} g/mol'
        await message.reply(reply_text)
    else:
        await state.set_state(SearchInfo.mol_mass)
        await message.answer('No such molecule')
