from loguru import logger

from aiogram import Router
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.types import (
    Message,
    ReplyKeyboardRemove
)

from keyboards import Keyboard
from chembl_search_engine import (
    retrieve_by_name,
    retrieve_by_id,
    retrieve_by_inchi,
    retrieve_by_smiles
)
from message_processing import display_mol_info


router = Router()
retrieve_dict = {
    'name': retrieve_by_name,
    'id': retrieve_by_id,
    'inchi': retrieve_by_inchi,
    'smiles': retrieve_by_smiles
}


@router.message(F.text)
async def search_by_any(message: Message, state: FSMContext) -> None:
    for func in retrieve_dict.values():
        result = await func(message.text)
        if result:
            await display_mol_info(message, state, result)
            break
    else:
        await message.answer(
            'No such molecule', 
            reply_markup=ReplyKeyboardRemove()
        )
    

