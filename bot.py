"""
- To run with debug mode use command:
poetry run python bot.py --debug

- If you don't want to see debug log level:
poetry run python bot.py
"""

from aiogram import Bot, Dispatcher, types, Router
import asyncio
import logging
import re
import os
import sys
from loguru import logger
from aiogram.filters.command import Command, CommandStart
from aiogram.fsm.context import FSMContext
from aiogram.fsm.state import State, StatesGroup
from aiogram import F
from chembl_webresource_client.new_client import new_client


form_router = Router()


class Form(StatesGroup):
    search_type = State()
    molecule_name = State()
    chembl_id = State()


class InterceptHandler(logging.Handler):
    """Class for intercepting log messages to use loguru instead of logging."""
    def emit(self, record) -> None:
        level = logger.level(record.levelname).name
        logger.log(level, record.getMessage())


@form_router.message(CommandStart())
async def start(message: types.Message, state: FSMContext) -> None:
    await state.set_state(Form.search_type)
    keyboard = [
        [types.KeyboardButton(text='By name'), 
         types.KeyboardButton(text='By Chembl ID')],
    ]
    reply_markup = types.ReplyKeyboardMarkup(keyboard=keyboard, 
                                             resize_keyboard=True, 
                                             input_field_placeholder="Choose option")
    await message.answer(
        text="Choose search type", 
        reply_markup=reply_markup
    )


@form_router.message(Form.search_type, F.text.casefold() == 'by name')
async def text(message: types.Message, state: FSMContext) -> None:
    await state.update_data(search_type=message.text)
    await state.set_state(Form.molecule_name)
    await message.reply(
        'Enter name of molecule', 
        reply_markup=types.ReplyKeyboardRemove()
    )


@form_router.message(Form.search_type, F.text.casefold() == 'by chembl id')
async def text(message: types.Message, state: FSMContext) -> None:
    await state.update_data(search_type=message.text)
    await state.set_state(Form.chembl_id)
    await message.reply(
        'Enter id', 
        reply_markup=types.ReplyKeyboardRemove()
    )


@form_router.message(Form.molecule_name)
async def text_by_name(message: types.Message, state: FSMContext) -> None:
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
        pref_name = result[0]['pref_name']
        await message.answer(
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL DB'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            eply_markup=types.ReplyKeyboardRemove()
        )
    else:
        await message.answer(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@form_router.message(Form.chembl_id)
async def text_by_id(message: types.Message, state: FSMContext) -> None:
    await state.update_data(chembl_id=message.text)
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
            f"Name: {pref_name if pref_name is not None else 'not provided in ChEMBL DB'} \n"
            f"ID: {result[0]['molecule_chembl_id']} \n"
            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
            reply_markup=types.ReplyKeyboardRemove()
        )
    else:
        await message.reply(
            'No such molecule', 
            reply_markup=types.ReplyKeyboardRemove()
        )


@form_router.message(Command('cancel'))
@form_router.message(F.text.casefold() == "cancel")
async def cancel_handler(message: types.Message, state: FSMContext) -> None:
    """Allow user to cancel any action."""
    current_state = await state.get_state()
    if current_state is None:
        return

    logger.info("Cancelling state %r", current_state)
    await state.clear()
    await message.answer(
        "Cancelled.",
        reply_markup=types.ReplyKeyboardRemove(),
    )


def add_logging() -> None:
    """Add debug mode option and replace logging with loguru."""
    if "--debug" not in sys.argv:
        """Used to add debug mode option.
        This part is also responsible for removing debug level 
        from terminal output. Insted debug level is written 
        into the log file 'mol_finder_bot.log'"""
        logger.remove()
        logger.add(sys.stdout, level="INFO")
    log_level = "DEBUG"
    log_format = "{time:YYYY-DD-MM > HH:mm:ss} | {level} | {message}"
    logger.add("mol_finder_bot.log", level=log_level, 
               format=log_format, colorize=False, 
               backtrace=True, diagnose=True)
    # replace logging with loguru
    logging.getLogger('aiogram').setLevel(logging.DEBUG)
    logging.getLogger('aiogram').addHandler(InterceptHandler())
    logging.getLogger('asyncio').setLevel(logging.DEBUG)
    logging.getLogger('asyncio').addHandler(InterceptHandler())


async def main() -> None:
    add_logging()

    API_TOKEN = os.getenv('MOL_FINDER_BOT_TOKEN')
    bot = Bot(token=API_TOKEN)
    dp = Dispatcher() # отслеживает все входящие события
    # Add command menu with set_my_commands
    await bot.set_my_commands([
        types.BotCommand(command="start", description="Запустить бота"),
        types.BotCommand(command="help", description="Помощь"),
        types.BotCommand(command="test", description="Тест"),
        types.BotCommand(command="form", description="Форма"),
        types.BotCommand(command="menu", description="Меню"),
    ])
    dp.include_router(form_router)
    await dp.start_polling(bot)


if __name__ == '__main__':    
    asyncio.run(main())

