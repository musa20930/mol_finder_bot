from aiogram import Bot, Dispatcher, types
import asyncio
import logging
import re
import os
import sys
from loguru import logger
from aiogram.filters.command import Command, CommandStart
from aiogram.types import BotCommand
from aiogram import F
from chembl_webresource_client.new_client import new_client


API_TOKEN = os.getenv('MOL_FINDER_BOT_TOKEN')
bot = Bot(token=API_TOKEN)
dp = Dispatcher() # отслеживает все входящие события


class InterceptHandler(logging.Handler):
    """Intercept log messages to use loguru instead of logging."""
    def emit(self, record):
        level = logger.level(record.levelname).name
        logger.log(level, record.getMessage())


@dp.message(CommandStart())
async def start(message: types.Message):
    keyboard = [
        [types.KeyboardButton(text='By name'), 
         types.KeyboardButton(text='By Chembl ID')],
    ]
    reply_markup = types.ReplyKeyboardMarkup(keyboard=keyboard, 
                                             resize_keyboard=True, 
                                             input_field_placeholder="Choose option")
    await message.answer(text="Choose search type", 
                         reply_markup=reply_markup)


@dp.message(F.text.casefold() == 'by name')
async def text(message: types.Message):
    await message.reply('Enter name of molecule', reply_markup=types.ReplyKeyboardRemove())


@dp.message(F.text.casefold() == 'by chembl id')
async def text(message: types.Message):
    await message.reply('Enter id', reply_markup=types.ReplyKeyboardRemove())


@dp.message(F.text)
async def text(message: types.Message):
    if re.search(r"^(\d+)$", message.text):
        await text_by_id(message)
    else:
        await text_by_name(message)


async def text_by_name(message: types.Message):
    name = message.text
    logger.debug(f'Searching \'{name}\' by name')
    molecule = new_client.molecule
    result = molecule.filter(pref_name__iexact=name).only(['molecule_chembl_id', 
                                                           'pref_name', 
                                                           'molecule_structures'])
    logger.debug(f'Search by name {result = }')
    if result:
        await message.answer(f"Name: {result[0]['pref_name']} \n"
                            f"ID: {result[0]['molecule_chembl_id']} \n"
                            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
                            reply_markup=types.ReplyKeyboardRemove())
    else:
        await message.answer('No such molecule', reply_markup=types.ReplyKeyboardRemove())


async def text_by_id(message: types.Message):
    chembl_id = message.text
    logger.debug(f'Searching by id \'{chembl_id}\'')
    molecule = new_client.molecule
    result = molecule.filter(molecule_chembl_id=f'CHEMBL{chembl_id}').only(['molecule_chembl_id', 
                                                                       'pref_name', 
                                                                       'molecule_structures']) 
    logger.debug(f'Search by id {result = }')
    if result:
        await message.reply(f"Name: {result[0]['pref_name']} \n"
                            f"ID: {result[0]['molecule_chembl_id']} \n"
                            f"SMILES: {result[0]['molecule_structures']['canonical_smiles']}", 
                            reply_markup=types.ReplyKeyboardRemove())
    else:
        await message.reply('No such molecule', reply_markup=types.ReplyKeyboardRemove())


@dp.message(Command('stop'))
async def stop(message: types.Message):
    await message.answer(f"Stop {message.from_user.id}")


def add_logging():
    if "--debug" not in sys.argv:
        """Used to add debug mode option 
        for logging while testing with Dash's debug mode.
        This part is also responsible for removing debug level 
        from terminal output. Insted debug level is written 
        into the log file 'chembot_app.log'"""
        logger.remove()
        logger.add(sys.stdout, level="INFO")
    log_level = "DEBUG"
    log_format = "{time:YYYY-DD-MM > HH:mm:ss} | {level} | {message}"
    logger.add("mol_finder_bot.log", level=log_level, 
               format=log_format, colorize=False, 
               backtrace=True, diagnose=True)
    logging.getLogger('aiogram').setLevel(logging.DEBUG)
    logging.getLogger('aiogram').addHandler(InterceptHandler())
    logging.getLogger('asyncio').setLevel(logging.DEBUG)
    logging.getLogger('asyncio').addHandler(InterceptHandler())


async def main():
    add_logging()
    await bot.set_my_commands([
        BotCommand(command="start", description="Запустить бота"),
        BotCommand(command="help", description="Помощь"),
        BotCommand(command="test", description="Тест"),
        BotCommand(command="form", description="Форма"),
        BotCommand(command="menu", description="Меню"),
    ])
    await dp.start_polling(bot)


if __name__ == '__main__':    
    asyncio.run(main())

