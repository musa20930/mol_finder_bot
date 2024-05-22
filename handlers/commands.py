from loguru import logger

from aiogram import Router
from aiogram.filters.command import Command, CommandStart
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.enums import ParseMode
from aiogram.types import (
    Message,
    ReplyKeyboardRemove
)


router = Router()


@router.message(CommandStart())
async def start(message: Message) -> None:
    await message.answer(
        text=f"Hi, @{message.from_user.username}! \nThis is MolFinder Bot!✨"
        f"Your go-to search for ChEMBL database. \n\n"
        f"Enter /search command to start. \n"
        f"Hope you find what you're looking for🙂\n\n"
        f"If you have any questions, try /help command or contact the author.\n"
        f"List of available commands can be found in the menu below👇",
        reply_markup=ReplyKeyboardRemove()
    )


@router.message(Command('help'))
async def help(message: Message) -> None:
    await message.answer(
        text=f"<b>Manual for MolFinder Bot</b>\n\n"
        f"<b>Where to start</b>\n"
        f"You can start by using /search command.\n"
        f"Choose any of the options below.\n"
        f"First row of buttons is used to get molecules by criteria, while second row is for retrieving by identifiers.\n"
        f"After finding molecule by one identifier, you can switch to any other by using /cancel command.\n\n"
        f"<b>Saving to file</b>\n"
        f"`Get Molecule Structure🧬` is for getting PNG or SVG image of molecule structure.\n"
        f"SDF, SMI AND CSV file formats are available as well.\n\n"
        f"<b>Canceling Tasks</b>\n"
        f"You can cancel any task, by typing /cancel command at any time.\n\n"
        f"<b>Contact</b>\n"
        f"If you still have some questions, "
        f"feel free to contanct the author at @musa_adham. 🙂\n"
        ,
        parse_mode=ParseMode.HTML, 
        reply_markup=ReplyKeyboardRemove()
    )


@router.message(Command('cancel'))
@router.message(F.text.casefold() == "cancel")
async def cancel_handler(message: Message, state: FSMContext) -> None:
    """Allow user to cancel any action."""
    current_state = await state.get_state()
    if current_state is None:
        return

    logger.info("Cancelling state %r", current_state)
    await state.clear()
    await message.answer(
        "Cancelled.",
        reply_markup=ReplyKeyboardRemove(),
    )