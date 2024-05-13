from loguru import logger

from aiogram import types, Router
from aiogram.filters.command import Command, CommandStart
from aiogram.fsm.context import FSMContext
from aiogram import F
from aiogram.enums import ParseMode


command_router = Router()


@command_router.message(CommandStart())
async def start(message: types.Message) -> None:
    await message.answer(
        text=f"Hi, @{message.from_user.username}! \nThis is MolFinder Bot!✨"
        f"Your go-to search for ChEMBL database. \n\n"
        f"Enter /search command to start. \n"
        f"Hope you find what you're looking for🙂\n\n"
        f"If you have any questions, try /help command or contact the author.\n"
        f"List of available commands can be found on the left in the menu below👇",
        reply_markup=types.ReplyKeyboardRemove()
    )


@command_router.message(Command('help'))
async def help(message: types.Message) -> None:
    await message.answer(
        text=f"<b>Manual for MolFinder Bot</b>\n\n"
        f"<b>Where to start</b>\n"
        f"You can start by using /search command.\n"
        f"After which, you will see 3 buttons:\n"
        f"\t\tSearch by Name\n"
        f"\t\tSearch by ChEMBL ID\n"
        f"\t\tSearch by Similarity\n\n"
        f"It is recommended to use first 2 options, as those are more stable. "
        f"Recommended search by similarity is mentioned below.\n\n"
        f"<b>Look up molecules and their info</b>\n"
        f"If you press one of the 3 buttons and type in your search query \n"
        f"buttons will change to the following:\n"
        f"\t\tSearch by Similarity❄️❄️\n"
        f"\t\tGet Molecule Structure🧬\n"
        f"\t\tSearch by Name\n\n"
        f"`Search by Similarity❄️❄️` "
        f"will bring you to the recommended way to search by similarity.\n"
        f"`Get Molecule Structure🧬` is for getting .png image of molecule structure.\n"
        f"`Search by Name` button is identical to the button mentioned above.\n\n"
        f"<b>Canceling Tasks</b>\n"
        f"You can cancel any task, by typing /cancel command at any time.\n\n"
        f"<b>Contact Information</b>\n"
        f"If you still have some questions, "
        f"feel free to contanct the author at @musa_adham.\n"
        f"Thank you for using MolFinder Bot! 🙂"
        ,
        parse_mode=ParseMode.HTML, 
        reply_markup=types.ReplyKeyboardRemove()
    )


@command_router.message(Command('cancel'))
@command_router.message(F.text.casefold() == "cancel")
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