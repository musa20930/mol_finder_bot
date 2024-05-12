from loguru import logger

from aiogram import types, Router
from aiogram.filters.command import Command, CommandStart
from aiogram.fsm.context import FSMContext
from aiogram import F


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
        text="Coming soon...",
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