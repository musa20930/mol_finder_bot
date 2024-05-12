"""
- To run with debug mode use command:
poetry run python bot.py --debug

- If you don't want to see debug log level:
poetry run python bot.py
"""

import asyncio
import os

from aiogram import Bot, Dispatcher, types

from handlers import commands, search, log_handler


async def main() -> None:
    log_handler.add_logging()

    API_TOKEN = os.getenv('MOL_FINDER_BOT_TOKEN')
    bot = Bot(token=API_TOKEN)
    dp = Dispatcher() # отслеживает все входящие события
    # Add command menu with set_my_commands
    await bot.set_my_commands([
        types.BotCommand(command="search", description="Search molecule"),
        types.BotCommand(command="start", description="Start bot"),
        types.BotCommand(command="cancel", description="Cancel running task"),
        types.BotCommand(command="help", description="Manual"),
    ])
    # routers useful to chain questions and write necessary data
    dp.include_routers(search.search_router, commands.command_router)
    await dp.start_polling(bot)


if __name__ == '__main__':    
    asyncio.run(main())

