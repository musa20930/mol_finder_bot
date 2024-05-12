import logging
import sys
from loguru import logger


class InterceptHandler(logging.Handler):
    """Class for intercepting log messages to use loguru instead of logging."""
    def emit(self, record) -> None:
        level = logger.level(record.levelname).name
        logger.log(level, record.getMessage())


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