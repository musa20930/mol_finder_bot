from aiogram.fsm.state import State, StatesGroup


class SearchInfo(StatesGroup):
    search_type = State()
    molecule_name = State()
    chembl_id = State()
    molecule_next_step = State()
    mol_info = State()
    similarity_percent = State()