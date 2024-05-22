from aiogram.fsm.state import State, StatesGroup


class SearchInfo(StatesGroup):
    molecule_name = State()
    chembl_id = State()
    inchi_key = State()
    smiles = State()
    mol_info = State()
    mol_descriptors = State()
    
    search_multiple = State()
    mol_series_info = State()
    similarity_percent = State()
    molecule_name_sim = State()
    chembl_id_sim = State()
    inchi_key_sim = State()
    smiles_sim = State()

    mol_mass = State()
    mol_mass_comparison = State()
    mol_logp = State()
    mol_logp_comparison = State()