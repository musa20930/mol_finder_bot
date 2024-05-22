from aiogram.types import (
    InlineKeyboardButton,
    InlineKeyboardMarkup
)

from dataclasses import dataclass


@dataclass
class Button:
    search_by_id = InlineKeyboardButton(text='Search by ChEMBL🆔', callback_data='chembl_id')
    search_by_name = InlineKeyboardButton(text='Search by Name', callback_data='name')
    search_by_inchi = InlineKeyboardButton(text='Search by INCHI 🔑', callback_data='inchi')
    search_by_smiles = InlineKeyboardButton(text='Search by SMILES', callback_data='smiles')
    get_mol_structure = InlineKeyboardButton(text='Get Molecule Structure🧬', callback_data='structure')
    search_by_similarity = InlineKeyboardButton(text='Search by Similarity❄️❄️', callback_data='similarity')
    calc_properties = InlineKeyboardButton(text='Calculate Properties', callback_data='properties')
    standardize = InlineKeyboardButton(text='Standardize Molecule', callback_data='standardize')

    save_as = InlineKeyboardButton(text='Save As...')
    save_as_csv = InlineKeyboardButton(text='CSV', callback_data='properties_csv')
    save_as_smi = InlineKeyboardButton(text='SMI', callback_data='smi')
    save_as_sdf = InlineKeyboardButton(text='SDF', callback_data='sdf')

    by_mass = InlineKeyboardButton(text='By Mass (g/mol)', callback_data='mass')
    by_logp = InlineKeyboardButton(text='By LogP', callback_data='logp')
    by_lipinski = InlineKeyboardButton(text='By Lipinski Rule5️⃣', callback_data='lipinski')
    drugs_by_date = InlineKeyboardButton(text='Approved Drugs by 📅', callback_data='drugs')
    by_similarity = InlineKeyboardButton(text='By Similarity❄️❄️', callback_data='similarity')
    by_name = InlineKeyboardButton(text='By Name', callback_data='name')
    by_id = InlineKeyboardButton(text='By ChEMBL🆔', callback_data='chembl_id')
    by_inchi = InlineKeyboardButton(text='By INCHI 🔑', callback_data='inchi')
    by_smiles = InlineKeyboardButton(text='By SMILES', callback_data='smiles')

    by_name_sim = InlineKeyboardButton(text='By Name', callback_data='name_sim')
    by_id_sim = InlineKeyboardButton(text='By ChEMBL🆔', callback_data='chembl_id_sim')
    by_inchi_sim = InlineKeyboardButton(text='By INCHI 🔑', callback_data='inchi_sim')
    by_smiles_sim = InlineKeyboardButton(text='By SMILES', callback_data='smiles_sim')


class Keyboard:
    @property
    def after_search(self):     # after finding molecule
        keyboard = [
            [Button.calc_properties],
            [Button.standardize],
            [Button.get_mol_structure],
            [Button.search_by_similarity],
        ]
        return InlineKeyboardMarkup(inline_keyboard=keyboard, 
                                    input_field_placeholder="Choose option")

    @property
    def search(self):
        keyboard = [
            [Button.by_mass, 
            Button.by_lipinski],
            [Button.by_logp,
            Button.drugs_by_date],
            [Button.by_similarity], 
            [Button.by_name, 
            Button.by_id],
            [Button.by_inchi,
            Button.by_smiles],
        ]
        return InlineKeyboardMarkup(inline_keyboard=keyboard, 
                                    input_field_placeholder="Choose option")

    @property
    def similarity_search(self):
        keyboard = [
            [Button.by_name_sim,
            Button.by_id_sim],
            [Button.by_inchi_sim,
            Button.by_smiles_sim],
        ]
        return InlineKeyboardMarkup(inline_keyboard=keyboard, 
                                    input_field_placeholder="Choose option")
    
    @property
    def properties(self):
        keyboard = [[Button.save_as_csv]]
        return InlineKeyboardMarkup(inline_keyboard=keyboard)

