from loguru import logger
import json
from io import StringIO, BytesIO
# import subprocess

from aiogram.fsm.context import FSMContext
from aiogram.enums import ParseMode
from aiogram.types import input_file
from aiogram.types import (
    Message,
    FSInputFile
)

from keyboards import Keyboard

from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
# from chython.files import smiles as chy_smiles
# from CGRtools.files import SMILESRead
# from svglib.svglib import svg2rlg
# from reportlab.graphics import renderPM


async def display_mol_info(message: Message, state: FSMContext, result: list) -> None:
    await state.update_data(mol_info=result[0])
    pref_name = result[0]['pref_name']
    await message.answer(
        f"<b>Name</b>: {pref_name if pref_name is not None else 'not provided'} \n"
        f"<b>ID</b>: {result[0]['molecule_chembl_id']} \n"
        f"<b>SMILES</b>: {result[0]['molecule_structures']['canonical_smiles']}",
        parse_mode=ParseMode.HTML
    )
    await message.answer(
        text='Choose option below or rerun /search', 
        reply_markup=Keyboard().after_search
    )


async def display_mol_info_for_sim(message: Message, state: FSMContext, result: list) -> None:
    await state.update_data(mol_info=result[0])
    pref_name = result[0]['pref_name']
    await message.answer(
        f"<b>Name</b>: {pref_name if pref_name is not None else 'not provided'} \n"
        f"<b>ID</b>: {result[0]['molecule_chembl_id']} \n"
        f"<b>SMILES</b>: {result[0]['molecule_structures']['canonical_smiles']}",
        parse_mode=ParseMode.HTML
    )


async def send_structure_png(mol_smiles: str, chembl_id):
    # TODO: chython to generate or cgrtools, cairo or inkscape to convert 
    logger.debug(f'Drawing structure of \'{mol_smiles}\'')
    pl = Chem.MolFromSmiles(mol_smiles)
    pl_image = Draw.MolToImage(pl, size=(1000, 1000))
    pl_image.save('assets/image.png')
    mol_img = FSInputFile('assets/image.png')
    return mol_img


async def send_structure_svg(mol_smiles: str, chembl_id):
    logger.debug(f'Drawing structure of \'{mol_smiles}\'')
    image = new_client.image
    image.set_format('svg')
    mol_img = image.get(chembl_id)
    # mol_img = await convert_svg_to_png(mol_img)

    file = StringIO()
    file.write(mol_img)
    file.seek(0)
    mol_svg = input_file.BufferedInputFile(BytesIO(file.read().encode('utf-8')).getbuffer(), 
                                           filename=f"{chembl_id}.svg")
    return mol_svg


# async def convert_svg_to_png(svg_str: str):
#     width, height = 1000, 1000
#     inkscape = '/Applications/Inkscape.app' # path to inkscape executable
#     # svg string -> png data
#     result = subprocess.run([inkscape, '--export-type=png', '--export-filename=-', f'--export-width={width}', f'--export-height={height}', '--pipe'], input=svg_str.encode(), capture_output=True)
#     #   (result.stdout will have the png data)
#     return result.stdout


async def save_to_smi(mol_series_smiles: str, chembl_id):
    logger.debug(f'Drawing structure of \'{mol_series_smiles}\'')
    file = StringIO()
    for i, mol_smiles in enumerate(mol_series_smiles):
        file.write(f'{mol_smiles}\n')
    file.seek(0)
    smi_file = input_file.BufferedInputFile(BytesIO(file.read().encode('utf-8')).getbuffer(), 
                                            filename=f"{chembl_id}.smi")
    return smi_file


async def save_to_sdf(mol_series_file: str, chembl_id):
    logger.debug(f'Saving to sdf \'{mol_series_file}\'')
    file = StringIO()
    for i, mol_file in enumerate(mol_series_file):
        file.write(f'{chembl_id[i]}.sdf')
        file.write(mol_file)
        file.write('\n$$$$\n')
    file.seek(0)
    sdf_file = input_file.BufferedInputFile(BytesIO(file.read().encode('utf-8')).getbuffer(), 
                                            filename=f"{chembl_id}.sdf")
    return sdf_file
