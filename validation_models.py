from pydantic import BaseModel, Field, field_validator, ValidationInfo
from typing import Union
from rdkit import Chem


class Smiles(BaseModel):
    smiles: Union[str, None] = None

    @field_validator('smiles')
    def validate_smiles(smi: Union[str, None], info: ValidationInfo) -> str:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            raise TypeError(f'Invalid smiles string: {smi}')
        return smi