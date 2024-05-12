from pydantic import BaseModel, Field, field_validator, ValidationInfo
from rdkit import Chem


class SmilesSimilarity(BaseModel):
    smiles_percent: str = Field(pattern=r"^(\w+)\s(\d+)$")


class Smiles(BaseModel):
    smiles: str = None

    @field_validator('smiles')
    def validate_smiles(smi: str, info: ValidationInfo) -> str:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            raise TypeError(f'Invalid smiles string: {smi}')
        return smi