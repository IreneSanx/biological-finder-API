from typing import List

from fastapi import APIRouter, Depends, HTTPException
from pydantic import ValidationError
from model.virtuoso import ejecutarConsulta

router = APIRouter()

@router.get("/return_query_omim", name="Devolver la consulta realizada", tags=["consultas"])
async def return_query_omim(disease_id: str = "", disease_name: str = "", disease_geneSymbol: str = "", disease_geneLocus: str = "", disease_symptom: str = ""):
    RQ = None
    try:
        RQ = ejecutarConsulta(id = disease_id, name = disease_name, geneSymbol = disease_geneSymbol, geneLocus = disease_geneLocus, symptom = disease_symptom)
    except:
        print("Something was wrong")
    
    return RQ
