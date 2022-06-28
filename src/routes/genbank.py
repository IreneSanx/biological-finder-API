from typing import List

from fastapi import APIRouter, Depends, HTTPException
from pydantic import ValidationError
from model.genbank import consultaGenBank

router = APIRouter()

@router.get("/return_query_genbank", name="Devolver la consulta realizada", tags=["consultas"])
async def return_query_genbank(gene_id: str = "", gene_locus: str = "", gene_name: str = "", gene_organism: str = "", gene_description: str = ""):
    RQG = None
    try:
        RQG = consultaGenBank(id = gene_id, locus = gene_locus, name = gene_name, organism = gene_organism, description = gene_description)
    except:
        print("Something was wrong")
    return RQG
