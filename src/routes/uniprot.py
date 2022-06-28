from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException
from pydantic import ValidationError
from model.uniprot import sparqlwrap

router = APIRouter()

@router.get("/return_query_uniprot", name="Devolver la consulta realizada", tags=["consultas"])
async def return_query_uniprot(protein_id: str = "", protein_name: str = "", protein_gene: str = "", protein_organism: str = "",
    protein_disease: str = "", protein_similarity: str = "", protein_location: str = "", protein_function: str = "", protein_pharmaceutical:str = ""):
    RQU = None
    try:
        RQU = sparqlwrap(proteinId = protein_id, proteinName = protein_name, geneName = protein_gene, organismName = protein_organism, diseaseAnnotation = protein_disease, domainName = "", similarityAnnotation = protein_similarity, locationAnnotation = protein_location, functionAnnotation = protein_function, pharmaceuticalAnnotation = protein_pharmaceutical)
    except:
        print("Something was wrong")
    return RQU