from typing import Optional
from fastapi import FastAPI
from routes.virtuoso import router as virtuoso_router
from routes.uniprot import router as uniprot_router
from routes.genbank import router as genbank_router
from starlette.middleware.cors import CORSMiddleware

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware, allow_origins=origins, allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"Hello": "World"}

@app.get("/items/{item_id}")
def read_item(item_id: int, q: Optional[str] = None):
    return {"item_id": item_id, "q": q}

app.include_router(virtuoso_router, prefix="/virtuoso/api")
app.include_router(uniprot_router, prefix="/uniprot/api")
app.include_router(genbank_router, prefix="/genbank/api")