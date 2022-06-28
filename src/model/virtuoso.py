from pathlib import Path

from pydantic import AnyUrl, BaseSettings

import re
from abc import ABC, abstractmethod
from typing import Dict
from functools import lru_cache

import httpx
from httpx import Response

COMMON_PREFIXES = """
	PREFIX umls: <http://bioportal.bioontology.org/ontologies/umls/>
    PREFIX omim: <http://purl.bioontology.org/ontology/OMIM/>
"""

DISEASE = "disease or syndrome"   

def buildQuery(id, name, geneSymbol, geneLocus, symptom):
    consulta = COMMON_PREFIXES
    consulta += "select distinct \n"
    consulta += "  ?id\n"
    consulta += "  ?name\n"
    consulta += "  ?geneSymbol\n"
    consulta += "  ?geneLocus\n" 
    consulta += "  ?symptom\n" 
    consulta += "where{{\n"

    consulta += " ?indv umls:hasSTY ?sty .\n"
    consulta += " ?sty skos:prefLabel ?labelSty.\n"
    consulta += " filter(regex(lcase(?labelSty)," + '"' + DISEASE + '"' + ")).\n"
    
    consulta += " ?indv skos:notation ?id.\n"
    if (id != ''):
        consulta += " filter(regex(lcase(?id), " + '"' + id + '"' + ")).\n"
    consulta += "\n"
    
    consulta += " ?indv skos:prefLabel ?name.\n"
    if (name != ''):
        consulta += " filter(regex(lcase(?name), " + '"' + name + '"' + ",\"i\" )).\n"
    consulta += "\n"

    consulta += " ?indv omim:has_manifestation ?hasMan.\n"
    consulta += " ?hasMan skos:prefLabel ?symptom.\n"
    if (symptom != ''):
        consulta += " filter(regex(lcase(?symptom), " + '"' + symptom + '"' + ",\"i\" )).\n"
    consulta += "\n"

    consulta += " ?indv omim:GENELOCUS ?geneLocus.\n" 
    if (symptom != ''):
        consulta += " filter(regex(lcase(?geneLocus), " + '"' + geneLocus + '"' + ",\"i\" )).\n"
    consulta += "\n"

    consulta += " ?indv omim:GENESYMBOL ?geneSymbol.\n"
    if (symptom != ''):
        consulta += " filter(regex(lcase(?geneSymbol), " + '"' + geneSymbol + '"' + ",\"i\" )).\n"
    consulta += "\n"
    
    consulta += "}}\n"

    return consulta

def printResultados(resu):
	list_ = []
	for elemen in resu: 
		claves = elemen.keys()
		for item in claves: 
			elemen[item] = elemen.get(item)["value"]
		list_.append(elemen)
	return list_

def combinarDiccionarios_(listaDic):
    listaIDs = []
    listaFinal = []
    for i in range(0, len(listaDic)):
        if listaDic[i]["id"] not in listaIDs:
            listaIDs.append(listaDic[i]["id"])
    for id in listaIDs:
        symptom = []
        genes = []
        for j in range(0, len(listaDic)):
            if listaDic[j]["id"] == id and not listaDic[j]["symptom"] in symptom: 
                if symptom == []:
                    symptom.append(listaDic[j]["symptom"])
                else:
                    symptom.append(", ")
                    symptom.append(listaDic[j]["symptom"])
            if not listaDic[j]["geneSymbol"] in genes: 
                genes.append(listaDic[j]["geneSymbol"])
        for g in range(0, len(listaDic)):
            if listaDic[g]["id"] == id:
                listaDic[g]["symptom"] = symptom
        for y in range(0, len(listaDic)):
            if listaDic[y]["id"] == id and listaDic[y]["geneSymbol"] in genes: 
                listaFinal.append(listaDic[y])
                genes.remove(listaDic[y]["geneSymbol"])
    return listaFinal

@lru_cache
def ejecutarConsulta(id, name, geneSymbol, geneLocus, symptom):

    # ontology RDF repository connection
    RDF_REPOSITORY_ENDPOINT = "http://localhost:8890/sparql-auth/"
    RDF_REPOSITORY_USERNAME = "dba"
    RDF_REPOSITORY_PASSWORD = "dba"
    RDF_REPOSITORY_DB = "http://localhost:8890/OMIM"

    rdf_connection={
            "database": RDF_REPOSITORY_DB,
            "endpoint": RDF_REPOSITORY_ENDPOINT,
            "username": RDF_REPOSITORY_USERNAME,
            "password": RDF_REPOSITORY_PASSWORD,
    }
    store = Virtuoso(**rdf_connection)

    consu = buildQuery(id, name, geneSymbol, geneLocus, symptom)

    print(consu)
    
    resu = store.query(consu)

    return resu


class RDFRepository(ABC):
    def __init__(self, endpoint: str, database: str, username: str = None, password: str = None):
        self.endpoint = endpoint
        self.database = database
        self.username = username
        self.password = password

    @abstractmethod
    def query(self, query: str, **params) -> dict:
        """
        Performs a query against the database.
        """
        pass

    @abstractmethod
    def update(self, query: str) -> None:
        """
        Performs an update query against the database.
        """
        pass

class Virtuoso(RDFRepository):
    """
     SPARQLWrapper for Virtuoso graph store.
    """

    COMMENTS_PATTERN = re.compile(r"(^|\n)\s*#.*?\n")

    def __init__(self, endpoint: str, database: str, username: str = None, password: str = None):
        super().__init__(endpoint, database, username, password)

        self.parameters: Dict[str, str] = {}
        self.headers: Dict[str, str] = {}

        self._setup_request()

    def _setup_request(self) -> None:
        self._add_parameter("default-graph-uri", self.database)
        self._add_parameter("User-Agent", "salon")

        self._add_header(
            "Accept",
            "application/sparql-results+json,application/json,text/javascript,application/javascript",
        )

    def query(self, query: str, **params) -> dict:
        """
        Run 'SELECT' query with http Auth DIGEST and return results in JSON format.
        Protocol details at http://www.w3.org/TR/sparql11-protocol/#query-operation
        """
        query_string = query.format(**params)
        query_string = re.sub(self.COMMENTS_PATTERN, "\n\n", query_string)
        print(query_string)

        self._add_parameter("query", query_string)
        req =  self._get()

        # convert to json and return bindings
        result = {}
        if not req.is_error:
            result = req.json()
            result = result["results"]["bindings"]

        self._remove_parameter("query")

        res = printResultados(result)

        resultadoFinal = combinarDiccionarios_(res)

        return resultadoFinal

    def update(self, query: str) -> None:
        """
        Run 'INSERT' update query with http Auth DIGEST.
        Protocol details at http://www.w3.org/TR/sparql11-protocol/#update-operation
        """
        self._add_header("Content-Type", "application/sparql-update")

        req =  self._post_directly(query)

        # convert to json and return bindings
        result = {}
        if not req.is_error:
            result = req.json()
            result = result["results"]["bindings"]

        self._remove_header("Content-Type")

        return result

    def _post_directly(self, query: str, **kwargs) -> Response:
        auth = httpx.DigestAuth(self.username, self.password)
        with httpx.Client(timeout=12000) as client:
            req =  client.post(
                self.endpoint,
                data=query,
                params=self.parameters,
                headers=self.headers,
                auth=auth,
            )
        if req.is_error:
            print(req.text, req.status_code)
        return req

    def _get(self, **kwargs) -> Response:
        auth = httpx.DigestAuth(self.username, self.password)
        with httpx.Client() as client:
            req =  client.get(
                self.endpoint,
                params=self.parameters,
                headers=self.headers,
                auth=auth,
            )
        return req

    def _add_header(self, param: str, value: str) -> None:
        """
        Adds new custom header to request.
        """
        self.headers[param] = value

    def _remove_header(self, param: str) -> None:
        """
        Deletes header from request.
        """
        try:
            del self.headers[param]
        except KeyError:
            pass

    def _add_parameter(self, param: str, value: str) -> None:
        """
        Adds new parameter to request.
        """
        self.parameters[param] = value

    def _remove_parameter(self, param: str) -> None:
        """
        Deletes parameter from request.
        """
        try:
            del self.parameters[param]
        except KeyError:
            pass


