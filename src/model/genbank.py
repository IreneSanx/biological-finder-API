import requests
from xml.dom import minidom
from collections import defaultdict
from functools import lru_cache
import rdflib
from rdflib import RDF, Graph, Literal

import re
from abc import ABC, abstractmethod
from typing import Dict
import httpx
from httpx import Response

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

RDF_REPOSITORY_ENDPOINT = "http://localhost:8890/sparql-auth/"
RDF_REPOSITORY_USERNAME = "dba"
RDF_REPOSITORY_PASSWORD = "dba"
RDF_REPOSITORY_DB = "http://localhost:8890/GENBANK"

rdf_connection={
        "database": RDF_REPOSITORY_DB,
        "endpoint": RDF_REPOSITORY_ENDPOINT,
        "username": RDF_REPOSITORY_USERNAME,
        "password": RDF_REPOSITORY_PASSWORD,
}

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

        return res

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

def _etree_to_dict(t):
    """ Transform an XML to a Python dictionary. """
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(_etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d

def xml_to_dict(archi):
    root = ET.parse(archi).getroot()
    root_as_dict = _etree_to_dict(root)
    diccionario = {}
    for conjunto in root_as_dict: 
        sequence = root_as_dict[conjunto]
        for element in sequence: 
            seq = sequence[element]
            for seq_entry in seq:
                seq_entry_ = seq[seq_entry]
                if "Seq-entry_set" in seq_entry_:
                    for seq_entry_set in seq_entry_:
                        seq_entry_set_ = seq_entry_[seq_entry_set]
                        for bioseq_set in seq_entry_set_:
                            bioseq_set_ = seq_entry_set_[bioseq_set]
                            bioseq_set_seq_set = bioseq_set_["Bioseq-set_seq-set"]
                            for seq_entry in bioseq_set_seq_set: 
                                seq_entry_ = bioseq_set_seq_set[seq_entry]
                                seq_entry_seq = seq_entry_[0]
                                for bioseq in seq_entry_seq:
                                    bioseq_ = seq_entry_seq[bioseq]
                                    for bioseq_id in bioseq_:
                                        bioseq_id_ = bioseq_[bioseq_id]
                                        bioseq_id_id = bioseq_id_["Bioseq_id"]
                                        for seq_id in bioseq_id_id: 
                                            seq_id_ = bioseq_id_id[seq_id]
                                            seq_id_gi = seq_id_[1]
                                            diccionario["ID"] = seq_id_gi["Seq-id_gi"]
                                            seq_id_embl = seq_id_[0]
                                            for textseq in seq_id_embl:
                                                textseq_id = seq_id_embl[textseq]
                                                textseq_id_ = textseq_id["Textseq-id"]
                                                diccionario["Locus"] = textseq_id_["Textseq-id_accession"]
                                        bioseq_descr = bioseq_id_["Bioseq_descr"]
                                        for seq_descr in bioseq_descr:
                                            seq_descr_ = bioseq_descr[seq_descr]
                                            seqdesc = seq_descr_["Seqdesc"]
                                            seqdesc_title = seqdesc[0]
                                            diccionario["Description"] = seqdesc_title["Seqdesc_title"]
                            bioseq_set_descr = bioseq_set_["Bioseq-set_descr"]   
                            for seqdesc in bioseq_set_descr: 
                                seqdesc_ = bioseq_set_descr[seqdesc]
                                for seqdesc_desc in seqdesc_:
                                    seq_desc = seqdesc_[seqdesc_desc]
                                    seqdesc_source = seq_desc[0]
                                    for biosource in seqdesc_source: 
                                        biosource_ = seqdesc_source[biosource]
                                        for biosource_source in biosource_: 
                                            biosour = biosource_[biosource_source]
                                            biosource_org = biosour["BioSource_org"]
                                            for org_ref in biosource_org: 
                                                orgref = biosource_org[org_ref]
                                                diccionario["Organism"] = orgref["Org-ref_taxname"]
                                                org_ref_orgname = orgref["Org-ref_orgname"]
                                                for orgname in org_ref_orgname:
                                                    orgname_ = org_ref_orgname[orgname]
                                                    diccionario["Taxonomy"] = orgname_["OrgName_lineage"]
                                        bioseq_annot = bioseq_id_["Bioseq_annot"]
                                        for seq_annot in bioseq_annot:
                                            seq_annot_ = bioseq_annot[seq_annot]
                                            for seq_annot_data in seq_annot_:
                                                seq_annot_data_ = seq_annot_[seq_annot_data]
                                                seq_annot_data_ftable = seq_annot_data_["Seq-annot_data_ftable"]
                                                for seq_feat in seq_annot_data_ftable:
                                                    seq_feat_ = seq_annot_data_ftable[seq_feat]
                                                    seq_feat_data = seq_feat_[0]
                                                    seq_feat_data_ = seq_feat_data["Seq-feat_data"]
                                                    for seqfeatdata in seq_feat_data_:
                                                        seqfeatdatagene = seq_feat_data_[seqfeatdata]
                                                        if "SeqFeatData_gene" in seqfeatdatagene:
                                                            seqfeatdata_gene = seqfeatdatagene["SeqFeatData_gene"]
                                                            generef = seqfeatdata_gene["Gene-ref"]
                                                            diccionario["Name"] = generef["Gene-ref_locus"]
                                                        else:
                                                            seq_feat_data1 = seq_feat_[1]
                                                            seq_feat_data_1 = seq_feat_data1["Seq-feat_data"]
                                                            for seqfeatdata1 in seq_feat_data_1:
                                                                seqfeatdatagene1 = seq_feat_data_1[seqfeatdata1]
                                                                if "SeqFeatData_gene" in seqfeatdatagene1:
                                                                    seqfeatdata_gene1 = seqfeatdatagene1["SeqFeatData_gene"]
                                                                    generef1 = seqfeatdata_gene1["Gene-ref"]
                                                                    diccionario["Name"] = generef1["Gene-ref_locus"]
                                                                else:
                                                                    if len(seq_feat_) >= 3:
                                                                        seq_feat_data2 = seq_feat_[2]
                                                                        seq_feat_data_2 = seq_feat_data2["Seq-feat_data"]
                                                                        for seqfeatdata2 in seq_feat_data_2:
                                                                            seqfeatdatagene2 = seq_feat_data_2[seqfeatdata2]
                                                                            if "SeqFeatData_gene" in seqfeatdatagene2:
                                                                                seqfeatdata_gene2 = seqfeatdatagene2["SeqFeatData_gene"]
                                                                                generef2 = seqfeatdata_gene2["Gene-ref"]
                                                                                diccionario["Name"] = generef2["Gene-ref_locus"]
                                                                            else:
                                                                                diccionario["Name"] = "No name given"
                else:
                    for seq_entry_seq in seq_entry_:
                        seq_entry_seq_ = seq_entry_[seq_entry_seq]
                        for bioseq in seq_entry_seq_:
                            bioseq_ = seq_entry_seq_[bioseq]
                            bioseq_id = bioseq_["Bioseq_id"]
                            seq_id = bioseq_id["Seq-id"]
                            seq_id_gi = seq_id[1]
                            if "Seq-id_gi" in seq_id_gi:
                                diccionario["ID"] = seq_id_gi["Seq-id_gi"]
                            else: 
                                if len(seq_id) > 2:
                                    seq_id_gi = seq_id[2]
                                    if "Seq-id_gi" in seq_id_gi:
                                        diccionario["ID"] = seq_id_gi["Seq-id_gi"]
                                    else:
                                        diccionario["ID"] = "No ID given"
                            seq_id_other = seq_id[0]
                            seq_id_other_ = seq_id_other["Seq-id_other"]
                            textseq_id = seq_id_other_["Textseq-id"]
                            textseq_accession = textseq_id["Textseq-id_accession"]
                            diccionario["Locus"] = textseq_accession
                            bioseq_descr = bioseq_["Bioseq_descr"]
                            seq_descr = bioseq_descr["Seq-descr"]
                            seqdesc = seq_descr["Seqdesc"]
                            seqdesc_title = seqdesc[12]
                            if "Seqdesc_title" in seqdesc_title:
                                diccionario["Description"] = seqdesc_title["Seqdesc_title"]
                            else:
                                diccionario["Description"] = "No title given"
                            seqdesc_source = seqdesc[0]
                            seqdescsource = seqdesc_source["Seqdesc_source"]
                            biosource = seqdescsource["BioSource"]
                            biosource_org = biosource["BioSource_org"]
                            org_ref = biosource_org["Org-ref"]
                            diccionario["Organism"] = org_ref["Org-ref_taxname"]
                            org_ref_orgname = org_ref["Org-ref_orgname"]
                            orgname = org_ref_orgname["OrgName"]
                            diccionario["Taxonomy"] = orgname["OrgName_lineage"]
                            diccionario["Name"] = "No name given"
    return diccionario


def xml_to_dict_name(archi):
    root = ET.parse(archi).getroot()
    root_as_dict = _etree_to_dict(root)
    diccionario = {}
    for conjunto in root_as_dict: 
        sequence = root_as_dict[conjunto]
        for entrez in sequence: 
            entrezgene_ = sequence[entrez]
            entreztrack = entrezgene_["Entrezgene_track-info"]
            entrezlocus = entrezgene_["Entrezgene_locus"]
            entrezgene = entrezgene_["Entrezgene_gene"]
            entrezsource = entrezgene_["Entrezgene_source"]
            diccionario["ID"] = entreztrack["Gene-track"]["Gene-track_geneid"]
            commentary = entrezlocus["Gene-commentary"]
            if "Gene-commentary_accession" in commentary:
                diccionario["Locus"] = commentary["Gene-commentary_accession"]
            else:
                diccionario["Locus"] = "No locus given"
            gen_ref = entrezgene["Gene-ref"]
            if "Gene-ref_desc" in gen_ref:
                diccionario["Description"] = entrezgene["Gene-ref"]["Gene-ref_desc"]
            else:
                diccionario["Description"] = "No description given"
            diccionario["Organism"] = entrezsource["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_taxname"]
            diccionario["Taxonomy"] = entrezsource["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_orgname"]["OrgName"]["OrgName_lineage"]
            diccionario["Name"] = entrezgene["Gene-ref"]["Gene-ref_locus"]
            
    return diccionario

def lista_ids(archi):
    root = ET.parse(archi).getroot()
    root_as_dict = _etree_to_dict(root)
    for conjunto in root_as_dict: 
        sequence = root_as_dict[conjunto]
        idList = sequence["IdList"]
        ids = idList["Id"]
    return ids

def descarga_ids_name(listaIds):
    lista = []
    for i in range(0,len(listaIds)):
        id = listaIds[i]
        busq = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=' + id + '&rettype=null&retmode=xml'
        myfile = requests.get(busq)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + id +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'wb').write(myfile.content)
        resultado = xml_to_dict_name('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml')
        lista.append(resultado)
    return lista

def descarga_ids(listaIds):
    lista = []
    for i in range(0,len(listaIds)):
        id = listaIds[i]
        busq = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + id + '&rettype=native&retmode=xml'
        myfile = requests.get(busq)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + id +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'wb').write(myfile.content)
        resultado = xml_to_dict('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml')
        lista.append(resultado)
    return lista

def dic_to_rdf(diccionario):
    """
    Transform a JSON-formatted workflow to RDF n-triples.
    """

    PREFIX_GENBANK = "http://www.ontologies.khaos.uma.es/GENBANK/"

    uriref_gene = rdflib.URIRef(PREFIX_GENBANK + "Gene")
    uriref_ID = rdflib.URIRef(PREFIX_GENBANK + "ID")
    uriref_locus = rdflib.URIRef(PREFIX_GENBANK + "Locus")
    uriref_description = rdflib.URIRef(PREFIX_GENBANK + "Description")
    uriref_organism = rdflib.URIRef(PREFIX_GENBANK + "Organism")
    uriref_taxonomy = rdflib.URIRef(PREFIX_GENBANK + "Taxonomy")
    uriref_geneName = rdflib.URIRef(PREFIX_GENBANK + "Name")

    g = Graph()

    name = rdflib.URIRef(PREFIX_GENBANK + "Gene_" + diccionario["ID"])
    g.add((name, RDF.type, uriref_gene))
    g.add((name, uriref_ID, Literal(diccionario["ID"])))
    g.add((name, uriref_locus, Literal(diccionario["Locus"])))
    g.add((name, uriref_description, Literal(diccionario["Description"])))
    g.add((name, uriref_organism, Literal(diccionario["Organism"])))
    g.add((name, uriref_taxonomy, Literal(diccionario["Taxonomy"])))
    g.add((name, uriref_geneName, Literal(diccionario["Name"])))
    
    return g

def printResultados(resu):
	list_ = []
	for elemen in resu: 
		claves = elemen.keys()
		for item in claves: 
			elemen[item] = elemen.get(item)["value"]
		list_.append(elemen)
	return list_

def build_query(id,locus):
    consulta = """
        PREFIX gen: <http://www.ontologies.khaos.uma.es/GENBANK/>
    """
    consulta += "select \n"
    consulta += "  ?ID\n"
    consulta += "  ?Locus\n"
    consulta += "  ?Description\n"
    consulta += "  ?Organism\n"
    consulta += "  ?Taxonomy\n"
    consulta += "  ?Name\n"
    consulta += "where{{\n"

    consulta += " ?gen rdf:type gen:Gene.\n"
    consulta += " ?gen gen:ID ?ID.\n"
    if (id != ''):
        consulta += " filter(regex(lcase(?ID), " + '"' + id + '"' + ")).\n"
    consulta += "\n"

    consulta += " ?gen gen:Locus ?Locus .\n"
    if (locus != ''):
        consulta += " filter(regex(lcase(?Locus), " + '"' + locus + '"' + ")).\n"
    consulta += "\n"

    consulta += " ?gen gen:Description ?Description.\n"
    consulta += " ?gen gen:Organism ?Organism.\n"
    consulta += " ?gen gen:Taxonomy ?Taxonomy.\n"
    consulta += " ?gen gen:Name ?Name.\n"
    consulta += "\n"
    
    consulta += "}}\n"

    return consulta

def ejecutarConsulta(id,locus):
    store = Virtuoso(**rdf_connection)
    consu = build_query(id,locus)
    print(consu)
    resu = store.query(consu)
    return resu

def insertElement(diccionario):
    rdf_ = dic_to_rdf(diccionario)
    rdf_as_nt = rdf_.serialize(format="nt")
    triples = str(rdf_as_nt.decode("UTF-8"))
    store = Virtuoso(**rdf_connection)
    query = "INSERT DATA { GRAPH <" + store.database + "> {" + triples + "} }"
    store.update(query)

@lru_cache
def consultaGenBank(id,locus,name,organism,description):
    if (id != ""): 
        lista = []
        """consulta_en_GenBank = ejecutarConsulta(id,"")
        if consulta_en_GenBank != []:
            lista.append(consulta_en_GenBank)
        else:
            urlID = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + id + '&rettype=native&retmode=xml'
            myfileID = requests.get(urlID)
            open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'w').close
            open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'wb').write(myfileID.content)
            resuID = xml_to_dict('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + id +'.xml')
            insertElement(resuID)
            lista.append(resuID)"""
        urlID = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + id + '&rettype=native&retmode=xml'
        myfileID = requests.get(urlID)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ id +'.xml', 'wb').write(myfileID.content)
        resuID = xml_to_dict('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + id +'.xml')
        lista.append(resuID)
    elif(locus != ""):
        lista = []
        """consulta_en_GenBank = ejecutarConsulta("",locus)
        if consulta_en_GenBank != []:
            lista.append(consulta_en_GenBank)
        else:
            urlLocus = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + locus + '&rettype=native&retmode=xml'
            myfileLocus = requests.get(urlLocus)
            open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ locus +'.xml', 'w').close
            open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ locus +'.xml', 'wb').write(myfileLocus.content)
            resuLocus = xml_to_dict('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + locus +'.xml')
            insertElement(resuLocus)
            lista.append(resuLocus)"""
        urlLocus = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + locus + '&rettype=native&retmode=xml'
        myfileLocus = requests.get(urlLocus)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ locus +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ locus +'.xml', 'wb').write(myfileLocus.content)
        resuLocus = xml_to_dict('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/' + locus +'.xml')
        lista.append(resuLocus)
    elif(name != ""):
        urlSearchName= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=' + name + '[Gene Name]&rettype=gb&retmode=xml'
        myfileSearchName = requests.get(urlSearchName)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ name +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ name +'.xml', 'wb').write(myfileSearchName.content)
        listaIds = lista_ids('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ name +'.xml')
        lista = descarga_ids_name(listaIds)
    elif (organism != ""):
        urlSearchOrganism = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term='+ organism +'[organism]&rettype=gb&retmode=xml'
        myfileSearchOrganism = requests.get(urlSearchOrganism)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ organism +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ organism +'.xml', 'wb').write(myfileSearchOrganism.content)
        listaIds = lista_ids('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ organism +'.xml')
        lista = descarga_ids(listaIds)
    elif (description != ""):
        urlSearchDescription = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term='+ description +'[title]&rettype=gb&retmode=xml'
        myfileSearchDescription = requests.get(urlSearchDescription)
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ description +'.xml', 'w').close
        open('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ description +'.xml', 'wb').write(myfileSearchDescription.content)
        listaIds = lista_ids('/Users/irenesanchez/Desktop/UNI/TFG/Virtuoso/src/'+ description +'.xml')
        lista = descarga_ids(listaIds)
    return lista
