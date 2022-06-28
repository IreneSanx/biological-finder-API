from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import json
from functools import lru_cache

# Constante con los prefijos comunes a usar en queries
COMMON_PREFIXES = """
	PREFIX up: <http://purl.uniprot.org/core/>
	PREFIX keywords: <http://purl.uniprot.org/keywords/>
	PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
	PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
	PREFIX ec: <http://purl.uniprot.org/enzyme/>
	PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
	PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
	PREFIX owl: <http://www.w3.org/2002/07/owl#>
	PREFIX bibo: <http://purl.org/ontology/bibo/>
	PREFIX dc: <http://purl.org/dc/terms/>
	PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
	PREFIX faldo: <http://biohackathon.org/resource/faldo#>
"""

# Lista con los nombres de las variables que obtenemos de la base de datos.
paramList = ['protein', 'proteinFullName', 'geneName', 'organismName', 'diseaseAnnotationText', 'domainFullName', 'similarityAnnotationText', 'locationAnnotationText', 'functionAnnotationText', 'pharmaceuticalAnnotationText'];

def buildQuery(proteinId, proteinName, geneName, organismName, diseaseAnnotation, domainName, similarityAnnotation, locationAnnotation, functionAnnotation, pharmaceuticalAnnotation):
	query = COMMON_PREFIXES
	query += "select distinct \n"
	query += "  ?protein\n"
	query += "	?proteinFullName\n"
	query += "	?geneName\n"
	query += "	?organismName\n"
	query += "	?diseaseAnnotationText\n"
	query += "	?domainFullName\n"
	query += "	?similarityAnnotationText\n"
	query += "	?locationAnnotationText\n"
	query += "	?functionAnnotationText\n"
	query += "	?pharmaceuticalAnnotationText\n"
	query += "where{\n"

	query += "	?protein a up:Protein .\n"

	if (proteinId != ''):
		query += "	VALUES ?protein {uniprotkb:"+ proteinId + "}\n"

	query += "\n"

	if (proteinName == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:recommendedName ?proteinName .\n"
	query += "	?proteinName up:fullName ?proteinFullName . \n"
	if (proteinName !=''):
		query += "	filter( regex(str(?proteinFullName), " + '"' + proteinName + '"' + ",\"i\" )) .\n"
	if (proteinName == ''):
		query += "	}\n"
	query += "\n"

	if (geneName == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:encodedBy ?gene .\n"
	query += "	?gene skos:prefLabel ?geneName .\n"
	if (geneName != ''):
		query += "	filter( regex(str(?geneName), " + '"' + geneName + '"' + ",\"i\" )) .\n"
	if (geneName == ''):
		query += "	}\n"

	query += "\n"

	if (organismName == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:organism ?organism .\n"
	query += "	?organism up:scientificName ?organismName .\n"
	if (organismName != ''):
		query += "	filter( regex(str(?organismName), " + '"' + organismName + '"' + ",\"i\" )) .\n"
	if (organismName == ''):
		query += "	}\n"

	query += "\n"

	if (diseaseAnnotation == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:annotation ?diseaseAnnotation .\n"
	query += "	?diseaseAnnotation a up:Disease_Annotation .\n"
	query += "	?diseaseAnnotation up:disease ?disease .\n"
	query += "	?disease rdfs:comment ?diseaseAnnotationText\n"
	if (diseaseAnnotation != ''):
		query += "	filter( regex(str(?diseaseAnnotationText), " + '"' + diseaseAnnotation + '"' + ",\"i\" )) .\n"
	if (diseaseAnnotation == ''):
		query += "	}\n"

	query += "\n"

	if (domainName == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:domain ?domain .\n"
	query += "	?domain up:recommendedName ?domainName .\n"
	query += "	?domainName up:fullName ?domainFullName .\n"
	if (domainName != ''):
		query += "	filter( regex(str(?domainFullName), " + '"' + domainName + '"' + ",\"i\" )) .\n"
	if (domainName == ''):
		query += "	}\n"

	query += "\n"

	if (similarityAnnotation == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:annotation ?similarityAnnotation .\n"
	query += "	?similarityAnnotation a up:Similarity_Annotation .\n"
	query += "	?similarityAnnotation rdfs:comment ?similarityAnnotationText .\n"
	if (similarityAnnotation != ''):
		query += "	filter( regex(str(?similarityAnotationText), " + '"' + similarityAnnotation + '"' + ",\"i\" )) .\n"
	if (similarityAnnotation == ''):
		query += "	}\n"

	query += "\n"

	if (locationAnnotation == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:annotation ?locationAnnotation .\n"
	query += "	?locationAnnotation a up:Subcellular_Location_Annotation .\n"
	query += "	?locationAnnotation up:locatedIn ?location .\n"
	query += "	?location up:cellularComponent ?cellComponent .\n"
	query += "	?cellComponent rdfs:comment ?locationAnnotationText .\n"
	if (locationAnnotation != ''):
		query += "	filter( regex(str(?locationAnnotationText), " + '"' + locationAnnotation + '"' + ",\"i\" )) .\n"
	if (locationAnnotation == ''):
		query += "	}\n"

	query += "\n"

	if (functionAnnotation == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:annotation ?functionAnnotation .\n"
	query += "	?functionAnnotation a up:Function_Annotation .\n"
	query += "	?functionAnnotation rdfs:comment ?functionAnnotationText .\n"
	if (functionAnnotation != ''):
		query += "	filter( regex(str(?functionAnnotationText), " + '"' + functionAnnotation + '"' + ",\"i\" )) .\n"
	if(functionAnnotation == ''):
		query += "	}\n"

	query += "\n"

	if (pharmaceuticalAnnotation == ''):
		query += "	OPTIONAL {\n"
	query += "	?protein up:annotation ?pharmaceuticalAnnotation .\n"
	query += "	?pharmaceuticalAnnotation a up:Pharmaceutical_Annotation .\n"
	query += "	?pharmaceuticalAnnotation rdfs:comment ?pharmaceuticalAnnotationText .\n"
	if (pharmaceuticalAnnotation != ''):
		query += "	filter( regex(str(?pharmaceuticalAnnotationText), " + '"' + pharmaceuticalAnnotation + '"' + ",\"i\" )) .\n"
	if (pharmaceuticalAnnotation == ''):
		query += "	}\n"
	query += "}\n"
	query += "limit 10\n"

	return query

def printResults(entrada, output):
	# Abrir fichero para escritura
	fileRes = open(output, 'w')
	resultados = entrada["results"]["bindings"]
	dump = json.dumps(resultados)
	fileRes.write(dump)

	"""# Imprimir cabecera
	for param in paramList:
		fileRes.write(param + "\t")
	fileRes.write("\n")

	# El formato json se puede recorrer de esta manera
	# para ir obteniendo valores de la respuesta.
	for entrada in json["results"]["bindings"]:
		for param in paramList:
			if (entrada.get(param)):
				fileRes.write(entrada.get(param)["value"] + "\t")
			else:
				fileRes.write("\t")
		fileRes.write("\n")"""
	fileRes.close()
	print(dump)

def printResultados(entrada):
	resu = entrada["results"]["bindings"]
	list_ = []
	for elemen in resu: 
		claves = elemen.keys()
		for item in claves: 
			elemen[item] = elemen.get(item)["value"]
		list_.append(elemen)
	return list_

@lru_cache
def sparqlwrap(proteinId, proteinName, geneName, organismName, diseaseAnnotation, domainName, similarityAnnotation, locationAnnotation, functionAnnotation, pharmaceuticalAnnotation):
	result = {}
	query = buildQuery(proteinId, proteinName, geneName, organismName, diseaseAnnotation, domainName, similarityAnnotation, locationAnnotation, functionAnnotation, pharmaceuticalAnnotation)
	print (query)

	# Creamos un objeto del tipo SPARQLWrapper indicando en que
	# direccion esta el servicio que recibe consultas en sparql
	# y responde a estas.
	sparql = SPARQLWrapper('https://sparql.uniprot.org/sparql')

	# Especificamos la consulta que queremos hacer en sparql.
	sparql.setQuery(query)

	# Indicamos en que formato queremos que nos devuelva
	# los resultados de la consulta. Puede ser json, xml,
	# rfd, turtle... Simplemente son distintos formatos
	# para representar los datos en ficheros de texto.
	sparql.setReturnFormat(JSON)

	# Esta es la instruccion que realiza la consulta a
	# uniprot. Devuelve un objeto de python que hay que
	# tratar.
	print ("Ejecutando query")
	results = sparql.query()

	# Con esto, convertimos el objeto devuelto por
	# el servicio al formato que especificamos antes.
	# En este caso, json.
	print ("Conviertiendo a json")
	json_ = results.convert()
	print ("Fin conversion a json")

	# Dentro de la variable results tenemos informacion
	# (metadatos) de lo que ha devuelto el servidor de
	# uniprot.
	# print (results.info())
	
	# Imprimir resultados
	return printResultados(json_)
