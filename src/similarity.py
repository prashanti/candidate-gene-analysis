def populate_query_profiles():
	datafile=open("../data/Swartz_phenotypic_profiles.txt")
	queryannotations=dict()
	#Character	Character Label	State Symbol	State Label	Entity ID	Entity Label	Quality ID	Quality Label	Related Entity ID	Related Entity Label	
	for line in datafile:
		if "Character" not in line:
			characternum=int(line.split("\t")[0].strip())
			RE='null'
			E=line.split("\t")[4].strip()
			Q=line.split("\t")[6].strip()
			RE=line.split("\t")[8].strip()
			if characternum not in queryannotations:
				queryannotations[characternum]=[]
			queryannotations[characternum].append([E,Q,RE])
	datafile.close()
	return queryannotations

def populate_gene_profiles():
	datafile=open("../data/zfin-gene-annotations.txt")
	geneannotations=dict()
	for line in datafile:
		gene=line.split("\t")[0].replace("http://zfin.org/","").strip()
		E=line.split("\t")[1].strip()
		Q=line.split("\t")[2].strip()
		if gene not in geneannotations:
			geneannotations[gene]=[]
		geneannotations[gene].append([E,Q])
	return geneannotations



def get_name(term):
	term=term.replace(" ","")
	term=term.replace("(","")
	term=term.replace(")","")
	term=term.replace("and","")
	term=term.replace("some","")
	term=term.replace("\t","")
	return(term.strip())

def prepare_corpus(geneprofiles,superclasses,cc):
	outfile=open("../data/Annotation-Corpus.txt",'w')
	annotationcount=dict()
	corpussize=0
	termic=dict()
	annotationgenecount=dict()
	for gene in geneprofiles:
		tempancestorset=set()
		for annotation in geneprofiles[gene]:
			E=annotation[0]
			Q=annotation[1]
			classname=get_name(E+Q)

			if classname not in superclasses:
				superclasses=get_super_classes(E,Q,'',superclasses,cc)
			superclasslist=superclasses[classname]
			for superclass in superclasslist:
				tempancestorset.add(superclass)
		outfile.write(gene+"\t"+','.join(tempancestorset)+"\n")
		for anc in tempancestorset:
			if anc not in annotationgenecount:
				annotationgenecount[anc]=0
			annotationgenecount[anc]+=1
	return superclasses,annotationgenecount,len(geneprofiles)
	
def computeIC(annotationgenecount,corpussize):
	termic=dict()
	for term in annotationgenecount:
		ic=-math.log(float(annotationgenecount[term])/float(corpussize))
		termic[term]=ic
	return termic


def get_super_classes(E,Q,RE,superclasses,cc):

	key=get_name(E+Q+RE)
	if E and E not in superclasses:
		superclasses=load_ancestors(E,superclasses,cc)
	if Q and Q not in superclasses:
		
		superclasses=load_ancestors(Q,superclasses,cc)
	if RE and RE not in superclasses:
		superclasses=load_ancestors(RE,superclasses,cc)

	if key not in superclasses:
		superclasses[key]=set()

	if Q and E and RE:
		for ancE in superclasses[E]:
			for ancQ in superclasses[Q]:
				for ancRE in superclasses[RE]:
						superclasses[key].add(get_name(ancE+ancQ+ancRE).strip())
						superclasses[key].add(get_name(ancE+ancQ).strip())
						
	elif E and Q:
		for ancE in superclasses[E]:
			for ancQ in superclasses[Q]:
				superclasses[key].add(get_name(ancE+ancQ).strip())
				
	elif Q and RE:
		for ancRE in superclasses[RE]:
			for ancQ in superclasses[Q]:
				superclasses[key].add(get_name(e2+q1).strip())

	elif E and RE:
		for ancE in superclasses[E]:
			for ancRE in superclasses[RE]:
					superclasses[key].add(get_name(ancE+ancRE).strip())
					superclasses[key].add(get_name(ancE).strip())					
	elif E:
		for ancE in superclasses[E]:
			superclasses[key].add(get_name(ancE).strip())
	elif RE:
		for ancRE in superclasses[RE]:
			superclasses[key].add(get_name(ancRE).strip())
	elif Q:
		for ancQ in superclasses[Q]:
			superclasses[key].add(get_name(ancQ).strip())
	else:
		1
	return superclasses

def getsimilarity(EQ1,EQ2,):
	1



def load_ancestors(term,superclasses,cc):
	term=term.strip()
	superclasses[term]=set()
	#NOTE: The term itself needs to be added to the list of ancestors
	query="select ancestor from tbl_allancestorscurationexperiment where term = "+ "\""+term+"\""
	cc.execute(query)	
	data = cc.fetchall()
	superclasses[term].add(term)
	for row in data:
		superclasses[term].add(row[0])
	return superclasses





def main():
	db = MySQLdb.connect("localhost","root","somya-103","ontologies")
	cursor = db.cursor()
	queryprofiles=populate_query_profiles()
	geneprofiles=populate_gene_profiles()
	superclasses=dict()
	# superclasses=get_super_classes('','PATO_0000911','',superclasses,cursor)
	
	superclasses,annotationgenecount,corpussize=prepare_corpus(geneprofiles,superclasses,cursor)
	termic=computeIC(annotationgenecount,corpussize)
	print termic
	# if no RE then send ''








if __name__ == "__main__":
	import math
	import MySQLdb
	import sys
	main()