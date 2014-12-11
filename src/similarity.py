def populate_query_profiles(quality):


	if quality == 0:
		datafile=open("../data/Phenotypic_Profiles_Subset.txt")
		queryannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			if "Character" not in line:
				characternum=int(line.split("\t")[0].strip())
				RE=''
				E=line.split("\t")[4].strip()
				if characternum not in queryannotations:
					queryannotations[characternum]=[]
				queryannotations[characternum].append([E])
		datafile.close()
		return queryannotations

	if quality ==1:
		datafile=open("../data/Phenotypic_Profiles_Subset.txt")
		queryannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
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


def populate_gene_profiles(quality):
	if quality ==1:
		datafile=open("../data/zfin-gene-annotations.txt")
		geneannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			gene=line.split("\t")[0].replace("http://zfin.org/","").strip()
			E=line.split("\t")[1].strip()
			Q=line.split("\t")[2].strip()
			if gene not in geneannotations:
				geneannotations[gene]=[]
			geneannotations[gene].append([E,Q])
		return geneannotations

	if quality ==0:
		datafile=open("../data/zfin-gene-annotations.txt")
		geneannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			gene=line.split("\t")[0].replace("http://zfin.org/","").strip()
			E=line.split("\t")[1].strip()
			if gene not in geneannotations:
				geneannotations[gene]=[]
			if E != '':	
				geneannotations[gene].append([E])
		return geneannotations


def get_name(term):
	term=term.replace(" ","")
	term=term.replace("(","")
	term=term.replace(")","")
	term=term.replace("and","")
	term=term.replace("some","")
	term=term.replace("\t","")
	term=term.replace("inheres_in","RO_0000052")
	return(term.strip())

def prepare_corpus(geneprofiles,superclasses,cc,quality):
	if quality ==1:
		outfile=open("../data/Annotation-Corpus_EQ1.txt",'w')
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
					superclasses=get_super_classes(E,Q,'',superclasses,cc,quality)
				superclasslist=superclasses[classname]
				for superclass in superclasslist:
					tempancestorset.add(superclass)
			outfile.write(gene+"\t"+','.join(tempancestorset)+"\n")
			for anc in tempancestorset:
				if anc not in annotationgenecount:
					annotationgenecount[anc]=0
				annotationgenecount[anc]+=1
		return superclasses,annotationgenecount,len(geneprofiles)
		

	if quality ==0:
		outfile=open("../data/Annotation-Corpus_E.txt",'w')
		annotationcount=dict()
		corpussize=0
		termic=dict()
		annotationgenecount=dict()
		for gene in geneprofiles:
			tempancestorset=set()
			for annotation in geneprofiles[gene]:
				E=annotation[0]
				classname=get_name(E)

				if classname not in superclasses:
					superclasses=get_super_classes(E,'','',superclasses,cc,quality)
				superclasslist=superclasses[classname]
					
				for superclass in superclasslist:
					tempancestorset.add(superclass)
			outfile.write(gene+"\t"+','.join(tempancestorset)+"\n")
			for anc in tempancestorset:
				if anc not in annotationgenecount:
					annotationgenecount[anc]=0
				annotationgenecount[anc]+=1
		return superclasses,annotationgenecount,len(geneprofiles)

def compute_IC(annotationgenecount,corpussize):
	termic=dict()
	maximum=round(-math.log(float(1)/float(corpussize)),2)
	for term in annotationgenecount:
		ic=-math.log(float(annotationgenecount[term])/float(corpussize))
		termic[term]=round(ic/maximum,2)
	return termic

def get_expression(E,Q,RE):
	if "Entity ID" in E:
		return ("null")

	expression="null"
	if len(Q)!=0:
		expression = Q
	
	if len(E)!=0:
		if expression != "null":
			expression=expression+" and inheres_in some ("+E+")"
		
		else:
			expression="inheres_in some ("+E+")"
		
	
	if len(RE)!=0:
		if expression != "null":
			expression=expression+" and towards some ("+RE+")"
		
		else:
			expression="towards some ("+RE+")"

	return(expression.strip())

def load_query_superclasses(queryprofiles,superclasses,cc,quality):
	if quality ==1:
		for characternum in queryprofiles:
			for annotation in queryprofiles[characternum]:
				superclasses=get_super_classes(annotation[0],annotation[1],annotation[2],superclasses,cc,quality)
		return superclasses

	if quality ==0:
		for characternum in queryprofiles:
			for annotation in queryprofiles[characternum]:
				superclasses=get_super_classes(annotation[0],'','',superclasses,cc,quality)
		return superclasses		

def get_super_classes(E,Q,RE,superclasses,cc,quality):

	key=get_name(E+Q+RE)


	if E and E not in superclasses:
		superclasses=load_ancestors(E,superclasses,cc)
		
	if Q and Q not in superclasses:
		superclasses=load_ancestors(Q,superclasses,cc)
	if RE and RE not in superclasses:
		superclasses=load_ancestors(RE,superclasses,cc)
	if key not in superclasses:
		superclasses[key]=set()
		
	if quality ==1:
		expression=get_expression(E,Q,RE)
		if expression and expression not in superclasses:
			superclasses=load_ancestors(expression,superclasses,cc)	
		for ancExpression in superclasses[expression]:
			# these ancestors are already in get_name format.
			superclasses[key].add(ancExpression)	


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


def get_EQ_similarity(EQ1,EQ2,superclasses,termic):
	annotation1="".join(EQ1)
	annotation2="".join(EQ2)
	class1=get_name(annotation1)
	class2=get_name(annotation2)
	commonancestors=set.intersection(superclasses[class1],superclasses[class2])
	maxic=0
	simj=0
	simj=float(len(set.intersection(superclasses[class1],superclasses[class2])))/float(len(set.union(superclasses[class1],superclasses[class2])))
	# print "annotations being compared",annotation1,annotation2
	# print "common ancestor set",set.intersection(superclasses[class1],superclasses[class2])
	# print "unionset",set.union(superclasses[class1],superclasses[class2])
	# print "\n\n\n"
	for anc in commonancestors:
		ic=termic[anc]
	
		if ic > maxic:
			maxic=ic

	return simj,maxic


def get_E_similarity(EQ1,EQ2,superclasses,termic):

	annotation1=EQ1[0]
	annotation2=EQ2[0]
	class1=get_name(annotation1)
	class2=get_name(annotation2)
	commonancestors=set.intersection(superclasses[class1],superclasses[class2])
	maxic=0
	simj=0
	simj=round(float(len(set.intersection(superclasses[class1],superclasses[class2])))/float(len(set.union(superclasses[class1],superclasses[class2]))),2)
	for anc in commonancestors:
		ic=termic[anc]
	
		if ic > maxic:
			maxic=ic

	return simj,maxic

def load_ancestors(term,superclasses,cc):
	term=term.strip()
	superclasses[term]=set()
	#NOTE: The term itself needs to be added to the list of ancestors
	query="select ancestor from tbl_candidategeneanalysis where term = "+ "\""+term+"\""
	cc.execute(query)	
	data = cc.fetchall()
	superclasses[term].add(get_name(term))
	for row in data:
		superclasses[term].add(row[0])
		
	return superclasses





def main():
	quality=int(sys.argv[1])
	cursor = db.cursor()
	superclasses=dict()


	queryprofiles=populate_query_profiles(quality)
	geneprofiles=populate_gene_profiles(quality)	
	superclasses,annotationgenecount,corpussize=prepare_corpus(geneprofiles,superclasses,cursor,quality)
	superclasses=load_query_superclasses(queryprofiles,superclasses,cursor,quality)



	termic=compute_IC(annotationgenecount,corpussize)
	if quality ==1:
		outfile=open("../data/EQComparisonScores1.tsv",'w')
	if quality ==0:
		outfile=open("../data/EComparisonScores.tsv",'w')

	outfile.write("Character Number\tGene\tMax SimJ\t Max nIC\n")


	# profile1=queryprofiles[1]
	# profile2=geneprofiles['http_//zfin.org/ZDB-GENE-070117-1345']


	# for annotation1 in profile1:
	# 	for annotation2 in profile2:
	# 		simj,ic=get_EQ_similarity(annotation1,annotation2,superclasses,termic)
	# 		print simj,ic

	
	for characternum in queryprofiles:
		for gene in geneprofiles:
			iclist=[]
			simjlist=[]
			profile1=queryprofiles[characternum]
			profile2=geneprofiles[gene]
			
			print "Comparing character",characternum," and gene ",gene

			if len(profile2)==0:
				maxsimj=0
				maxic=0
			else:
				for annotation1 in profile1:
					for annotation2 in profile2:
						if quality ==1:

							simj,ic=get_EQ_similarity(annotation1,annotation2,superclasses,termic)
						if quality==0:
							simj,ic=get_E_similarity(annotation1,annotation2,superclasses,termic)
						simjlist.append(simj)
						iclist.append(ic)
				maxsimj=numpy.max(simjlist)
				maxic=numpy.max(iclist)
			outfile.write(str(characternum)+"\t"+gene+"\t"+str(maxsimj)+"\t"+str(maxic)+ "\n")
	outfile.close()







if __name__ == "__main__":
	import math
	import MySQLdb
	import sys
	import numpy
	main()