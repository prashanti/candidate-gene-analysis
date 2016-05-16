from __future__ import division
def load_UBERONPATO_superclasses():
	infile=open("../data/Subsumers.txt")
	superclasses=dict()
	for line in infile:
		term,ancestor=line.split("\t")
		if term not in superclasses:
			superclasses[term]=set()
			superclasses[term].add(term)
		superclasses[term].add(ancestor.strip())
	return superclasses

def populate_gene_EQ_subsumers():
	superclasses=load_UBERONPATO_superclasses()
	infile=open("../data/Annotations_KBGenes2016.txt")
	
	annotations=set()
	eqsubsumersdict=dict()
	for line in infile:
		gene,E,Q=line.split("\t")
		E=E.strip()
		Q=Q.strip()
		if E and E not in superclasses:
			superclasses[E]=set()
			superclasses[E].add(E)
		if Q and Q not in superclasses:
			superclasses[Q]=set()
			superclasses[Q].add(Q)

		eqterm=E+Q
		if eqterm not in eqsubsumersdict:
			if E and Q:
				eqsubsumers = [''.join(x) for x in list(itertools.product(superclasses[E],superclasses[Q]))]
				eqsubsumersdict[eqterm]=list(set.union(set(eqsubsumers),superclasses[Q]))
			elif E and Q=="":
				eqsubsumersdict[eqterm]=list(superclasses[E])

			elif E=="" and Q:
				eqsubsumersdict[eqterm]=list(superclasses[Q])

	
	populate_profile_EQ_subsumers(eqsubsumersdict,superclasses,1)
	
def populate_E_subsumers():
	
	superclasses=load_UBERONPATO_superclasses()
	infile=open("../data/Annotations_KBGenes2016.txt")
	annotations=set()
	esubsumersdict=dict()
	for line in infile:
		gene,E,Q=line.split("\t")
		E=E.strip()
		if E and E not in superclasses:
			superclasses[E]=set()
			superclasses[E].add(E)
	
	for E in superclasses:
		superclasses[E]=list(superclasses[E])
	outfile=open("../data/ESubsumers.txt",'w')
	
	json.dump(superclasses,outfile)


def populate_profile_EQ_subsumers(eqsubsumers,superclasses,quality):
	datafile=open("../data/Swartz_phenotypic_profiles.txt")
	outfile=open("../data/EQSubsumers.txt",'w')
	datafile.next()
	for line in datafile:
		line=line.replace(":","_")
		characternum=int(line.split("\t")[0].strip())
		E=line.split("\t")[4].strip()
		Q=line.split("\t")[6].strip()
		RE=line.split("\t")[8].strip()


		key=E+Q+RE


		if key not in superclasses:
			eqsubsumers[key]=set()
			
		if quality ==1:
			expression=get_expression(E,Q,RE)
			eqsubsumers[key]=set.union(eqsubsumers[key],superclasses[expression])

		if Q and E:		
			eqsubsumers[key] = [''.join(x) for x in list(itertools.product(superclasses[E],superclasses[Q]))]		
			
		elif E:
			eqsubsumers[key]=set.union(eqsubsumers[key],superclasses[E])
		elif Q:
			eqsubsumers[key]=set.union(eqsubsumers[key],superclasses[Q])

	datafile.close()
	json.dump(eqsubsumers,outfile)




def populate_query_profiles(quality):
	if quality == 0:
		datafile=open("../data/Swartz_phenotypic_profiles.txt")
		queryannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			if "Character" not in line:
				characternum=int(line.split("\t")[0].strip())
				E=line.split("\t")[4].strip()
				if characternum not in queryannotations:
					queryannotations[characternum]=[]
				queryannotations[characternum].append(E)
		datafile.close()
		return queryannotations

	if quality ==1:
		datafile=open("../data/Swartz_phenotypic_profiles.txt")
		queryannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			if "Character" not in line:
				characternum=int(line.split("\t")[0].strip())
				E=line.split("\t")[4].strip()
				Q=line.split("\t")[6].strip()
				RE=line.split("\t")[8].strip()
				if characternum not in queryannotations:
					queryannotations[characternum]=[]
				queryannotations[characternum].append(E+Q+RE)
		datafile.close()
		return queryannotations


def populate_gene_profiles(quality):
	if quality ==1:
		datafile=open("../data/Annotations_KBGenes2016.txt")
		geneannotations=dict()
		for line in datafile:
			gene=line.split("\t")[0].strip()
			E=line.split("\t")[1].strip()
			Q=line.split("\t")[2].strip()
			if gene not in geneannotations:
				geneannotations[gene]=[]
			geneannotations[gene].append(E+Q)
		return geneannotations


	if quality ==0:
		datafile=open("../data/Annotations_KBGenes2016.txt")
		geneannotations=dict()
		for line in datafile:
			gene=line.split("\t")[0].strip()
			E=line.split("\t")[1].strip()
			if gene not in geneannotations:
				geneannotations[gene]=[]
			if E!='':
				geneannotations[gene].append(E)
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

def compute_ic(geneprofiles,queryprofiles,subsumers,quality):
	
	if quality ==1:
		icout=open("../data/AnnotationIC.txt",'w')
		annotationlist=[]
		for geneid in geneprofiles:
			for annotation in geneprofiles[geneid]:
				annotationlist+=list(subsumers[annotation])
		for profileid in queryprofiles:
			for annotation in queryprofiles[profileid]:
				annotationlist+=list(subsumers[annotation])

		corpussize=len(annotationlist)
		maxic=-math.log(1/corpussize,2)
		freq=dict()
		for annotation in annotationlist:
			if annotation not in freq:
				freq[annotation]=0
			freq[annotation]+=1
		
		icdict=dict()
		for annotation in set(annotationlist):
			p=freq[annotation]/corpussize
			print freq[annotation]
			ic=-math.log(p,2)
			icdict[annotation]=round(ic/maxic,2)


		json.dump(icdict,icout)
		icout.close()
		
		

	elif quality ==0:
		icout=open("../data/E_AnnotationIC.txt",'w')
		annotationlist=[]
		for geneid in geneprofiles:
			for annotation in geneprofiles[geneid]:

				annotationlist+=list(subsumers[annotation])
		for profileid in queryprofiles:
			for annotation in queryprofiles[profileid]:
				annotationlist+=list(subsumers[annotation])

		corpussize=len(annotationlist)
		maxic=-math.log(1/corpussize,2)
		freq=dict()
		for annotation in annotationlist:
			if annotation not in freq:
				freq[annotation]=0
			freq[annotation]+=1
		
		icdict=dict()
		for annotation in set(annotationlist):
			p=freq[annotation]/corpussize
			ic=-math.log(p,2)
			icdict[annotation]=round(ic/maxic,2)


		json.dump(icdict,icout)
		icout.close()



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

def getmicaic(term1,term2,ancestors,icdict):
    micaic=0
    mica=""
    commonancestors=set.intersection(set(ancestors[term1]),set(ancestors[term2]))
    lcslist=[icdict[anc] for anc in commonancestors]
    
    
    if len(lcslist)>0:
        micaic=np.max(lcslist)     
        return micaic
    else:
        return 0

def calculate_bestpairs_asymmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist=[]
    termmatchic=[]
    matchdata=[]
    for term1 in profile1:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        bestmatchiclist.append(np.max(termmatchic))

    return bpsimilaritydict,np.median(bestmatchiclist)

def calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,ancestors,bpsimilaritydict):
    finalsim=0
    bestmatchiclist1=[]
    bestmatchiclist2=[]
    termmatchic=[]
    matchdata=[]
    for term1 in profile1:
        termmatchic=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in bpsimilaritydict:
                termmatchic.append(bpsimilaritydict[termtuple])
            
            else:
                micaic=getmicaic(term1,term2,ancestors,icdict)
                termmatchic.append(micaic)
                bpsimilaritydict[termtuple]=micaic
        
        bestmatchiclist1.append(np.max(termmatchic))
    

    termmatchic=[]
    matchdata=[]
    for term1 in profile2:
        termmatchic=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            termmatchic.append(bpsimilaritydict[termtuple])
        bestmatchiclist2.append(np.max(termmatchic))
            
    return bpsimilaritydict,np.mean((np.median(bestmatchiclist1),np.median(bestmatchiclist2)))


def get_evidence(filename):
	eqresults=open("../data/EQResults_Symm_Resnik_CandidateGenes.tsv")
	eqresults.next()
	evidencedict=dict()
	genesymboldict=dict()
	
	for line in eqresults:
		data=line.strip().split("\t")
		characternum,geneid,genesymbol,evidence=data[0],data[1],data[2],data[4]
		evidencedict[characternum+"\t"+geneid]=evidence
		genesymboldict[characternum+"\t"+geneid]=genesymbol
	eqresults.close()
	outfile=open("../data/EResults_Symm_Resnik_CandidateGenes.tsv",'w')
	outfile.write("Character Number	GeneID	Genename	E IC	Evidence\n")
	eresults=open(filename)
	eresults.next()
	for line in eresults:
		data=line.strip().split("\t")
		characternum,geneid,score=data[0],data[1],data[2]
		evidence=evidencedict[characternum+"\t"+geneid]
		genesymbol=genesymboldict[characternum+"\t"+geneid]
		outfile.write(characternum+"\t"+geneid+"\t"+genesymbol+"\t"+str(score)+"\t"+evidence+"\n")
	outfile.close()
	sys.exit()


def main():
	granularity=sys.argv[1]
	get_evidence("../data/EComparisonScores_Symm_Resnik.tsv")

	if granularity=="EQ":
		# run only once. Subsumer dict is written to "../data/2016/EQSubsumers.txt"
		#populate_gene_EQ_subsumers()
		eqsubsumers=json.load(open("../data/EQSubsumers.txt"))
		icdict=json.load(open("../data/AnnotationIC.txt"))
		queryprofiles=populate_query_profiles(quality)
		geneprofiles=populate_gene_profiles(quality)	
		#compute_ic(geneprofiles,queryprofiles,eqsubsumers,quality)

		
		outfile=open("../data/EQComparisonScores_Symm_Resnik.tsv",'w')
		
		outfile.write("Character Number\tGene\tMedian nIC\n")
		bpsimilaritydict=dict()
		for characternum in queryprofiles:
			for gene in geneprofiles:
				profile1=queryprofiles[characternum]
				profile2=geneprofiles[gene]
				print "Comparing character",characternum," and gene ",gene
				print len(profile1),len(profile2)
				if len(profile2)==0:
					maxsimj=0
					maxic=0
				else:
					#bpsimilaritydict,medianic=calculate_bestpairs_asymmetric_resnik(profile1,profile2,icdict,eqsubsumers,bpsimilaritydict)

					bpsimilaritydict,medianic=calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,eqsubsumers,bpsimilaritydict)
					
				outfile.write(str(characternum)+"\t"+gene+"\t"+str(medianic)+ "\n")
		outfile.close()

	elif quality==0:
		queryprofiles=populate_query_profiles(quality)
		geneprofiles=populate_gene_profiles(quality)	
		#populate_E_subsumers()
		subsumers=json.load(open("../data/ESubsumers.txt"))
		compute_ic(geneprofiles,queryprofiles,subsumers,quality)
		icdict=json.load(open("../data/E_AnnotationIC.txt"))

		outfile=open("../data/EComparisonScores_Symm_Resnik.tsv",'w')

		outfile.write("Character Number\tGene\tMedian nIC\n")
		bpsimilaritydict=dict()
		for characternum in queryprofiles:
			for gene in geneprofiles:
				profile1=queryprofiles[characternum]
				profile2=geneprofiles[gene]
				print "Comparing character",characternum," and gene ",gene
				if len(profile2)==0:
					maxsimj=0
					maxic=0
				else:
					bpsimilaritydict,medianic=calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,subsumers,bpsimilaritydict)
					
				outfile.write(str(characternum)+"\t"+gene+"\t"+str(medianic)+ "\n")
		outfile.close()

	







if __name__ == "__main__":
	import math
	import json
	import itertools
	import MySQLdb
	import sys
	import numpy as np
	main()