from __future__ import division
def load_UBERONPATO_superclasses(granularity):
	if granularity=="EA":
		infile=open("../data/EA_Superclasses.txt")
	else:
		infile=open("../data/EQ_Superclasses.txt")
	superclasses=dict()
	for line in infile:
		term,ancestor=line.split("\t")
		if term not in superclasses:
			superclasses[term]=set()
			superclasses[term].add(term)
		superclasses[term].add(ancestor.strip())
	return superclasses


def populate_EA_subsumers():
	superclasses=load_UBERONPATO_superclasses('EA')
	infile=open("../data/EA_Annotations_KBGenes2016.txt")
	outfile=open("../data/EASubsumers.txt",'w')
	annotations=set()
	easubsumersdict=dict()
	for line in infile:
		gene,E,A=line.split("\t")
		E=E.strip()
		A=A.strip()
		if E and E not in superclasses:
			superclasses[E]=set()
			superclasses[E].add(E)
		if A and A not in superclasses:
			superclasses[A]=set()
			superclasses[A].add(A)

		eaterm=E+A
		if eaterm not in easubsumersdict:
			if E and A:
				easubsumers = [''.join(x) for x in list(itertools.product(superclasses[E],superclasses[A]))]
				easubsumersdict[eaterm]=list(set.union(set(easubsumers),superclasses[A]))
			elif E and A=="":
				easubsumersdict[eaterm]=list(superclasses[E])

			elif E=="" and A:
				easubsumersdict[eaterm]=list(superclasses[A])

	infile.close()
	
	datafile=open("../data/EA_Swartz_phenotypic_profiles.txt")
	
	datafile.next()
	for line in datafile:
		line=line.replace(":","_")
		characternum=int(line.split("\t")[0].strip())
		E=line.split("\t")[4].strip()
		A=line.split("\t")[6].strip()
		RE=line.split("\t")[8].strip()
		key=E+A+RE
		if key not in superclasses:
			easubsumersdict[key]=set()
			easubsumersdict[key].add(key)
			
		expression=get_expression(E,A,RE)
		if expression in superclasses:
			easubsumersdict[key]=set.union(easubsumersdict[key],superclasses[expression])

		if A and E:		
			easubsumersdict[key] = [''.join(x) for x in list(itertools.product(superclasses[E],superclasses[A]))]		
			
		elif E:
			easubsumersdict[key]=set.union(easubsumersdict[key],superclasses[E])
		elif A:
			easubsumersdict[key]=set.union(easubsumersdict[key],superclasses[A])

	datafile.close()
	json.dump(easubsumersdict,outfile)



def populate_EQ_subsumers():
	superclasses=load_UBERONPATO_superclasses('EQ')
	infile=open("../data/Annotations_KBGenes2016.txt")
	outfile=open("../data/EQSubsumers.txt",'w')
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

	infile.close()
	
	datafile=open("../data/Swartz_phenotypic_profiles.txt")
	
	datafile.next()
	for line in datafile:
		line=line.replace(":","_")
		characternum=int(line.split("\t")[0].strip())
		E=line.split("\t")[4].strip()
		Q=line.split("\t")[6].strip()
		RE=line.split("\t")[8].strip()
		key=E+Q+RE
		if key not in superclasses:
			eqsubsumersdict[key]=set()
			eqsubsumersdict[key].add(key)
		expression=get_expression(E,Q,RE)
		if expression in superclasses:
			eqsubsumersdict[key]=set.union(eqsubsumersdict[key],superclasses[expression])

		if Q and E:		
			eqsubsumersdict[key] = [''.join(x) for x in list(itertools.product(superclasses[E],superclasses[Q]))]		
			
		elif E:
			eqsubsumersdict[key]=set.union(eqsubsumersdict[key],superclasses[E])
		elif Q:
			eqsubsumersdict[key]=set.union(eqsubsumersdict[key],superclasses[Q])

	datafile.close()
	json.dump(eqsubsumersdict,outfile)
	
def populate_E_subsumers():
	
	superclasses=load_UBERONPATO_superclasses('E')
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




def populate_query_profiles(granularity):
	if granularity=="E":
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
	elif granularity =="EA":
		datafile=open("../data/EA_Swartz_phenotypic_profiles.txt")
		queryannotations=dict()
		for line in datafile:
			line=line.replace(":","_")
			if "Character" not in line:
				characternum=int(line.split("\t")[0].strip())
				E=line.split("\t")[4].strip()
				A=line.split("\t")[6].strip()
				RE=line.split("\t")[8].strip()
				if characternum not in queryannotations:
					queryannotations[characternum]=set()
				queryannotations[characternum].add(E+A+RE)
		datafile.close()
	elif granularity =="EQ":
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
					queryannotations[characternum]=set()
				queryannotations[characternum].add(E+Q+RE)
		datafile.close()
	return queryannotations


def populate_gene_profiles(granularity):
	if granularity=='EQ':
		datafile=open("../data/Annotations_KBGenes2016.txt")
		geneannotations=dict()
		for line in datafile:
			gene=line.split("\t")[0].strip()
			E=line.split("\t")[1].strip()
			Q=line.split("\t")[2].strip()
			if gene not in geneannotations:
				geneannotations[gene]=set()
			geneannotations[gene].add(E+Q)
		
	if granularity=='EA':
		datafile=open("../data/EA_Annotations_KBGenes2016.txt")
		geneannotations=dict()
		for line in datafile:
			gene=line.split("\t")[0].strip()
			E=line.split("\t")[1].strip()
			A=line.split("\t")[2].strip()
			if gene not in geneannotations:
				geneannotations[gene]=set()
			geneannotations[gene].add(E+A)
		
	elif granularity =="E":
		datafile=open("../data/Annotations_KBGenes2016.txt")
		geneannotations=dict()
		for line in datafile:
			gene=line.split("\t")[0].strip()
			E=line.split("\t")[1].strip()
			if gene not in geneannotations:
				geneannotations[gene]=set()
			if E!='':
				geneannotations[gene].add(E)
	return geneannotations


def compute_profile_ic(geneprofiles,queryprofiles,subsumers,granularity):
	icout=open("../data/"+granularity+"_ProfileIC.txt",'w')
	icdict=dict()
	corpusprofiles=dict()
	freq=dict()
	for profileid in geneprofiles:
		corpusprofiles[profileid]=set()
		for annotation in geneprofiles[profileid]:
			corpusprofiles[profileid].add(annotation)
			
			for subsumer in subsumers[annotation]:
				corpusprofiles[profileid].add(subsumer)

	for profileid in queryprofiles:
		corpusprofiles[profileid]=set()
		for annotation in queryprofiles[profileid]:
			corpusprofiles[profileid].add(annotation)
			for subsumer in subsumers[annotation]:
				corpusprofiles[profileid].add(subsumer)

	for profileid in corpusprofiles:
		for annotation in corpusprofiles[profileid]:
			if annotation not in freq:
				freq[annotation]=0
			freq[annotation]+=1
	corpussize=len(corpusprofiles)
	print "corpussize",corpussize
	maxic=round(-math.log(1/corpussize),2)
	for annotation in freq:
		ic=round((-math.log(freq[annotation]/corpussize))/maxic,2)
		icdict[annotation]=ic
	outfile=open("../data/"+granularity+"_ProfileIC.txt",'w')
	json.dump(icdict,outfile)




def compute_annotation_ic(geneprofiles,queryprofiles,subsumers,granularity):
	
	icout=open("../data/"+granularity+"_AnnotationIC.txt",'w')
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

def calculate_groupwise_jaccard(profile1,profile2,ancestors):
    
    ancestors1=set()
    ancestors2=set()
    for term in profile1:
        ancestors1=set.union(ancestors1,ancestors[term])
    for term in profile2:
        ancestors2=set.union(ancestors2,ancestors[term])
    common=set.intersection(ancestors1,ancestors2)
    
    if len(common) > 0:
        union=set.union(ancestors1,ancestors2)  
        simj=len(common)/len(union)

    else:
        simj=0  
    return simj

def calculate_bestpairs_symmetric_jaccard(profile1,profile2,ancestors,similaritydict):
    finalsim=0
    bestmatchsimj1=[]
    bestmatchsimj2=[]
    termmatchsimj=[]
    for term1 in profile1:
        termmatchsimj=[]
        for term2 in profile2:
            termtuple=tuple(sorted((term1,term2)))
            if termtuple in similaritydict:
                simj=similaritydict[termtuple]
                termmatchsimj.append(simj)
            else:
                simj=getsimj(term1,term2,ancestors)
                similaritydict[termtuple]=simj
            termmatchsimj.append(simj)
        bestmatchsimj1.append(max(termmatchsimj))
    
    
    termmatchsimj=[]
    for term1 in profile2:
        termmatchsimj=[]
        for term2 in profile1:
            termtuple=tuple(sorted((term1,term2)))
            simj=similaritydict[termtuple]
            termmatchsimj.append(simj)
        bestmatchsimj2.append(max(termmatchsimj))
    
    
    return similaritydict,np.mean((np.median(bestmatchsimj1),np.median(bestmatchsimj2)))

def calculate_simgic(profile1,profile2,ancestors,icdict):
    
    ancestors1=set()
    ancestors2=set()
    commonic=0
    unionic=0
    for term in profile1:
        ancestors1=set.union(ancestors1,ancestors[term])
    for term in profile2:
        ancestors2=set.union(ancestors2,ancestors[term])
    common=set.intersection(ancestors1,ancestors2)
    
    if len(common) > 0:
        union=set.union(ancestors1,ancestors2) 
        for term in common:
            commonic=commonic+icdict[term]
        for term in union:
            unionic=unionic+icdict[term] 
        simgic=commonic/unionic    
    else:
        simgic=0  
    return simgic
    


def getsimj(term1,term2,ancestors):
	ancestors[term1]=set(ancestors[term1])
	ancestors[term2]=set(ancestors[term2])
	if len(set.union(ancestors[term1],ancestors[term2])) >0:
		simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
	else:
		simj=0
	return simj



def compile_evidence():
	infile=open("../data/Candidate_Genes_Paula.tsv")
	evidencedict=dict()
	infile.next()
	for line in infile:
		data=line.strip().split("\t")
		evidence=int(data[7])
		species=data[3].lower()
		data[1]=data[1].lower().strip()
		data[0]=data[0].lower().strip()
		if species=="zebrafish" and data[1]!="":
			genename=data[1]
		else:
			genename=data[0]
		names=genename.split(",")
		for name in names:
			name=name.strip()
			if name!="":
				if name not in evidencedict:
					evidencedict[name]=dict()
				if species not in evidencedict[name]:
					evidencedict[name][species]=evidence
				elif evidence < evidencedict[name][species]:
					evidencedict[name][species]=evidence
	outfile=open("../data/Compiled_Candidate_Evidence.txt",'w')
	json.dump(evidencedict,outfile)

def load_id_2_name():
	id2name=dict()
	datafile=open("../data/KBGenes_2016.tsv")
	datafile.next()
	for line in datafile:
		data=line.split("\t")
		geneid=data[0].replace("<","").replace(">","")

		genename=data[1].replace("\"","").replace("^^<http://www.w3.org/2001/XMLSchema#string>","")
		id2name[geneid]=genename.strip().lower()
	datafile.close()
	return id2name	

def populate_all_subsumers():
	populate_E_subsumers()
	populate_EA_subsumers()
	populate_EQ_subsumers()


def main():
	
	granularity=sys.argv[1]
	metric=sys.argv[2]
	id2name=load_id_2_name()
	# compile_evidence()
	# sys.exit()

	
	evidencedict=json.load(open("../data/Compiled_Candidate_Evidence.txt"))
	#populate_all_subsumers()
	
	queryprofiles=populate_query_profiles(granularity)
	geneprofiles=populate_gene_profiles(granularity)
	subsumers=json.load(open("../data/"+granularity+"Subsumers.txt"))
	#compute_profile_ic(geneprofiles,queryprofiles,subsumers,granularity)
	#compute_annotation_ic(geneprofiles,queryprofiles,subsumers,granularity)
	if "pic" in metric:
		icdict=json.load(open("../data/"+granularity+"_ProfileIC.txt"))
	else:
		icdict=json.load(open("../data/"+granularity+"_AnnotationIC.txt"))

		
	outfilename="../results/"+granularity+"_"+metric+"_Results.tsv"

	outfile=open(outfilename,'w')
	outfile.write("Character Number\tGeneID\tGenename\tMedian nIC\tEvidence\n")
	bpsimilaritydict=dict()
	for characternum in queryprofiles:
		for gene in geneprofiles:
			profile1=queryprofiles[characternum]
			profile2=geneprofiles[gene]
			print "Character",characternum
			if "zfin" in gene:
				species="zebrafish"
			elif "MGI" in gene:
				species="mouse"
			elif "ncbi" in gene:
				species="human"
			elif "xenbase" in gene:
				species="xenopus"
			else:
				print "No species",gene
			if len(profile2)==0:
				mediansim=0
			else:
				if metric =="bp_sym_aic_resnik" or metric =="bp_sym_pic_resnik":
					bpsimilaritydict,mediansim=calculate_bestpairs_symmetric_resnik(profile1,profile2,icdict,subsumers,bpsimilaritydict)
				elif metric=="aic_simGIC" or metric=="pic_simGIC":
						# variable is mediansim but no median is calculated for groupwise
						mediansim=calculate_simgic(profile1,profile2,subsumers,icdict)
				elif metric=="bp_sym_jaccard":
					bpsimilaritydict,mediansim=calculate_bestpairs_symmetric_jaccard(profile1,profile2,subsumers,bpsimilaritydict)
				elif metric=="groupwise_jaccard":
					mediansim=calculate_groupwise_jaccard(profile1,profile2,subsumers)
			genename=id2name[gene].strip()
			
			if genename in evidencedict and species in evidencedict[genename]:
				evidence=evidencedict[genename][species]
			else:
				evidence="N/A"
			outfile.write(str(characternum)+"\t"+gene+"\t"+genename+"\t"+str(mediansim)+"\t"+str(evidence)+"\n")

	outfile.close()

if __name__ == "__main__":
	import math
	import os
	import json
	import copy
	import itertools
	import MySQLdb
	import sys
	import numpy as np
	main()