def main():
	genes=open('../data/KBGenes_2016.tsv','r')
	genelist=[]
	startline=0
	linecount=0
	error=[]
	genes.next()
	for line in genes:
		line=line.replace("<","").replace(">","")
		genelist.append(data[0].strip())
	nodata=[]
	print len(genelist)
	for geneid in genelist:
		try:
			url="http://kb.phenoscape.org/api/gene/eq?id="+geneid.strip()
			content = urllib2.urlopen(url).read()
			data = json.loads(content)
			if data:
				for line in data:
					if line['entity']:
						for entity in line['entity']:
							if line['quality']:
								for quality in line['quality']:
									temp=entity.split("http://purl.obolibrary.org/obo/")[1]
									quality=quality.split("http://purl.obolibrary.org/obo/")[1]
									outfile.write(geneid+"\t"+temp+"\t"+quality+"\n")
							else:
								outfile.write(geneid+"\t"+entity+"\t"+""+"\n")
					else:
						for quality in line['quality']:
							quality=quality.split("http://purl.obolibrary.org/obo/")[1]
							outfile.write(geneid+"\t"+""+"\t"+quality+"\n")
			else:
				nodata.append(geneid)	
		except:
			error.append(geneid)
			pass			
			
	print "No data", nodata
	print "Error", error
	outfile.close()
if __name__ == "__main__":
	import json
	import sys
	import urllib2
	import nltk
	main()