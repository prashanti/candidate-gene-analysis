def main():
	genes=open('../data/zfin-genes.txt','r')
	genelist=[]
	startline=1
	linecount=0
	for line in genes:
		linecount+=1
		data=line.split("\t")
		if linecount>=startline:
			genelist.append(data[0].replace("ZFIN:","http://zfin.org/").strip())
	nodata=[]
	for geneid in genelist:
		url="http://pkb-new.nescent.org/kb/gene/eq?id="+geneid.strip()
		content = urllib2.urlopen(url).read()
		raw = nltk.clean_html(content)
		data = json.loads(raw)
		if data:
			for line in data:
				if line['entity']:
					for entity in line['entity']:
						if line['quality']:
							for quality in line['quality']:
								temp=entity.split("http://purl.obolibrary.org/obo/")[1]
								quality=quality.split("http://purl.obolibrary.org/obo/")[1]
								print geneid+"\t"+temp+"\t"+quality
						else:
							print geneid+"\t"+entity+"\t"+""
				else:
					for quality in line['quality']:
						quality=quality.split("http://purl.obolibrary.org/obo/")[1]
						print geneid+"\t"+""+"\t"+quality
		else:
			nodata.append(geneid)				
	#print "nodata",nodata		
	

if __name__ == "__main__":
	import json
	import urllib2
	import nltk
	main()