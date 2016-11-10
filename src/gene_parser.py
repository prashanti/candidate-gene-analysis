import requests
import sys

"""
gene_parser.py:

    Parses genes from a source file (ZDB, MGI, NIH, or XB type) and finds a corresponding identifier in a target file.
    Genes and identifiers are printed to an output file.
    If a specific gene is not found in the target file, program will print "Not found" in the output file.
    If a specific gene is not found in the target file, program will ALSO check rest.ensembl.org in case a particular
        target file is missing entries that can be found on ensembl.
"""

__author__ = "Alex Hahn"
__version__ = "1.1.2"
__email__ = "ashahn@uncg.edu"

"""
    Change Log

        1.0
            - Program parses genes from a source file (KBGenes) and gets ensembl IDs from 4 text files (one for each
                species type).
            - If a gene is not found, then output "Not found" in the output file.
            - Output is tab-separated with a URL id followed by its ensembl ID.

        1.1
            - If a gene is not found in a text file, then search rest.ensembl.org for its ID. If it is found on
                ensembl.org, then output its correct ID. If it is still not found, output "Not found".
            - rest.ensembl.org is queried using the 'urllib2' and 'json' packages.

        1.1.1
            - Version 1.1 had stability issues with urllib2 and urlopen(). Some executions ran smoothly, while others
                would timeout and fail.
            - rest.ensembl.org requests are now performed using the 'requests' package.
                - So far, this has resolved the timeout issue.

        1.1.2
            - Changes how the program outputs results
            - For a given input file <source> stored at ../data/<source>, the output is now stored at
                ../results/<source>_EnsemblIDs.tsv
            - Each output file is now a copy of the input file with an additional column at the end of each line that
                with the ensembl IDs.
"""


class Parser:

    def __init__(self):
        # Ensembl server used for json requests
        self.SERVER = "https://rest.ensembl.org"

        # Location of search files for each gene type
        self.ZDB_DIR = "../data/ensembl_1_to_1.txt"
        self.MGI_DIR = "../data/MGI_Gene_Model_Coord.rpt"
        self.NIH_DIR = "../data/gene2ensembl.tsv"
        self.XB_DIR = "../data/GenePageEnsemblModelMapping.txt"

    # This method reads each line of the source file and determines which type of gene the current line contains. Then,
    # it calls the search_file method to get the results
    def read_source(self, source_file):
        source = open(source_file, "r")

        target_file = "../data/" + source_file[source_file.rfind("/") + 1:]
        target_file = target_file[:target_file.rfind(".")] + "_EnsemblIDs.tsv"
        target = open(target_file, "w+")

        # Read each line in the input file
        # Search its corresponding gene file
        # Output the results to the output file
        for line in source:
            short_gene = None
            ext = None
            result = None

            if line.find("zfin.org/ZDB-GENE") != -1 or line.find("zebrafish") != -1:
                if line.find("zebrafish") != -1:
                    short_gene = line[:line.find('\t')]
                    search_id = short_gene
                else:
                    short_gene = line[line.find("\"") + 1:line.find("\"", line.find("\"") + 1)]
                    search_id = line[line.find("ZDB"):line.find(">", line.find("ZDB"))]

                ext = "/xrefs/symbol/danio_rerio/" + short_gene + "?"
                result = self.search_file(self.ZDB_DIR, search_id, 3)

            elif line.find("www.informatics.jax.org/marker/MGI") != -1 or line.find("mouse") != -1:
                if line.find("mouse") != -1:
                    short_gene = line[:line.find('\t')]
                    search_id = short_gene
                else:
                    short_gene = line[line.find("\"") + 1:line.find("\"", line.find("\"") + 1)]
                    search_id = line[line.find("MGI"):line.find(">", line.find("MGI"))]

                ext = "/xrefs/symbol/mus_musculus/" + short_gene + "?"
                result = self.search_file(self.MGI_DIR, search_id, 10)

            elif line.find("www.ncbi.nlm.nih.gov") != -1 or line.find("human") != -1:
                if line.find("human") != -1:
                    short_gene = line[:line.find('\t')]
                    search_id = short_gene
                else:
                    short_gene = line[line.find("\"") + 1:line.find("\"", line.find("\"") + 1)]
                    search_id = short_gene

                ext = "/xrefs/symbol/homo_sapiens/" + short_gene + "?"
                result = self.search_file(self.NIH_DIR, search_id, 2)

            # For XB-Gene, use the gene name in col 1 from the source file. Using the ID from col 0 did not yield any
            # results in the search file.
            elif line.find("xenbase.org/XB-GENEPAGE") != -1 or line.find("xenopus") != -1:
                if line.find("xenopus") != -1:
                    short_gene = line[:line.find('\t')]
                    search_id = short_gene
                else:
                    short_gene = line[line.find("\"") + 1:line.find("\"", line.find("\"") + 1)]
                    search_id = line[line.find("\"", line.find("XB-GENE")) + 1:line.find("\"", line.find("\"") + 1)]

                ext = "/xrefs/symbol/xenopus/" + short_gene + "?"
                result = self.search_file(self.XB_DIR, search_id, 3)

            # If result was not found in the target files, then check the ensembl database. If it is found on ensembl,
            # then set result to be the correct ID. If not found, then result is still "Not found".
            if result == "Not found" and ext is not None:
                r = requests.get(self.SERVER+ext, headers={"Content-Type": "application/json"})

                if not r.ok:
                    r.raise_for_status()
                    sys.exit()

                data = r.json()
                if data:
                    result = data[0]["id"]

            # Print results to the output file.
            if short_gene is not None and result is not None:
                target.write(line.rstrip() + "\t" + result)

                if result.rfind("\n") == -1:
                    target.write("\n")

                target.flush()

        source.close()
        target.close()

    # This method can search any of the gene search files. It takes the file path, search key, and a column number as
    # parameters. The column number specifies which file to input from the search file.
    def search_file(self, path, key, col):
        s_file = open(path, "r")

        for line in s_file:
            if line.find(key) != -1:
                return line.split('\t')[col]

        return "Not found"

# Main method
if __name__ == '__main__':

    # Make a Parser object and complete the gene parsing.
    parser = Parser()
    parser.read_source("../data/KBGenes_2016.tsv")
    parser.read_source("../data/preprocessed_candidatelist_09_17_16.txt")


