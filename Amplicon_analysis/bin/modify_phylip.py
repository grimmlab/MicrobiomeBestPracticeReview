from Bio import SeqIO        
import re            
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="fasta file to be modified")
args = parser.parse_args()
print(args.fasta_file)
fasta_file = args.fasta_file

new_names = []

fasta = SeqIO.parse(handle=fasta_file, format='fasta')

for seq in fasta:
    # print(seq.description)
    otu = re.search('(Otu[0-9]+)\|[0-9]+', seq.description)
    # print(otu.group(1))  
    seq.id = otu.group(1)
    new_names.append(seq)

newfasta_name = str(fasta_file)
nu = newfasta_name.find('fasta')
_new_fastafile = newfasta_name[0:nu] + 'otu_modified.fasta'
#print (self._new_fastafile)
SeqIO.write(sequences=new_names, format='fasta', handle=_new_fastafile)  
