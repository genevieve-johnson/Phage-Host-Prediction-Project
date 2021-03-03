import sys
import os
import csv
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
import numpy


##THIS SECTION WAS DONE BY GENEVIEVE JOHNSON UNTIL LINE BELOW
'''
# TEST SEQUENCES
# will code in getting viral sequence from input
virusSeq = 'ATCGATCGATC'
bactSeqs = ["ATGTCTGTGCGGATGCAGTCTCCCTGT", "AAATCGGGTGAAGTCCCTGTGTCAGACTG"]
x=kmer_comparison(virusSeq,bactSeqs)
'''

# if you were reading in from fasta files...
path = 'C:\\Users\\genev\\Desktop\\PUTONTI_LAB\\'
virus_file='Spike_seq.fasta'
bact_file='bact.fasta'

v=list(SeqIO.parse(virus_file,'fasta')) # list of viral sequences
b=list(SeqIO.parse(bact_file,'fasta')) # list of bacterial sequences


def crispr(qry_file):
    path='C:\\Users\\genev\\Documents\\PUTONTI_LAB\\'
    db_file='spacer_database'
    output_file='blast_trial_SRR5535767'
    command='blastn -query '+path+qry_file+' -db '+path+db_file+' -max_target_seqs 2 -perc_identity 100 -outfmt="10 qseqid sseqid qcovs pident length bitscore" -out '+path+output_file
    print(os.system(command))


    bacteria = ""
    with open(path+output_file, 'r') as f:
        reader = csv.reader(f)
        your_list = list(reader)
        results = {}
        for z in your_list:
            if z[0] in results.keys():
                results[z[0]].append(z[1])
            else:
                results[z[0]]=[]
                results[z[0]].append(z[1])    
    return results


#need list of all the bacteria that are exact matches, not just the first
#one since there will be several exact matches

#don't need to save the scores though so it's just stored in a list and then
#the list is printed

Entrez.email = 'cputonti@luc.edu'

def get_taxa_for_acc(accessions): #accessions is a list of accession numbers
    handle=Entrez.esummary(db="nucleotide",id=','.join(accessions),retmode="xml")
    records=Entrez.read(handle)
    handle.close()

    tax_id=[]
    for r in records:
        tax_id.append(str(r['TaxId']))
    handle=Entrez.efetch(db="Taxonomy",id=','.join(tax_id),retmode="xml")
    records=Entrez.read(handle)
    handle.close()

    taxa_names=[]
    for r in records:
        taxa_names.append(r['ScientificName'])
    return taxa_names  #returns a list of genus+species+strain in order of the accession numbers passed in

def get_taxa(header):
    acc=[]
    accessions=header.split("|")
    for i in accessions:
        a=i.split("_")
        if len(a)==4:   #it's a refseq NC_
            acc.append(a[0]+"_"+a[1])
        elif len(a)==3:
            acc.append(a[0])
        else:
            print("I do not know what to do with "+i)
    return get_taxa_for_acc(acc)


z = crispr(virus_file)
if z != "":
    for a in z.keys():
        for b in z[a]:
            x = get_taxa(b)
            if x == "":
                print("Sorry no CRISPR results for this sequence!")
            else:
                almost_there = set(x)
                final = str(almost_there).replace("{","").replace("}","").replace("'","")
        print("CRISPR RESULTS FOR " + a + ":")
        print("Your virus has matched to " + str(final) + " as its host bacteria.")

print("\n")


print("GENETIC HOMOLOGY RESULTS: ")

def homology(qry_file_2, user_threshold):
    fasta_query = open(qry_file_2).read()
    result_handle = NCBIWWW.qblast(program = "blastn", database = "nt", entrez_query = "Viruses[ORGN]", sequence = fasta_query)

    with open("my_blast2.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    blast = open("my_blast2.xml")

    blast_records= NCBIXML.parse(blast)
    blast_records = list(blast_records)

    Entrez.email = "gjohnson6@luc.edu"

    for i in blast_records:
        for alignment in i.alignments:
            for hsp in alignment.hsps:
                if hsp.bits > user_threshold:
                    titles = alignment.title
                    splits = titles.split("|")
                    handle = Entrez.efetch(db="nucleotide", id=splits[3], rettype="gb", retmode="text")
                    record = SeqIO.read(handle, "genbank")
                    print("Results for " + splits[3] + ":")
                    for feet in record.features:
                        if feet.type == 'source':
                            if 'host' in feet.qualifiers.keys():
                                print("possible host bacteria:", feet.qualifiers['host'][0])
                            if 'lab_host' in feet.qualifiers.keys():
                                print("lab host bacteria strain:", feet.qualifiers['lab_host'][0]) 
                    print("related phage:", splits[4])
                    print("bitscore:", hsp.bits)
                    print("\n")
                    done = ("End of Genetic Homology Results!")
    return done

w = homology(virus_file, 1000)
print(w)

#we are blasting the virus against all known viruses and then finding the matches' host bacteria
#need to save the bacteria into a data structure similar to a map but only the
#bacteria that reach above a certain threshold score and that score also
#needs to printed alongside the bacteria so that the user can decide


#--------------------------------------------------------------------------------------------------------#

##THIS SECTION BELOW WAS DONE BY TAYLOR MILLER-ENSMINGER


print("K-MER FREQUENCY RESULTS:")

# k-mer code for larger project
# import for rev comp
# calculate the % of the sequence seq for each 4^k k-mers
# returns dictionary of k-mer frequencies

def create_kmer_usage_profile(seq, k):

    mer_counts = {}

    my_rev_compl = Seq(str(seq))  # str(Seq(seq).reverse_complement())

    my_rev_compl = str(my_rev_compl.reverse_complement())

    for kk in range(0, len(seq) - k + 1):

        mer = seq[kk:kk + k]

        merC = my_rev_compl[kk:kk + k]

        if mer in mer_counts.keys():

            mer_counts[mer] += 1

        else:

            mer_counts[mer] = 1

        if merC in mer_counts.keys():

            mer_counts[merC] += 1

        else:

            mer_counts[merC] = 1

    totalMers = sum(mer_counts.values())

    for i in mer_counts.keys():

        mer_counts[i] /= float(totalMers)

    return mer_counts





# calculate the correlation of a bacteria and a virus

def correlationCalc(bact, virus):

    BactList = []

    VirusList = []

    x = set()

    x=set(list(bact.keys())+list(virus.keys()))

    for i in x:

        if i in bact.keys():

            BactList.append(bact[i])

        else:

            BactList.append(0)

        if i in virus.keys():

            VirusList.append(virus[i])

        else:

            VirusList.append(0)

    #print(numpy.corrcoef(BactList, VirusList)[0, 1])

    return numpy.corrcoef(BactList, VirusList)[0, 1]





# calculates the best correlation for a viral sequence and set of bacterial sequences

# input single viral sequence and list of bacterial sequences

def kmer_comparison(viral_seq, bacterial_seqs):

    # k-mer range considered 3 through 6



    # list of tuples for best correlations for each k

    # index of bacteria in bacterial_seqs and correlation value

    bestCorrK = []



    for k in range(3, 6):

        vKmers = create_kmer_usage_profile(viral_seq, k)

        corrValues = []

        for i in bacterial_seqs:

            bKmers = create_kmer_usage_profile(i, k)

            corrValues.append(correlationCalc(vKmers,bKmers))

        bestCorrK.append((corrValues.index(max(corrValues)),max(corrValues)))

    return bestCorrK

bactSeqs=[]
for i in b:
    bactSeqs.append(str(i.seq))

for i in v:
    virusSeq=str(i.seq)
    x=kmer_comparison(virusSeq,bactSeqs)
    output_string=i.id

    for j in x:
        Entrez.email = "cputonti@luc.edu"
        handle = Entrez.efetch(db="nucleotide", id=b[j[0]].id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        for feet in record.features:
            if feet.type == 'source':
                if 'organism' in feet.qualifiers.keys():
                    print("possible host bacteria:", feet.qualifiers['organism'][0])
                    print("kmer frequency: " + str(j[1]))
                    print('\n')

