# Blasts a fasta file against a Reference Sequences Library, parses the results and creates a result Excel file

# FIRST STEP: create the local BLAST indices with this command line in Python terminal:
# makeblastdb -in ReferenceSequencesDatabase_16S_test.fasta -parse_seqids -out ref_library -dbtype nucl
# Then obtain the BLAST result in .txt format

# Packages
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from Bio import SeqIO
import xlsxwriter
import pandas as pd
from operator import itemgetter
import os.path

# USER VARIABLES
input_file = input("Insert file name (without file extension):")  # file has to be in the same folder of the script
ref_library = input("Insert Reference Sequences Library name (without file extension):")
cutoff_pident = input("Percent identity cut-off in BLAST:")  # does not save in .txt file results below a certain cut-off, e.g. <90%

# SECOND STEP: run the local BLAST
cline = NcbiblastnCommandline(query=input_file + ".fasta", db="ref_library",
                              evalue=0.001, max_target_seqs=150,
                              perc_identity=cutoff_pident, out=input_file + ".txt",
                              outfmt="6 std qlen slen gaps qcovs")
# for specific fields check blastn -help in outfmt, otherwise 12 standard fields
# STD = "qseqid" "sseqid" "pident" "length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore"
# ADDED = "qlen" "slen" "gaps" "qcovs"
print(cline)
cline()  # Execute blast command line


# THIRD STEP: parse the local blast result
best_blasts = []  # Create an empty blast result list
my_fields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps qcovs"
QueryResults = SearchIO.parse(input_file + ".txt", "blast-tab", fields=my_fields)

for QueryResult in QueryResults:
    to_sort = []  # Create an empty list to load the hits for each blast search
    for hit in QueryResult.hits:
        for hsp in hit.hsps:
            qcov = (hsp.query_end - hsp.query_start) / QueryResult.seq_len * 100  # no +1 (Python starts from 0)
            if (qcov >= 75) or (hsp.aln_span > 0.9 * hit.seq_len):
                # Check for aln length, in case of very short queries/subjects
                to_sort.append((QueryResult.id, hsp.ident_pct, hit.id, qcov, QueryResult.seq_len, hit.seq_len,
                                hsp.aln_span, hsp.evalue))
    # "query", "%id", "best_match", "%qcov", "qlen", "slen", "aln_len", "evalue"
    sorted_by_perc_id = sorted(to_sort, reverse=True, key=itemgetter(1))  # Sort hits by %id by descending order
    if len(sorted_by_perc_id) != 0:
       best_blasts.append(sorted_by_perc_id[0])  # Put the first result of the list in the blast result list

# Load the blast result list in a dataframe
df = pd.DataFrame(best_blasts, columns=("query", "%id", "best_match", "%qcov", "qlen", "slen", "aln_len", "evalue"))

print("Parsing completed!")


# FOURTH STEP: checking results. Additional info from Check_tags.csv file (Columns: "Label", "Check_tag")

intersp_div = input("Percent interspecific sequence divergence:")
path = "./" + "Check_tags.csv"
if os.path.isfile(path): # check if the file exist 
  check_tags = pd.read_csv("Check_tags.csv", sep=None, engine="python")  

df["notes"] = ""  # Add the empty notes column to the results dataframe
df_notes = []
qry_to_check = []
for qry, pid, sp, qc, note in zip(df["query"], df["%id"], df["best_match"], df["%qcov"], df["notes"]):
    if (pid != "NA" and pid < (100 - float(intersp_div))) or (qc < 75):
        note = "to_check"
        qry_to_check.append(qry)
    else:
        if os.path.isfile(path):
            for lab, tag in zip(check_tags["Label"], check_tags["Check_tag"]):
                if lab in sp:
                    note = tag
    df_notes.append(note)
df["notes"] = df_notes

# Saving blast results for only the results labelled to_check
df_to_check = pd.DataFrame(columns = ["query","%id","best_match","%qcov","qlen","slen","aln_len","evalue"])
all_blast = pd.read_csv(input_file + ".txt", sep="\t", names=my_fields.split(), index_col=None)
for check in qry_to_check:
   check_qry = all_blast[all_blast["qseqid"] == check]
   check_qry= check_qry.sort_values(by="pident", ascending=False)
   df_to_check = pd.concat([df_to_check, check_qry])

if not(df_to_check.empty):
  df_to_check = df_to_check[["qseqid","pident","sseqid","qcovs","qlen","slen","length","evalue"]]


#General statistics and counts

data_lineage = []  # Create an empty data statistics list
keys = []  # List of RSL taxa

# Create an empty dictionary with all the species lineages from the Reference Sequences Library + Outgroup
with open(ref_library + ".fasta", "r") as fasta:
    for k in SeqIO.parse(fasta, "fasta"):
        keys.append(k.id)
    keys.append("Results %id <" + cutoff_pident)
d = {key: [] for key in keys}

# Store all the codes in the input fasta to count the outgroups
samples = []
with open(input_file + ".fasta", "r") as fasta:
    for sam in SeqIO.parse(fasta, "fasta"):
        samples.append(sam.id)
outgroup = [sample for sample in samples if
            sample not in df["query"].values.tolist()]  # counts how many samples without blast result

# Search for the keys in the blast result list and append the sample identifier (qr) to that key
for qr, pi, tax, cov, qln, sln, aln, evl in best_blasts:
    # "query", "%id", "best_match", "%qcov", "qlen", "slen", "aln_len", "evalue"
    d.get(tax, []).append(qr)
d["Results %id <" + cutoff_pident] = outgroup
data_lineage = dict(sorted(d.items()))  # Sort data by taxon number in ascending order

n_exc = 0
n_missing = 0
for n in data_lineage.values():
    if not n:
        n.append("NA")  # Add NA if there are no samples for that lineage
        n_missing = n_missing + 1
    else:
        n.sort()
        n_exc = n_exc + len(n)  # Counts how many samples are available for that lineage
        n[0:] = [";".join(n[0:])]  # Merge all samples into one element with ";" as separator

df_lineage = pd.DataFrame.from_dict(data_lineage, orient="index", columns=["Sequences"])
df_lineage.reset_index(inplace=True)
df_lineage.rename(columns={"index": "Lineages"}, inplace=True)

# Counts how many sequences (lineages) are in the Reference Sequences Library
n_fasta = len([1 for line in open(input_file + ".fasta") if line.startswith(">")])

if data_lineage["Results %id <" + cutoff_pident][0] != "NA":
    outg = len(data_lineage["Results %id <" + cutoff_pident][0].split(";"))  # Counts of results below cut-off
else:
    outg = 0
n_taxa = len(data_lineage) - 1  # Counts of taxa. Remove 1 because of "Results %id<cut-off"
df_stats = pd.DataFrame({"Stats": ["Total taxa RSL", "Taxa without samples", "Taxonomic coverage (%)",
                                   "Total input sequences", "Sequences processed (%)",
                                   "N results %id <" + cutoff_pident],
                         "n": [n_taxa, n_missing, ((n_taxa - n_missing) / n_taxa*100), n_fasta, (n_exc / n_fasta)*100, outg]})


# FIFTH STEP: Excel with results creation

with pd.ExcelWriter(input_file + ".xlsx") as writer:
    df.to_excel(writer, index=False, sheet_name="Best_blast")
    df_to_check.to_excel(writer, index=False, sheet_name="Results_to_check")
    df_lineage.to_excel(writer, index=False, sheet_name="Lineages")
    df_stats.to_excel(writer, index=False, sheet_name="Summary")

    # Formatting
    workbook = writer.book
    format1 = workbook.add_format({"num_format": "0.0"})
    format2 = workbook.add_format({"num_format": "0.00E+00"})
    format3 = workbook.add_format({"bg_color": "red"})
    format4 = workbook.add_format({"align": "left"})

    worksheet1 = writer.sheets["Best_blast"]  # set red color for values under the interspecific divergence
    worksheet1.set_column("B:B", None, format1)  # %id and qcov in float format .1f
    worksheet1.set_column("D:D", None, format1)
    worksheet1.set_column("H:H", None, format2)  # evalue in scientific format
    worksheet1.conditional_format("B1:B" + str(len(df)+1), {"type": "cell",
                                                          "criteria": "<",
                                                          "value": (100 - float(intersp_div)),
                                                          "format": format3})
    worksheet1.conditional_format("D1:D" + str(len(df)+1), {"type": "cell",
                                                          "criteria": "<",
                                                          "value": 75,
                                                          "format": format3})
    worksheet1.autofit()

    worksheet2 = writer.sheets["Results_to_check"]
    worksheet2.set_column("B:B", None, format1)  # %id and qcov in float format .1f
    worksheet2.set_column("D:D", None, format1)
    worksheet2.set_column("H:H", None, format2)  # evalue in scientific format
    worksheet2.autofit()

    worksheet3 = writer.sheets["Lineages"]
    worksheet3.set_column("A:A", 50, format4)  # format the first column

    worksheet4 = writer.sheets["Summary"]
    worksheet4.set_column("A:A", 25, format4)
    worksheet4.set_row(3, None, format1)
