# PARSID.py

PARSID.py (PArser for Rapid Species IDentification) is a script implemented in Python. It provides an automated solution for the molecular species identification of sequences of a chosen genetic marker. It covers the workflow of the generation from the creation of the custom database for the local BLAST, to parsing the results and filtering them, and provides an easy-to-consult output in EXCEL format. The filtering is controlled by the user, by inputting cut-off variables, and by providing if necessary, a supplementary file to integrate specific check warnings.

## Packages needed:
```
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from Bio import SeqIO
import xlsxwriter
import pandas as pd
from operator import itemgetter
import os.path
```

## Files needed
1. **Input file**: FASTA file containing the sequences with unique identifiers (e.g. extraction codes, vouchers, accession numbers). **Note**: the unique sequence identifiers should not include special characters (e.g. spaces, vertical bars). The sequences should not include gaps. Example file:`Input16S_test.fasta`
2. **Reference Sequences Library (RSL)**: single-marker multi-FASTA file against which the local BLAST of the input file is performed. For each lineage of the study group, it contains a representative sequence for each nominal species bearing taxonomic updated nomenclature, and each candidate species (if any) bearing unique working names. This multi-FASTA file needs to be structured with unique identifiers for each sequence, which are strings containing either taxonomic relevant information (e.g. species name, voucher etc.) or accession numbers if the sequence is obtained from GenBank. **Note**: the unique sequence identifiers should not include special characters (e.g. spaces, vertical bars) that can hamper the functioning of the code (especially Step 1). The sequences should not include gaps. Example file: `ReferenceSequencesDatabase_16S_test.fasta`
3. **Check_tags file (not mandatory)**: CSV format file (delimiters accepted: comma, semicolon, tab, or space) containing the list of the taxa to tag and manually check later. **Note**: the file should be organized into two columns, “Label” and “Check_tag”. "Label" can include the whole species name or just a portion (e.g. genus, species, lineage unique code identifier, accession number, etc.), as it is spelled in the RSL. "Check_tag" can include whichever check warning associated to that specific label, that will show next to the molecular identification of the lineage (e.g. check sampling locality, check specimen pictures, etc.). For the tagged lineages the molecular identification alone is not enough to reliably assign a sample to a lineage, and in those cases further data need to be checked at a later stage (e.g. collection locality, pictures, etc.). Example file:`Check_tags.csv`


## Step 1. Create the custom database for the local BLAST from the RSL
This is the only step that needs to be run from the Python terminal.
```
makeblastdb -in ReferenceSequencesDatabase_16S_test.fasta -parse_seqids -out ref_library -dbtype nucl
```
At the end of this step, the local BLAST indices should appear in the folder.

## Step 2. Local BLAST
The script asks the user the file names of the `input_file` and the RSL `ref_library`, to be provided without file extension. At this stage, the user needs to provide the `cutoff_pident`, which is the cut-off in the BLAST percent sequence identity (under that value the results will not be saved in the local BLAST result).
The BLAST result will include the 12 standard BLAST fields (check `blastn -help` in `outfmt`), plus query length `qlen`, subject length `slen`, number of gaps `gaps`, and query coverage `qcovs`.
```
cline = NcbiblastnCommandline(query=input_file + ".fasta", db="ref_library",
                              evalue=0.001, max_target_seqs=150,
                              perc_identity=cutoff_pident, out=input_file + ".txt",
                              outfmt="6 std qlen slen gaps qcovs")
cline()
```
The result is a .txt file containing the local BLAST hits.

## Step 3. Parsing the local BLAST result
At the beginning, the parse object is created.
```
QueryResults = SearchIO.parse(input_file + ".txt", "blast-tab", fields=my_fields)
```
For each query, the query coverage is calculated. **Note**: in Python the index starts from 0.
`hsp.query_end` is the ending position of the HSP
`hsp.query_start` is the starting position of the HSP
`QueryResult.seq_len` is the length of the query
```
qcov = (hsp.query_end - hsp.query_start) / QueryResult.seq_len * 100
```
During the parsing process, the hits from each query are filtered by (query coverage >= 75%) and (subject length > 90 % of the alignment length) and sorted by percentage of identity, so only the first result is saved. 
At the end of the parsing, this message is printed:
```
# Parsing completed!
```

## Step 4. Checking the parsed result
The script checks if the file "Check_tags.csv" exists, and integrate the information to the parsed result in the "notes" column.
At this stage, the user needs to provide the `intersp_div`, which is the percentage of interspecific sequence divergence. Under the sequences similarity threshold (calculated as 100 - `intersp_div`), the results will be tagged as "to check" (for example, if the `intersp_div` is 3%, the results below 97% will be tagged as "to check"). The results with a `qcov`<75% will be also tagged as "to check".
The results tagged as "to check" are stored to be printed in a separate sheet "Results_to_check" in the final EXCEL file.

## Step 5. Saving final output in EXCEL format
The final output is saved in a easy-to-consult spreadsheet. The final part of the code includes output file formatting.
```
with pd.ExcelWriter(input_file + ".xlsx") as writer:
    df.to_excel(writer, index=False, sheet_name="Best_blast")
    df_to_check.to_excel(writer, index=False, sheet_name="Results_to_check")
    df_lineage.to_excel(writer, index=False, sheet_name="Lineages")
    df_stats.to_excel(writer, index=False, sheet_name="Summary")
```
The sheets are:
1. Best_blast: includes the list of samples with the respective percentage of identity, molecular species identification, percentage of query coverage, query length, subject length, alignment length, evalue, notes. In this sheet, the values of percentage identity below the similarity threshold are highlighted in red. 
2. Results_to_check: includes the full list of hits for the samples that were tagged as "to check" in the notes. The fields are the same as the sheet "Best_blast", but without the notes.
3. Lineages: includes the list of all the lineages in the RSL with the samples found for each. The last field is "Results %id <" + `cutoff_pident`, that lists all the samples below the user's value.
4. Summary: includes general information about the inputs. The stats are:

   | Stat  | Definition |
   | --- | --- |
   | Total taxa  | total number of taxa in the RSL |
   | Taxa without samples  | number of taxa in the RSL without sequences |
   | Taxonomic coverage  | percentage of the RSL covered by the sequences contained by the input file |
   | Total input sequences  | total number of sequences in the input file |
   | Sequences processed (%)  | percentage of the sequences analyzed by the script (ideally is 100%) |
   | N results %id<`cutoff_pident` | number of results below the user's value `cutoff_pident` |
 
