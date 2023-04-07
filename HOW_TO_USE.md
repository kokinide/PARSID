# PARSID.py

PARSID.py (PArser for Rapid Species IDentification) is a script implemented in Python. It provides an automated solution for the routine process of identifying the species of unidentified sequences of a chosen molecular marker. It covers all the workflow from the generation of the custom database for the local BLAST, to parsing the results and filtering them, and provides an easy-to-consult output in EXCEL format. The filtering is done by the user, by inserting cut-off variables, and by providing if necessary, a supplementary file to integrate specific check warnings.

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

## Files needed:
1. **Input file**: FASTA file containing the sequences with unique identifier (e.g. extraction codes, vouchers, accession numbers). **Note**: the unique sequence identifiers should not include special characters (e.g. spaces, vertical bars). The sequences should not include gaps. Example file:
2. **Reference Sequences Library (RSL)**: single-marker reference database in FASTA format against which the local BLAST is performed. For each lineage of the study group, it contains a representative sequence for each nominal species, bearing taxonomic updated nomenclature, and each candidate species, bearing unique working names. This multi-FASTA file needs to be structured with unique identifiers for each sequence, which are strings containing either taxonomic relevant information (e.g. species name, voucher etc.) or accession numbers if the sequence is obtained from GenBank. **Note**: the unique sequence identifiers should not include special characters (e.g. spaces, vertical bars) that can hamper the functioning of the code. The sequences should not include gaps. Example file:`input_file` `ReferenceSequencesDatabase_16S_test.fasta`
3. **Check_tags file (not mandatory)**: CSV format file (delimiters accepted: comma, semicolon, tab, or space) containing the list of the taxa to tag and manually check later. **Note**: the file should be organized into two columns, “Label” and “Check_tag”. "Labels" can include the whole species name or just a portion (e.g. genus, species, lineage unique code identifier, accession number, etc.), as it is spelled in the RSL. "Check_tag" can include whichever check warning associated to that specific label that is relevat to show next to the molecular identification of the lineage (e.g. check sampling locality, check specimen pictures, etc.). Usually the tagged lineages are the ones for which the sole molecular identification is not enough to reliably assign a sample to a lineage, and in those cases further data need to be checked at a later stage (e.g. collection locality, pictures, etc.). Example file: 


## Step 1. Create the custom database for the local BLAST from the RSL
This is the only step that needs to be run from the Python terminal
```
makeblastdb -in ReferenceSequencesDatabase_16S_test.fasta -parse_seqids -out ref_library -dbtype nucl
```

## Step 2. Local BLAST
The script asks the user the file names of the `input_file` and the RSL `ref_library`, to be provided without file extension. At this stage, the user needs to provide the `cutoff_pident`, which is the cut-off in the BLAST percent sequence identity (under that cut-off the results will not be saved in the local BLAST result).
The BLAST result will include the 12 standard BLAST fields (check `blastn -help` in `outfmt`), plus query length `qlen`, subject length `slen`, number of gaps `gaps`, and query cover `qcovs`.
```
cline = NcbiblastnCommandline(query=input_file + ".fasta", db="ref_library",
                              evalue=0.001, max_target_seqs=150,
                              perc_identity=cutoff_pident, out=input_file + ".txt",
                              outfmt="6 std qlen slen gaps qcovs")
cline()
```
