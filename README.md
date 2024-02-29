# VirHost
## Desription 
Target: metagenomic-assembled RNA viruses, 
host label lists include 4 eukaryotic kingdom and the prokayotic domain, allow order level classification...

model structure
Outstanding performance in .. experiments ...

## Usage
VirHost is designed to predict the hosts of RNA viruses within known virus orders. To enable VirHost's prediction, the viral taxonomic information is required and its format (.csv) is as below. The first column is the accession numbers in the query fasta file, and the second column show the virus orders of the query viruses.
|   | y\|virus order |
| ------------- | ------------- |
| NC_000858.1 | Ortervirales  |
| NC_019922.1  | Norzivirales  |
| ...  | ...  |

(Optional) For friendly usage, we provided a simple alignment-based method to classify the virus sequences at the order level by BLASTn.
```
python determine_order.py -i test.fasta
```

The input files include the fasta file and the corresponding taxonomic information table. The prediction will be output to {otuput_dir}/result.csv

```
python VirHost.py -i test.fasta --taxa VH_taxa.csv -o VH_result
```

The output format is as:
|   | y\|virus order | pred\|L1 | pred\|L2 | evidence |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| NC_000858.1 | Ortervirales  | Chordata | Primates | pred_high_confidence |
| NC_001426.1  | Norzivirales  | Bacteria |  | assign |
| ...  | ...  | ...  | ...  | ...  |
