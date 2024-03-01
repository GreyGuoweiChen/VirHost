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

The evidence list has four labels, including **pred_high_confidence**, **pred_low_confidence**, **assign**, **BLASTn**, and **unclassified**. The "pred" prefix represent that the prediction is made by the learning model. If the prediction score is hgiher than the built-in score cutoff, we regard it as a high-confidence prediction and thus give "pred_high_confidence"; otherwise, we will give "pred_low_confidence". The "assign" evidence means that the reference viruses in the same order with the query virus infect a specific host group. So we assign the host group to the virus order without prediction.

Besides, VirHost encodes the query sequences at the protein level. To obtain the protein level representation, we translate the sequences into proteins by prodigal first. For those sequences do not encode proteins, it is challenging to generate a reliable result, and VirHost does not accept these non-coding sequences for higher confidence. However, for user convenience, we try to predict the host of these sequences adopting the best alignment stratedy by BLASTn. We aligned the query virus against our reference database, assign the host as the label of the best hit result, and give "BLASTn" evidence. Finally, we give the "unclassified" evidence to query viruses which do not have order information.
