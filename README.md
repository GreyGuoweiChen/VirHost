# VirHost: a machine learning-based method for predicting reservoir hosts of RNA viruses through viral genomes
Target: metagenomic-assembled RNA viruses, 
host label lists include 4 eukaryotic kingdom and the prokayotic domain, allow order level classification...

model structure
Outstanding performance in .. experiments ...

Viruses are obligate intracellular parasites that depend on living organisms for their replication and survival. While the application of metagenomic high-throughput sequencing technologies have facilitated the discovery of the viral dark matter, how to determine the reservoir hosts of the metagenome-originated viruses remains challenging owing to the complex composition of the metagenomic sequencing samples. The high genetic diversity of RNA viruses poses a great challenge to the alignment-based methods.

Here, we introduce VirHost, a machine learning-based tool that predicts the reservoirs of RNA viruses solely based on viral genomes. It takes complete RNA viral genomes as input and predicts the natural reservoir host group from kingdom level to order level. 

VirHost is designed as a two-layer classification framework to hierarcically predict the host groups of query viruses. Layer 1 contains 5 branches (Chordata, Invertebrate, Viridiplantae, Fungi, Bacteria), which are categorized into kingdom and phylum level. Layer 2, designed for Chordata subtree, has 10 leaves, which are at the class and order level. VirHost hierarchically predict the host lineage along the tree. To achieve an accuracte host prediction method, it combines the virustaxonomic information, viral genomic traits, and sequences homology. Among 

## Dependency:
* python 3.x
* BLAST 2.12.0
* Prodigal 2.6.3
* xgboost 2.0.3
* pandas 2.0.3
* biopython 1.83
* numpy 1.23.5

### Quick install
We highly recommend using conda to install all the dependencies. To install, please download VirHost by "git clone"
```
git clone https://github.com/GreyGuoweiChen/VirHost.git
cd VirHost

# create the environment and install the dependencies using conda
conda create -n virhost -c bioconda python=3.8 blast=2.12.0 xgboost pandas numpy prodigal biopython numpy

# activate the environment
conda activate virhost
```


## Usage:
To enable VirHost's prediction, the viral taxonomic information is required and its format (.csv) is as below. The first column is the accession numbers in the query fasta file, and the second column show the virus orders of the query viruses.
|   | y\|virus order |
| ------------- | ------------- |
| NC_000858.1 | Ortervirales  |
| NC_019922.1  | Norzivirales  |
| ...  | ...  |

(Optional) For friendly usage, we provided a simple alignment-based method to classify the virus sequences at the order level by BLASTN. However, we recommend users to generate their own taxonomic classification result at the order level, which is expected to improve the prediction confidence.
```
python determine_order.py [-i INPUT_CONTIG] [-o TAXONOMIC_RESULT]
```

The input files include the fasta file and the corresponding taxonomic information table. The prediction will be output to {otuput_dir}/result.csv

```
python VirBot.py [-i INPUT_CONTIG] [--taxa TAXONOMIC_INFORMATION_OF_INPUT] [-o OUTPUT_DIRECTORY]
```

The output format is as:
|   | y\|virus order | pred\|L1 | pred\|L2 | evidence |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| NC_000858.1 | Ortervirales  | Chordata | Primates | pred_high_confidence |
| NC_001426.1  | Norzivirales  | Bacteria |  | assign |
| ...  | ...  | ...  | ...  | ...  |

The evidence list has 5 labels, including **pred_high_confidence**, **pred_low_confidence**, **assign**, **BLASTn**, and **unclassified**. The "pred" prefix represent that the prediction is made by the learning model. If the prediction score is hgiher than the built-in score cutoff, we regard it as a high-confidence prediction and thus give "pred_high_confidence"; otherwise, we will give "pred_low_confidence". The "assign" evidence means that the reference viruses in the same order with the query virus infect a specific host group. So we assign the host group to the virus order without prediction.

Besides, VirHost encodes the query sequences at the protein level. To obtain the protein level representation, we translate the sequences into proteins by prodigal first. For those sequences do not encode proteins, it is challenging to generate a reliable result, and VirHost does not accept these non-coding sequences for higher confidence. However, for user convenience, we try to predict the host of these sequences adopting the best alignment stratedy by BLASTn. We aligned the query virus against our reference database, assign the host as the label of the best hit result, and give "BLASTn" evidence. Finally, we give the "unclassified" evidence to query viruses which do not have order information.

### Example
```
# create the taxonomic file first
python determine_order.py -i test.fasta -o VH_taxa.csv

# predict the host of query viruses
python virhost.py -i test.fasta --taxa VH_taxa.csv -o VH_result
```



