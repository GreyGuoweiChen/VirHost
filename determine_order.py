import os
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO

VirHost_path = str(os.path.dirname(os.path.abspath(__file__)))

def parse_cmd():
    parser = argparse.ArgumentParser(description="ARGUMENTS")
    parser.add_argument('-i', '--input', type = str, help="The name of query fasta file.")
    parser.add_argument('-o', '--output', default='.', type = str,
                        help="The output directory，including the output and intermediate file.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise Exception("The fasta file does not exist")
    if not os.path.exists(args.output):
        raise Exception("The output directory does not exist")
    return args

def make_blast_db(arg):
    try:
        res = subprocess.check_call(f"which blastn",
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                    shell=True)
    except Exception as e:
        raise Exception("Lack of dependency package in python environment: blastn")

    virus_zip = VirHost_path + "/virus/virus.fasta.zip"
    virus_fasta = VirHost_path + "/virus/virus.fasta"
    if not os.path.exists(virus_fasta):
        print("Unzipping the compressed fasta file ...")
        res = subprocess.check_call(f"unzip {virus_zip} -d {VirHost_path}/virus/", shell=True)

    print("Aligning ...")
    res = subprocess.check_call(
        f"makeblastdb "
        f"-in {VirHost_path}/virus/virus.fasta "
        f"-dbtype nucl "
        f"-out {VirHost_path}/virus/blastn",
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)

def alignment(arg):
    res = subprocess.check_call(  # f"source ~/.bashrc; "
        f"blastn -task blastn -max_target_seqs 10 -max_hsps 1 -evalue 1 -outfmt 6 "
        f"-query {arg.input} "
        f"-db {VirHost_path}/virus/blastn "
        f"-out {arg.output}/VH_blastn.tmp.txt", shell=True)
    print("Finished aligning.")

def read_blast(X, df_train, blast_out):
    print("Parsing the alignment result ...")
    with open(blast_out, 'r') as f:
        for line in f.readlines():
            t = line.strip().split("\t")
            if t[0] == '':  # filter blank row
                break
            score = float(t[-1])
            # if t[0] in X.index and t[1] in df_train.index:  #query in X.index; hit in df_train
            order = df_train.loc[t[1], "y|virus order"]
            if score > X.loc[t[0], "score"] and score > 42:   # score > 42;
                X.loc[t[0], "y|virus order"] = order
                X.loc[t[0], "score"] = score
                # X.loc[t[0], "E"] = float(t[-2])

    return X


if __name__ == "__main__":
    arg = parse_cmd()
    fasta_file, file_output = arg.input, arg.output

    make_blast_db(arg)
    alignment(arg)

    acc_oredered_list = []
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            acc_oredered_list.append(record.id)
    df = pd.DataFrame(index = acc_oredered_list, columns=["y|virus order", "score"])
    df["score"] = 0
    # df["E"] = 1

    df_train = pd.read_csv(f"{VirHost_path}/virus/virus_label.csv", index_col=0)
    df = read_blast(df, df_train, f"{arg.output}/VH_blastn.tmp.txt")
    res = subprocess.check_call( f"rm {arg.output}/VH_blastn.tmp.txt", shell=True)

    df["y|virus order"] = df["y|virus order"].fillna("Unclassified")
    df = df.drop(columns=["score"])
    df.to_csv(f"{arg.output}/VH_taxa.csv")






