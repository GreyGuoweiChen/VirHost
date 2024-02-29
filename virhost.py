import os
import subprocess
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle as pkl
from collections import defaultdict
from copy import deepcopy
import warnings
warnings.filterwarnings('ignore')

VirHost_path = str(os.path.dirname(os.path.abspath(__file__)))
order_pred_list = pd.read_csv(f"{VirHost_path}/model/order_pred.csv", index_col=0).index.to_list()
df_order_assign = pd.read_csv(f"{VirHost_path}/model/order_assign.csv", index_col=0)
order_assign_list = df_order_assign.index.to_list()

########################################################################################################
nu_counting_dict = {a: 0 for a in ["A", "T", "C", "G"]}
dinu_counting_dict = {a + b: 0 for a in ["A", "T", "C", "G"] for b in ["A", "T", "C", "G"]}  # 3*16 =48
codon_table = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
               "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
               "AAT": "N", "AAC": "N",
               "GAT": "D", "GAC": "D",
               "TGT": "C", "TGC": "C",
               "CAA": "Q", "CAG": "Q",
               "GAA": "E", "GAG": "E",
               "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
               "CAT": "H", "CAC": "H",
               "ATG": "M",
               "ATT": "I", "ATC": "I", "ATA": "I",
               "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L",
               "AAA": "K", "AAG": "K",
               "TTT": "F", "TTC": "F",
               "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
               "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
               "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
               "TGG": "W",
               "TAT": "Y", "TAC": "Y",
               "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
               "TAA": " ", "TGA": " ", "TAG": " "}
reversed_codon_table = defaultdict(list)
for codon, aa in codon_table.items():
    reversed_codon_table[aa].append(codon)
codon_counting_dict = {key: 0 for key, value in codon_table.items()}  # 64
aa_counting_dict = {aa: 0 for aa, codon_list in reversed_codon_table.items()}
feature_count_dict = {"nu": nu_counting_dict,
                      "dinu": deepcopy(dinu_counting_dict), "dinu_bdg": deepcopy(dinu_counting_dict), "dinu_non_bdg": deepcopy(dinu_counting_dict),
                      "codon": codon_counting_dict,
                      "aa": aa_counting_dict}
feature_dict = {"nu": nu_counting_dict,
                "dinu": deepcopy(dinu_counting_dict), "dinu_bdg": deepcopy(dinu_counting_dict), "dinu_non_bdg": deepcopy(dinu_counting_dict),
                "codon_bias": codon_counting_dict,
                "aa_bias": aa_counting_dict} ###### wrong: amino_acid_pair_counting_dict
feat_name = []
for key in feature_dict.keys():
    feat_name += [key + "|" + name for name in list(feature_dict[key].keys())]
# feat_names = 137 genomic traits;

########################################################################################################
def check_dependancy():
    try:
        res = subprocess.check_call(f"which blastn",
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                    shell=True)  # f"source ~/.bashrc; "
    except Exception as e:
        raise Exception("Lack of dependency package in python environment: blastn")

    try:
        res = subprocess.check_call(f"which prodigal",
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                    shell=True)  # f"source ~/.bashrc; "
    except Exception as e:
        raise Exception("Lack of dependency package in python environment: prodigal")

def parse_cmd():
    parser = argparse.ArgumentParser(description = "ARGUMENTS")
    parser.add_argument('-i', '--input', type = str, help = "The query fasta file.")
    parser.add_argument('--taxa', default="VH_taxa.csv", type = str, help = "The virus order taxa of query sequences. (.csv)")
    parser.add_argument('-o', '--output', default = "VH_result", type = str,
                        help="The output directory，including the output and intermediate file.")
    # parser.add_argument('--layer', default = 2, type = int,
    #                     help="The deepest layer of the host phylogenetic tree.")
    args = parser.parse_args()
    return args


def check_cmd(args):
    file_seq = args.input
    if not os.path.exists(file_seq):
        raise Exception("The input fasta file do not exist.")

    file_taxa = args.taxa
    if not os.path.exists(file_taxa):
        raise Exception("The taxonomic information is required.")

    try:    # check seq number of fasta and file_taxa
        num = 0
        with open(file_seq) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                num +=1
        df_taxa = pd.read_csv(file_taxa, index_col=0)
        if num != len(df_taxa):
            raise Exception(f"The number of input sequences do not equal to the number of taxonomic information\n"
                            f"Fasta file includes {num}\n"
                            f"Taxonomic file includes {len(df_taxa)} records")
    except Exception as e:
        print(e)
        exit()

    output_dir = args.output
    if os.path.exists(output_dir):
        raise Exception(f"The output directory ({output_dir}) already exists.")


########################################################################################################
def encoding(file_seq, df_taxa, output_dir):
    def write_acc2fasta(fasta_seq, fasta_out, acc):
        all_records = []
        with open(fasta_seq) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in acc:
                    all_records.append(record)

        SeqIO.write(all_records, fasta_out, "fasta")

    def align_query(file_blastn_db, fileout, blastn_out):
        res = subprocess.check_call( #f"source ~/.bashrc; "
                                     f"blastn -task blastn -max_target_seqs 5000 -max_hsps 1 -evalue 0.001 -outfmt 6 "
                                     f"-query {fileout} "
                                     f"-db {file_blastn_db} "
                                     f"-out {blastn_out}", shell=True)

    def read_blast(X, df_train, blast_out, taxa2onehot, layerid):
        taxa_list = list(taxa2onehot.keys())

        with open(blast_out, 'r') as f:
            for line in f.readlines():
                t = line.strip().split("\t")
                if t[0] == '':  # filter blank row
                    break
                # elif t[0] == t[1]:   # query != hit # do not used when predicting!!
                #     continue
                score = float(t[-1])
                if t[0] in X.index and t[1] in df_train.index:  #query in X.index; hit in df_train
                    host = df_train.loc[t[1], f"y|L{layerid}"]
                    if score >= 10 and score  > X.loc[t[0], host]:   # score >= 10; best score for each group
                        X.loc[t[0], host] = score

        # X[taxa_list] = 1
        X = X.apply(np.log10)
        taxa_dict = {host: "blast|"+host for host in taxa_list}
        X = X.rename(columns = taxa_dict)

        return X

    def SH_encode(order, fasta_order, temp_dir, blastn_out):

        taxa2layerid = {"Biota": 0, "Chordata": 1}
        host_tree_dict, leaves = pkl.load(open(VirHost_path + f"/model/{order}/host_info.pkl", "rb"))
        df_train = pd.read_csv(VirHost_path + f"/model/{order}/data_label.csv", index_col=0)

        M, acc_oredered_list = [], []
        with open(fasta_order, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                acc_oredered_list.append(record.id)

        if 1: # predicte to layer 1
            clf_nodes = (0, "Biota", 1, host_tree_dict["Biota"])
            taxa_dict = {host: "blast|" + host for host in clf_nodes[3]}
            taxa2onehot = {taxa: onehot for onehot, taxa in enumerate(clf_nodes[3])}
            X_test_tmp = pd.DataFrame(1, index=acc_oredered_list, columns=clf_nodes[3])
            X_test_tmp = read_blast(X_test_tmp, df_train, blastn_out, taxa2onehot, clf_nodes[2])
            X_test_tmp.to_csv(temp_dir + f"/X_SH_L1.csv")

        if 2 and "Chordata" in host_tree_dict: # predicte to layer 2

            clf_nodes = (1, "Chordata", 2, host_tree_dict["Chordata"])
            taxa_dict = {host: "blast|" + host for host in clf_nodes[3]}
            taxa2onehot = {taxa: onehot for onehot, taxa in enumerate(clf_nodes[3])}
            df_train = df_train.loc[df_train[f"y|L{clf_nodes[2]}"].isin(clf_nodes[3])]
            X_test_tmp = pd.DataFrame(1, index=acc_oredered_list, columns=clf_nodes[3])
            X_test_tmp = read_blast(X_test_tmp, df_train, blastn_out, taxa2onehot, clf_nodes[2])
            X_test_tmp.to_csv(temp_dir + f"/X_SH_L2.csv")

    def GT_encode(fasta_order, temp_dir):
        def count_coding_seq(coding_seq, featurn_count):
            for i in range(0, len(coding_seq), 3):
                codon = coding_seq[i:i + 3]
                for nu in codon:
                    if nu in nu_counting_dict:
                        featurn_count["nu"][nu] += 1
                if codon in codon_table:
                    featurn_count["codon"][codon] += 1
                    featurn_count["aa"][codon_table[codon]] += 1
                if coding_seq[i: i + 2] in dinu_counting_dict:
                    featurn_count["dinu"][coding_seq[i: i + 2]] += 1
                    featurn_count["dinu_non_bdg"][coding_seq[i: i + 2]] += 1
                if coding_seq[i + 1: i + 3] in dinu_counting_dict:
                    featurn_count["dinu"][coding_seq[i + 1: i + 3]] += 1
                    featurn_count["dinu_non_bdg"][coding_seq[i + 1: i + 3]] += 1
                if i + 6 <= len(coding_seq):
                    if coding_seq[i + 2: i + 4] in dinu_counting_dict:
                        featurn_count["dinu"][coding_seq[i + 2: i + 4]] += 1
                        featurn_count["dinu_bdg"][coding_seq[i + 2: i + 4]] += 1
            return featurn_count

        def transform_feature(feature_count):
            feature = deepcopy(feature_dict)

            # frequency of nucleotide [0:3]
            sum_nu = sum([count for nu, count in feature_count["nu"].items()])
            feature["nu"] = {nu: count / sum_nu for nu, count in feature_count["nu"].items()}

            # prob_i_j / (prob_i * prob_j) in the coding region (16 [4:19] 15:CG)
            sum_dinu = sum([count for dinu, count in feature_count["dinu"].items()])

            feature["dinu"] = {dinu: count / sum_dinu / (feature["nu"][dinu[0]] * feature["nu"][dinu[1]])
            if (feature["nu"][dinu[0]] * feature["nu"][dinu[1]]) != 0 else 0
                               for dinu, count in feature_count["dinu"].items()}

            # prob_i_j / (prob_i * prob_j) in the codon bridge region (16 [20: 35] 31:CG)
            sum_dinu_bdg = sum(feature_count["dinu_bdg"].values())
            feature["dinu_bdg"] = {dinu: count / sum_dinu_bdg / (feature["nu"][dinu[0]] * feature["nu"][dinu[1]])
            if (feature["nu"][dinu[0]] * feature["nu"][dinu[1]]) != 0 else 0
                                   for dinu, count in feature_count["dinu_bdg"].items()}

            # prob_i_j / (prob_i * prob_j) in the codon non_bridge region (16 [36: 41])
            sum_dinu_non_bdg = sum(feature_count["dinu_non_bdg"].values())
            feature["dinu_non_bdg"] = {
                dinu: count / sum_dinu_non_bdg / (feature["nu"][dinu[0]] * feature["nu"][dinu[1]])
                if (feature["nu"][dinu[0]] * feature["nu"][dinu[1]]) != 0 else 0
                for dinu, count in feature_count["dinu_non_bdg"].items()}

            # codon bias among the coding amino acid (64 codon [42: 105])
            feature["codon_bias"] = {codon: count / feature_count["aa"][codon_table[codon]]
            if feature_count["aa"][codon_table[codon]] != 0 else 0
                                     for codon, count in feature_count["codon"].items()}

            # AA frequency (21 amino acid [106:136])
            sum_aa = sum(feature_count["aa"].values())
            feature["aa_bias"] = {aa: count / sum_aa for aa, count in feature_count["aa"].items()}

            feature_array = []
            for domain, domain_feature in feature.items():
                feature_array += list(domain_feature.values())
            feature_array = np.array(feature_array)
            feature_array[feature_array == 0] = 0.0001
            feature_array = np.log2(feature_array)
            # print(len(feature_array))
            return feature_array

        def parse_seq(record, cds_list):
            # coding each nucleotide record
            feature_count = deepcopy(feature_count_dict)
            rc_record = record.reverse_complement()
            l = len(record.seq)

            for start, stop, strand in cds_list:
                if strand == '+':
                    # pr(start, stop, strand , Seq(record.seq[start-1:stop]))
                    cds = record.seq[start - 1:stop]
                    feature_count = count_coding_seq(cds, feature_count)
                elif strand == '-':
                    # pr(start, stop, strand, rc_record.seq[L-stop:L-start+1])
                    cds = rc_record.seq[l - stop: l - start + 1]
                    feature_count = count_coding_seq(cds, feature_count)
            if len(cds_list) > 0:
                feature = transform_feature(feature_count)
            else:
                feature = np.full([137], np.nan)
                # print(record.id, "No protein")
            return feature

        file_cds = temp_dir + "/seq_cds.gff"
        res = subprocess.check_call(f"prodigal -p meta -f gff -i {fasta_order} -o {file_cds}",
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)  # f"source ~/.bashrc; "

        pdg_dict = defaultdict(list)
        with open(file_cds, "r") as cds_handle:
            for line in cds_handle.readlines():
                if line.startswith('#'):
                    continue
                id, source, type, start, end, score, strand, phase, attributes = line.strip().split()
                pdg_dict[id].append((int(start), int(end), strand))

        M, acc_oredered_list, acc_without_protein = [], [], []
        with open(fasta_order, "r") as fasta_handle:
            i = 0
            for record in SeqIO.parse(fasta_handle, "fasta"):
                acc_oredered_list.append(record.id)
                try:
                    feature = parse_seq(record, pdg_dict[record.id])
                except Exception as e:  # report error
                    print(i, record.id, e)
                    feature = np.full([137], np.nan)
                if np.isnan(feature[0]):  # no coding region
                    acc_without_protein.append(record.id)
                M.append(feature)
                i += 1
            M = np.array(M)
        # print(M.shape, len(feat_name))  # 137, 137
        virus_feat_df = pd.DataFrame(M, columns=feat_name, index=acc_oredered_list)
        virus_feat_df.to_pickle(temp_dir + f"/X_GT.pkl")
        return acc_without_protein

    def encode_by_order(file_seq, order, acc):
        temp_dir = f"{output_dir}/tmp/{order}"
        os.makedirs(temp_dir)

        # #write fasta for the order
        fasta_order = temp_dir + f"/virus.fasta"
        write_acc2fasta(file_seq, fasta_order, acc)

        # #align sequences
        file_blastn_db = VirHost_path + f"/model/{order}/blast_db"
        blastn_out = temp_dir + f"/blastn.txt"
        align_query(file_blastn_db, fasta_order, blastn_out)

        # #encode sequences homology
        SH_encode(order, fasta_order, temp_dir, blastn_out)

        # #encode genomic traits
        acc_without_protein = GT_encode(fasta_order, temp_dir)

        return acc_without_protein


    order_list = df_taxa["y|virus order"].unique()
    acc_without_protein = []
    unaccepted_order = []
    for order in order_list:
        if order in order_pred_list: # # encode for predicting
            tmp_acc_without_protein = encode_by_order(file_seq, order, df_taxa[(df_taxa[f"y|virus order"] == order)].index.tolist())
                # align # parse_alignment # feature encoding # counting seq without protein # save
            acc_without_protein.extend(tmp_acc_without_protein)
        elif order in order_assign_list:   # # assign >> do not need any encoding
            pass
        else:
            unaccepted_order.append(order)

    return acc_without_protein, unaccepted_order


########################################################################################################
def predict(df_taxa, output_dir, acc_without_protein):
    def predict_order(order, output_dir, acc_without_protein):
        temp_dir = f"{output_dir}/tmp/{order}"
        host_tree_dict, leaves = pkl.load(open(VirHost_path + f"/model/{order}/host_info.pkl", "rb"))

        X_GT = pd.read_pickle(temp_dir + f"/X_GT.pkl")

        # y_pred
        acc_without_prots = X_GT[X_GT.index.isin(acc_without_protein)].index
        y_without_prots = pd.DataFrame("Unclassified", columns=["pred|L1", "pred|L2", "evidence"], index = acc_without_prots)
        acc_with_prots = X_GT[~(X_GT.index.isin(acc_without_protein))].index
        y_pred = pd.DataFrame(columns=["pred|L1", "pred|L2", "evidence"], index = acc_with_prots)

        # X consists of acc in the corresponding order, including sequencing encoding proteins and do not.
        if 1:
            clf_nodes = (0, "Biota", 1, host_tree_dict["Biota"])
            onehot2taxa = {onehot: taxa for onehot, taxa in enumerate(clf_nodes[3])}

            if len(host_tree_dict["Biota"]) == 1:
                y_pred[["pred|L1","evidence"]] = [onehot2taxa[0], "pred_high_confidence"]
                if len(y_without_prots) > 0:
                    y_without_prots[["pred|L1", "evidence"]] = [onehot2taxa[0], "pred_high_confidence"]
            else:
                X_SH = pd.read_csv(temp_dir + f"/X_SH_L{clf_nodes[2]}.csv", index_col=0)
                X = X_GT.merge(X_SH, left_index=True, right_index=True, suffixes=(None, None))
                X = X.loc[X.index.isin(acc_with_prots)]

                clf = pkl.load(open(VirHost_path + f"/model/{order}/xgb_sub_bl_L{clf_nodes[2]}.model", 'rb'))
                y_prob = clf.predict_proba(X)
                y_hat = np.argmax(y_prob, axis=1)
                y_pred[f"pred|L{clf_nodes[2]}"] = list(map(lambda i: onehot2taxa[i], y_hat))
                # y_hat = clf.predict(X)

                y_prob = np.max(y_prob, axis=1)
                cutoff = pkl.load(open(VirHost_path + f"/model/{order}/cutoff.pkl", "rb"))
                y_pred["evidence"] = list(map(lambda p: "pred_high_confidence" if p > cutoff else "pred_low_confidence", y_prob))

                if len(y_without_prots) > 0:
                    y_without_prots["pred|L1"] = X_SH.loc[X_SH.index.isin(acc_without_prots)].idxmax(axis = 1,)
                    y_without_prots["pred|L1"] = y_without_prots["pred|L1"].map(lambda label: label[6:]) # remove "blast|"
                    y_without_prots["evidence"] = ["BLASTn"]

        if 2 and "Chordata" in host_tree_dict:
            clf_nodes = (1, "Chordata", 2, host_tree_dict["Chordata"])
            onehot2taxa = {onehot: taxa for onehot, taxa in enumerate(clf_nodes[3])}
            acc_L2 = y_pred.loc[y_pred["pred|L1"] == "Chordata"]

            if len(host_tree_dict["Biota"]) == 1:
                y_pred.loc[(y_pred["pred|L1"] == "Chordata"), f"pred|L{clf_nodes[2]}"] = onehot2taxa[0]
            else:
                X_SH = pd.read_csv(temp_dir + f"/X_SH_L{clf_nodes[2]}.csv", index_col=0)
                X = X_GT.merge(X_SH, left_index=True, right_index=True, suffixes=(None, None))
                X = X.loc[(X.index.isin(acc_with_prots))] # shape: [all_seq in this order, 137] > [acc_with_prots, 137]
                X = X.loc[(y_pred["pred|L1"] == "Chordata")] # shape: [acc_with_prots, 137] > ["Chordata", 137]
                if len(X) > 0:
                    clf = pkl.load(open(VirHost_path + f"/model/{order}/xgb_sub_bl_L{clf_nodes[2]}.model", 'rb'))
                    y_hat = clf.predict(X)
                    y_pred.loc[(y_pred["pred|L1"] == "Chordata"), f"pred|L{clf_nodes[2]}"] = list(map(lambda i: onehot2taxa[i], y_hat))


        return y_pred, y_without_prots

    taxa2layerid = {"Biota": 0, "Chordata": 1}

    df_taxa = pd.concat([df_taxa, pd.DataFrame(columns=["pred|L1", "pred|L2", "evidence"])], sort=False)
    order_list = df_taxa["y|virus order"].unique()

    for order in order_list:
        if order in order_pred_list:  # predict
            y_pred,  y_without_prots = predict_order(order, output_dir, acc_without_protein)
            if len(y_pred) > 0:
                df_taxa.loc[df_taxa.index.isin(y_pred.index), ["pred|L1", "pred|L2", "evidence"]] = \
                    y_pred[["pred|L1", "pred|L2", "evidence"]]
            if len(y_without_prots) > 0:
                df_taxa.loc[df_taxa.index.isin(y_without_prots.index), ["pred|L1", "pred|L2", "evidence"]] = \
                    y_without_prots[["pred|L1", "pred|L2", "evidence"]]
        elif order in order_assign_list: #order in [""]: # assign
            df_taxa.loc[(df_taxa["y|virus order"] == order), "pred|L1"] = df_order_assign.loc[order, "host"]
            df_taxa.loc[(df_taxa["y|virus order"] == order), "evidence"] = "assign"
        else:
            pass
    return df_taxa


########################################################################################################
if __name__ == "__main__":
    check_dependancy()
    args = parse_cmd()
    check_cmd(args)


    file_seq, file_taxa, output_dir = args.input, args.taxa, args.output
    os.makedirs(output_dir)
    temp_dir = f"{output_dir}/tmp"
    os.makedirs(temp_dir)

    print("Encoding the feature...")
    df_taxa = pd.read_csv(file_taxa, index_col=0)
    acc_without_protein, unaccepted_order = encoding(file_seq, df_taxa, output_dir)
    print("Done!")

    if len(acc_without_protein) > 0:
        print(f"There are {len(acc_without_protein)} non-coding sequences")
        with open(output_dir + "/acc_without_prot.txt", "w") as output:
            output.write('\n'.join(acc_without_protein))
    if len(unaccepted_order) > 0:
        print(f"The following orders are not accpeted :", unaccepted_order)
        num = len(df_taxa[df_taxa["y|virus order"].isin(unaccepted_order)])
        print(f"There are {num} sequences outside the known virus orders")

    print("Predicting the hosts...")
    df_pred = predict(df_taxa, output_dir, acc_without_protein)
    df_pred.to_csv(f"{output_dir}/result.csv")
    print("Done!")
    print(f"Output the prediction to {output_dir}/result")


