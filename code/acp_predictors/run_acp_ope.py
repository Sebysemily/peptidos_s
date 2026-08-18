#!/usr/bin/env python3
import os
import re
import sys
import argparse
import joblib
import pandas as pd
import numpy as np
from numpy import array, argmax, linalg as la
from keras.preprocessing.sequence import pad_sequences

def readFasta(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return myFasta

def OE(seq):
    chars = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']
    fea = []
    for i in range(len(seq)):
        if seq[i] =='A': tem_vec = [1]
        elif seq[i]=='C': tem_vec = [2]
        elif seq[i]=='D': tem_vec = [3]
        elif seq[i]=='E' or seq[i]=='U': tem_vec = [4]
        elif seq[i]=='F': tem_vec = [5]
        elif seq[i]=='G': tem_vec = [6]
        elif seq[i]=='H': tem_vec = [7]
        elif seq[i]=='I': tem_vec = [8]
        elif seq[i]=='K': tem_vec = [9]
        elif seq[i]=='L': tem_vec = [10]
        elif seq[i]=='M' or seq[i]=='O': tem_vec = [11]
        elif seq[i]=='N': tem_vec = [12]
        elif seq[i]=='P': tem_vec = [13]
        elif seq[i]=='Q': tem_vec = [14]
        elif seq[i]=='R': tem_vec = [15]
        elif seq[i]=='S': tem_vec = [16]
        elif seq[i]=='T': tem_vec = [17]
        elif seq[i]=='V': tem_vec = [18]
        elif seq[i]=='W': tem_vec = [19]
        elif seq[i]=='X' or seq[i]=='B' or seq[i]=='Z': tem_vec = [20]    
        elif seq[i]=='Y': tem_vec = [21]
        else: tem_vec = [20]
        fea.append(tem_vec)
    return fea

def generateGroupPairs(groupKey):
    gPair = {}
    for key1 in groupKey:
        for key2 in groupKey:
            gPair[key1+'.'+key2] = 0
    return gPair

def CKSAAGP(fastas, gap = 5):
    group = {'alphaticr': 'GAVLMI', 'aromatic': 'FYW', 'postivecharger': 'KRH', 'negativecharger': 'DE', 'uncharger': 'STCPNQ'}
    AA = 'ARNDCQEGHILKMFPSTWYV'
    groupKey = group.keys()
    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key
    gPairIndex = []
    for key1 in groupKey:
        for key2 in groupKey:
            gPairIndex.append(key1+'.'+key2)
    encodings = []
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = []
        for g in range(gap + 1):
            gPair = generateGroupPairs(groupKey)
            sum_pairs = 0
            for p1 in range(len(sequence)):
                p2 = p1 + g + 1
                if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
                    gPair[index[sequence[p1]]+'.'+index[sequence[p2]]] += 1
                    sum_pairs += 1
            if sum_pairs == 0:
                for gp in gPairIndex: code.append(0)
            else:
                for gp in gPairIndex: code.append(gPair[gp] / sum_pairs)
        encodings.append(code)
    return encodings

def TransDict_from_list(groups):
    transDict = dict()
    tar_list = ['0', '1', '2', '3', '4', '5', '6']
    result = {}
    for index, group in enumerate(groups):
        for c in sorted(group):
            result[c] = str(tar_list[index])
    return result

def translate_sequence(seq, TranslationDict):
    from_list = list(TranslationDict.keys())
    to_list = list(TranslationDict.values())
    return seq.translate(str.maketrans("".join(from_list), "".join(to_list)))

def get_3_protein_trids():
    nucle_com = []
    chars = ['0', '1', '2', '3', '4', '5', '6']
    base = len(chars)
    end = base ** 3
    for i in range(end):
        n = i
        ch0 = chars[n % base]
        n = int(n / base)
        ch1 = chars[n % base]
        n = int(n / base)
        ch2 = chars[n % base]
        nucle_com.append(ch0 + ch1 + ch2)
    return nucle_com

def get_4_nucleotide_composition(tris, seq):
    seq_len = len(seq)
    tri_feature = [0] * len(tris)
    k = len(tris[0])
    note_feature = [[0 for cols in range(len(seq) - k + 1)] for rows in range(len(tris))]
    for x in range(len(seq) + 1 - k):
        kmer = seq[x:x + k]
        if kmer in tris:
            ind = tris.index(kmer)
            note_feature[ind][x] += 1
    if len(seq) - k + 1 > 0:
        try:
            u, s, v = la.svd(note_feature)
            for i in range(len(s)):
                tri_feature = tri_feature + u[:, i] * s[i] / seq_len # Changed list append to vector add
        except:
            pass # fallback if SVD fails on empty or weird shapes
    return list(tri_feature)

def prepare_feature_kmer(protein_seqs):
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    kmer = []
    for seq in protein_seqs:
        protein_seq = translate_sequence(seq, group_dict)
        protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq)
        kmer.append(protein_tri_fea)
    return np.array(kmer)

def extract_features(fasta_file):
    fastas = readFasta(fasta_file)
    if not fastas: return [], [], [], []
    
    ids = [f[0] for f in fastas]
    seqs = [f[1] for f in fastas]
    
    # 1. OE
    x_test_oe = [OE(s) for s in seqs]
    # The models in ACP-OPE were trained on sequences up to length 50. 
    # We must pad to 50 for the CNN/BiLSTM to work, as maxlen_test was 50.
    maxlen = 50
    x_test_pad = np.array(pad_sequences(x_test_oe, padding='post', maxlen=maxlen))
    
    # 2. AAC
    hc_AAC_test = np.zeros((len(x_test_oe), 21))
    for row, seq in enumerate(x_test_oe):
        for i in seq:
            hc_AAC_test[row][i[0]-1] += 1/len(seq)
            
    # 3. DPC
    comb = []
    for i in range(1,22):
        for j in range(i,22):
            comb.append((i,j))
    comb_index = {c: idx for idx, c in enumerate(comb)}
    hc_DPC_test = np.zeros((len(x_test_oe), len(comb)))
    for row, seq in enumerate(x_test_oe):
        if len(seq) > 1:
            for i in range(len(seq)-1):
                a = tuple(sorted([seq[i][0], seq[i+1][0]]))
                if a in comb_index:
                    hc_DPC_test[row][comb_index[a]] += 1/(len(seq)-1)
                    
    # 4. CKSAAGP
    hc_CKS_test = np.array(CKSAAGP(fastas))
    
    # 5. K-mer
    kmer_test = prepare_feature_kmer(seqs)
    
    # Combine HC features
    hc_test = np.c_[hc_AAC_test, hc_DPC_test, hc_CKS_test, kmer_test]
    
    # Combine all for LGBM
    x_test_1 = np.squeeze(x_test_pad) if x_test_pad.ndim > 2 else x_test_pad
    # But wait! x_test_pad is shape (N, 50, 1). squeeze makes it (N, 50).
    if x_test_pad.shape[-1] == 1:
        x_test_1 = x_test_pad[:, :, 0]
    else:
        x_test_1 = x_test_pad
        
    X_test_lgbm = np.c_[hc_test, x_test_1]
    
    return ids, seqs, x_test_pad, hc_test, X_test_lgbm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--csv-output", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--mapping", required=False)
    args = parser.parse_args()
    
    ids, seqs, x_test_pad, hc_test, X_test_lgbm = extract_features(args.fasta)
    
    if len(ids) == 0:
        df = pd.DataFrame(columns=["peptide_id", "sequence", "prediction_model1", "score_model1", "prediction_model2", "score_model2"])
        df.to_csv(args.csv_output, index=False)
        return
        
    # Load models
    model1 = joblib.load(os.path.join(args.model_dir, 'bilstm_main7836.joblib'))
    model2 = joblib.load(os.path.join(args.model_dir, 'lgbm_main7865.joblib'))
    
    # Predict Model 1 (BiLSTM)
    pred1_prob = model1.predict([x_test_pad, hc_test]).ravel()
    
    # Fix for LGBM NotFittedError due to version mismatches
    if not hasattr(model2, 'fitted_'):
        model2.fitted_ = True
    if hasattr(model2, '_classes') and not hasattr(model2, 'classes_'):
        model2.classes_ = model2._classes
        
    # Predict Model 2 (LGBM)
    pred2_prob = model2.predict_proba(X_test_lgbm)[:, 1]
    
    records = []
    for i in range(len(ids)):
        score1 = float(pred1_prob[i])
        score2 = float(pred2_prob[i])
        pred1 = "AntiCP" if score1 >= 0.5 else "Non AntiCP"
        pred2 = "AntiCP" if score2 >= 0.5 else "Non AntiCP"
        records.append({
            "indexed_id": ids[i],
            "sequence": seqs[i],
            "score_model1": score1,
            "prediction_model1": pred1,
            "score_model2": score2,
            "prediction_model2": pred2
        })
        
    df = pd.DataFrame(records)
    
    if args.mapping and os.path.exists(args.mapping):
        mapping_df = pd.read_csv(args.mapping)
        df = df.merge(mapping_df[["indexed_id", "peptide_id"]], on="indexed_id", how="left")
        df = df.drop(columns=["indexed_id"])
        cols = ["peptide_id"] + [c for c in df.columns if c != "peptide_id"]
        df = df[cols]
    else:
        df = df.rename(columns={"indexed_id": "peptide_id"})
        cols = ["peptide_id"] + [c for c in df.columns if c != "peptide_id"]
        df = df[cols]
        
    df.to_csv(args.csv_output, index=False)

if __name__ == "__main__":
    main()
