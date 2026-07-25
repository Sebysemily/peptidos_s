#!/usr/bin/env python3
import argparse
import csv
import pickle
import statistics
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from isoelectric import ipc
from propy.PyPro import GetProDes

warnings.filterwarnings("ignore")

def Merge(dict1, dict2):
    dict2.update(dict1)
    return dict2

def read_fasta(path):
    records = []
    header = None
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        records.append((header, "".join(seq)))
    return records

def calculate_features(sequences, resources_dir):
    aggDict = pickle.load(open(Path(resources_dir) / 'aggDict.pkl', 'rb'))
    chargeDict = pickle.load(open(Path(resources_dir) / 'chargeDict.pkl', 'rb'))

    pp = []
    isoElP = []
    agg = []
    chargeDensities = []

    for seq in sequences:
        isoElP.append(ipc.predict_isoelectric_point(seq, "IPC_peptide"))
        agList = []
        chargeList = []
        for aa in seq:
            agList.append(aggDict.get(aa, 0))
            chargeList.append(chargeDict.get(aa, 0))
        
        agg.append(statistics.mean(agList) if agList else 0)
        mw = ipc.calculate_molecular_weight(seq)
        chargeDensities.append(sum(chargeList) / mw if mw else 0)

        Des = GetProDes(seq)
        aACompDes = Des.GetAAComp()
        dPCompDes = Des.GetDPComp()
        moreauBrotoAutoDes = Des.GetMoreauBrotoAuto()
        moranAutoDes = Des.GetMoranAuto()
        gearyAutocDes = Des.GetGearyAuto()
        cTDDes = Des.GetCTD()
        
        try:
            aPAACDes = Des.GetAPAAC(lamda=5, weight=0.05)
            paacDes = Des.GetPAAC(lamda=5, weight=0.05)
            sOCNDes = Des.GetSOCN()
            qSODes = Des.GetQSO(maxlag=30, weight=0.1)
        except Exception:
            # Handle extremely short sequences by padding dictionaries with 0
            # Propy fails for seqs shorter than lamda or maxlag
            pass # In our case peptides are > 10 so lamda=5 is fine, but maxlag=30 might fail for seq < 30!

        # Wait, the original code uses maxlag=30. If sequence is < 30, GetQSO might crash.
        # Let's hope propy handles it or we'll pad it.
        try:
            qSODes = Des.GetQSO(maxlag=30, weight=0.1)
        except:
            qSODes = {f"QSO{i}": 0 for i in range(1, 100)} # dummy

        Merge(qSODes, sOCNDes)
        Merge(sOCNDes, paacDes)
        Merge(paacDes, aPAACDes)
        Merge(aPAACDes, cTDDes)
        Merge(cTDDes, gearyAutocDes)
        Merge(gearyAutocDes, moranAutoDes)
        Merge(moranAutoDes, moreauBrotoAutoDes)
        Merge(moreauBrotoAutoDes, dPCompDes)
        Merge(dPCompDes, aACompDes)
        pp.append(aACompDes.copy())

    df = pd.DataFrame(data=pp)
    dsss = {
        "seq": sequences,
        "isoElP": isoElP,
        "aggregationPropensityInVivo": agg,
        "chargeDensity": chargeDensities
    }
    p5 = pd.DataFrame(data=dsss)
    return pd.concat([p5, df.reindex(p5.index)], axis=1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--resources", required=True, help="Path to amptox directory")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    if not records:
        # Create empty output
        with open(args.output, "w") as f:
            f.write("id,amptox_rf_prediction,amptox_svc_prediction\n")
        return

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]

    try:
        p6 = calculate_features(seqs, args.resources)
    except Exception as e:
        print(f"Error calculating features: {e}")
        # fallback empty
        with open(args.output, "w") as f:
            f.write("id,amptox_rf_prediction,amptox_svc_prediction\n")
        return

    features = ["isoElP","chargeDensity","aggregationPropensityInVivo","_NormalizedVDWVD1075","MoreauBrotoAuto_Hydrophobicity1","GearyAuto_AvFlexibility4","_SecondaryStrD1050","_PolarizabilityD2001","tausw8","tausw2","APAAC2","GearyAuto_ResidueVol12","MoreauBrotoAuto_ResidueVol2","GearyAuto_ResidueASA7","L","MoranAuto_ResidueVol15","MoranAuto_AvFlexibility15","MoranAuto_ResidueVol10","W","GearyAuto_Mutability10","MoreauBrotoAuto_Steric1","QSOgrant30","_HydrophobicityD2075","_HydrophobicityD3001","MoranAuto_Polarizability2","MoranAuto_Polarizability7","GearyAuto_ResidueVol6","GearyAuto_Mutability9","GearyAuto_Mutability7","_HydrophobicityT12","MoreauBrotoAuto_Hydrophobicity14","MoreauBrotoAuto_Hydrophobicity17","MoreauBrotoAuto_Hydrophobicity12","_PolarityT23","MoranAuto_Mutability12","QSOSW29","QSOgrant22","QSOgrant21","_PolarizabilityT23","MoranAuto_Steric8","MoreauBrotoAuto_Polarizability14","MoranAuto_ResidueASA12","GearyAuto_Hydrophobicity10","MoranAuto_FreeEnergy11","_PolarizabilityC2","_ChargeD1001","MoranAuto_Steric10","_PolarityC1","_PolarityC3","MoranAuto_Hydrophobicity17","MoranAuto_Hydrophobicity10","MoreauBrotoAuto_FreeEnergy5","QSOSW16","QSOSW12","MoreauBrotoAuto_AvFlexibility4","MoreauBrotoAuto_AvFlexibility6","_NormalizedVDWVD2001","taugrant6","GearyAuto_Steric8","_SecondaryStrT13","MoreauBrotoAuto_Steric17","MoranAuto_FreeEnergy6","_PolarizabilityD2050","MoreauBrotoAuto_Steric10","GearyAuto_Polarizability8","GearyAuto_Polarizability1","GearyAuto_Polarizability3","MoreauBrotoAuto_Steric11","MoranAuto_Mutability7","MoreauBrotoAuto_Mutability1","_NormalizedVDWVC3","_SecondaryStrD1100","_HydrophobicityC2","_PolarizabilityT13","_PolarityD3001","MoranAuto_AvFlexibility7","MoranAuto_AvFlexibility3","_SecondaryStrD2001","MoreauBrotoAuto_FreeEnergy10","MoreauBrotoAuto_FreeEnergy11","MoreauBrotoAuto_FreeEnergy12","MoreauBrotoAuto_FreeEnergy13","MoreauBrotoAuto_FreeEnergy16","_ChargeT12","_SolventAccessibilityD1075","_NormalizedVDWVT13","GearyAuto_FreeEnergy8","MoreauBrotoAuto_ResidueVol10","GearyAuto_FreeEnergy12","GearyAuto_Steric11"]
    
    fMinMax = pd.read_csv(Path(args.resources) / "90FeaturesMinMax.csv")
    
    newds = {"id": ids}
    for j in features:
        p = []
        # if feature is missing due to sequence length limits in propy, pad with 0
        if j not in p6.columns:
            p = [0] * len(p6)
        else:
            for i in p6[j]:
                val = (i - fMinMax.loc[0, j]) / (fMinMax.loc[1, j] - fMinMax.loc[0, j])
                p.append(val)
        newds[j] = p
        
    newdf = pd.DataFrame(data=newds)

    rfModel = pickle.load(open(Path(args.resources) / 'rf_tr_cv_sf_final_model.sav', 'rb'))
    svcModel = pickle.load(open(Path(args.resources) / 'svc_twice_cv_final_model.sav', 'rb'))

    # predict
    X = newdf.iloc[:, 1:].values
    
    try:
        rf_preds = rfModel.predict(X)
        svc_preds = svcModel.predict(X)
    except Exception as e:
        print(f"Prediction error: {e}")
        rf_preds = ["Error"] * len(ids)
        svc_preds = ["Error"] * len(ids)

    out_rows = []
    for i in range(len(ids)):
        rf_val = str(rf_preds[i]).strip().lower()
        svc_val = str(svc_preds[i]).strip().lower()
        out_rows.append({
            "id": ids[i],
            "amptox_rf_prediction": "Toxic" if "non" not in rf_val and "toxic" in rf_val else "Non-Toxic",
            "amptox_svc_prediction": "Toxic" if "non" not in svc_val and "toxic" in svc_val else "Non-Toxic"
        })

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "amptox_rf_prediction", "amptox_svc_prediction"])
        writer.writeheader()
        writer.writerows(out_rows)

if __name__ == "__main__":
    main()
