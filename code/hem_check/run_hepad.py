#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


DATASET_MODEL = {
    "Hmp1": "1",
    "Hmp2": "2",
    "Hmp3": "3",
    "Hmpm": "4",
}


def ensure_hepad_resource_links(hepad_root, workdir):
    for subdir in ("iFeaturedata", "PFeaturedata"):
        source = hepad_root / "devPackage" / subdir
        target = workdir / subdir
        if not source.is_dir():
            raise SystemExit(f"Missing HEPAD resource directory: {source}")
        if target.exists():
            continue
        target.symlink_to(source, target_is_directory=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run HEPAD hemolytic peptide prediction on a FASTA file."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "--dataset",
        required=True,
        choices=sorted(DATASET_MODEL),
        help="HEPAD training dataset/model set",
    )
    parser.add_argument(
        "--hepad-root",
        required=True,
        help="Path to the HEPAD repository checkout",
    )
    parser.add_argument(
        "--workdir",
        required=True,
        help="Working directory for intermediate and raw outputs",
    )
    args = parser.parse_args()

    dataset = args.dataset
    hepad_root = Path(args.hepad_root).resolve()
    workdir = Path(args.workdir).resolve()
    input_fasta = Path(args.input).resolve()

    main_program = hepad_root / "mainProgram"
    if not (main_program / "main_predict.py").exists() and not (
        hepad_root / "userPackage" / "Package_HEPAD.py"
    ).exists():
        raise SystemExit(f"HEPAD checkout is incomplete: {hepad_root}")

    encoded_dir = workdir / "mlData" / "new_data"
    output_dir = workdir / "output"
    encoded_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    path_dict = {
        "paramPath": str(hepad_root / "data" / "param" / dataset),
        "saveCsvPath": str(encoded_dir),
        "modelPath": str(hepad_root / "data" / "finalModel" / dataset),
        "outputPath": str(output_dir),
    }

    for label, path in path_dict.items():
        if label == "saveCsvPath":
            continue
        if not Path(path).exists():
            raise SystemExit(f"Missing HEPAD {label} for {dataset}: {path}")

    ensure_hepad_resource_links(hepad_root, workdir)

    sys.path.insert(0, str(hepad_root))
    sys.path.insert(0, str(workdir))
    os.chdir(main_program)

    from userPackage.Package_HEPAD import HEPAD_Predict

    class FixedHEPADPredict(HEPAD_Predict):
        def featureEncode(self):
            from userPackage.Package_Encode import EncodeAllFeatures
            import pandas as pd

            encode_obj = EncodeAllFeatures()
            encode_obj.dataEncodeSetup(
                loadJsonPath=(
                    f'{self.pathDict["paramPath"]}/{self.featureTypeDictJson}'
                )
            )
            encode_obj.dataEncodeOutput(dataList=self.dataList)
            test_df = encode_obj.dataNormalization(
                loadNmlzScalerPklPath=(
                    f'{self.pathDict["paramPath"]}/{self.nmlzPkl}'
                )
            )
            feature_df = pd.read_csv(self.featureRankCsv)
            feature_list = feature_df["feature name"].to_list()
            if self.model_use == "1":
                test_df["MotifBitVec_KKG"] = 0
            elif self.model_use == "2":
                test_df["y"] = 0
            test_df = test_df[feature_list]
            test_df.to_csv(
                f'{self.pathDict["saveCsvPath"]}/test_F{self.featureNum}.csv'
            )

        def doPredict(self):
            import os

            import pandas as pd
            from pycaret.classification import load_model

            data_test_df = pd.read_csv(
                f'{self.pathDict["saveCsvPath"]}/test_F{self.featureNum}.csv',
                index_col=[0],
            ).fillna(0.0)
            pred_vector_list = []
            prob_vector_list = []
            for model_name in self.modelNameList:
                model_path = os.path.join(
                    self.pathDict["modelPath"],
                    f"{model_name}_final",
                )
                pipeline = load_model(model_path)
                if model_name == "lightgbm":
                    booster = pipeline.named_steps["trained_model"]._Booster
                    probabilities = booster.predict(data_test_df.values)
                    if getattr(probabilities, "ndim", 1) > 1:
                        probabilities = probabilities[:, 1]
                    predictions = (probabilities >= 0.5).astype(int)
                else:
                    predictions = pipeline.predict(data_test_df)
                    probabilities = pipeline.predict_proba(data_test_df)[:, 1]
                pred_vector_list.append(predictions)
                prob_vector_list.append(probabilities)

            self.predVectorListIndp = pred_vector_list
            self.probVectorListIndp = prob_vector_list
            self.predVectorDf = pd.DataFrame(
                pred_vector_list,
                index=self.modelNameList,
                columns=data_test_df.index,
            ).T
            self.probVectorDf = pd.DataFrame(
                prob_vector_list,
                index=self.modelNameList,
                columns=data_test_df.index,
            ).T
            self.predVectorDf.to_csv(
                f'{self.pathDict["outputPath"]}/binary_vector.csv'
            )
            self.probVectorDf.to_csv(
                f'{self.pathDict["outputPath"]}/probability_vector.csv'
            )

    predictor = FixedHEPADPredict(
        model_use=DATASET_MODEL[dataset],
        pathDict=path_dict,
    )
    predictor.loadData(inputDataList=[str(input_fasta)])
    predictor.featureEncode()
    predictor.doPredict()

    binary_report = output_dir / "binary_vector.csv"
    probability_report = output_dir / "probability_vector.csv"
    if not binary_report.exists() or not probability_report.exists():
        raise SystemExit(
            f"HEPAD did not produce expected outputs in {output_dir}"
        )

    # Cleanup massive intermediate files to save disk space
    import shutil
    mldata_dir = workdir / "mlData"
    if mldata_dir.exists():
        shutil.rmtree(mldata_dir, ignore_errors=True)

if __name__ == "__main__":
    main()
