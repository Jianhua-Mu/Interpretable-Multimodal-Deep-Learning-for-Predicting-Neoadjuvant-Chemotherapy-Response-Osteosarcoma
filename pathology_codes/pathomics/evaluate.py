"""Command-line evaluation workflow for pathology-only CLAM and MIL checkpoints."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pathomics.data.dataset_generic import Generic_MIL_Dataset
from pathomics.support.eval_utils import evaluate_dataset
from pathomics.train import parse_label_map


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate pathology classification checkpoints.")
    parser.add_argument("--dataset-csv", type=Path, required=True, help="CSV with case_id, slide_id, and label columns")
    parser.add_argument("--feature-root", type=Path, required=True, help="Folder containing pt_files/<slide_id>.pt")
    parser.add_argument("--model-directory", type=Path, required=True, help="Directory containing s_<fold>_checkpoint.pt files")
    parser.add_argument("--split-directory", type=Path, default=None, help="Directory containing split CSVs; defaults to model-directory")
    parser.add_argument("--output-directory", type=Path, default=Path("evaluation"))
    parser.add_argument("--evaluation-name", default="pathology_evaluation")
    parser.add_argument("--label-map", type=parse_label_map, default={0: 0, 1: 1}, help='JSON label map, e.g. "{\\"negative\\": 0, \\"positive\\": 1}"')
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--model-type", choices=["clam_sb", "clam_mb", "mil"], default="clam_sb")
    parser.add_argument("--model-size", choices=["small", "big"], default="small")
    parser.add_argument("--dropout", type=float, default=0.25)
    parser.add_argument("--fold-count", type=int, default=10)
    parser.add_argument("--fold-start", type=int, default=0)
    parser.add_argument("--fold-end", type=int, default=None)
    parser.add_argument("--fold", type=int, default=None, help="Evaluate only this fold")
    parser.add_argument("--split", choices=["train", "val", "test", "all"], default="test")
    parser.add_argument("--micro-average", action="store_true")
    return parser


def build_dataset(arguments: argparse.Namespace) -> Generic_MIL_Dataset:
    return Generic_MIL_Dataset(
        csv_path=str(arguments.dataset_csv),
        data_dir=str(arguments.feature_root),
        shuffle=False,
        print_info=True,
        label_dict=arguments.label_map,
        patient_strat=False,
        label_col=arguments.label_column,
        ignore=[],
    )


def evaluation_folds(arguments: argparse.Namespace) -> list[int]:
    if arguments.fold is not None:
        return [arguments.fold]
    final_fold = arguments.fold_end if arguments.fold_end is not None else arguments.fold_count
    if arguments.fold_start < 0 or final_fold <= arguments.fold_start or final_fold > arguments.fold_count:
        raise ValueError("fold range must satisfy 0 <= fold-start < fold-end <= fold-count")
    return list(range(arguments.fold_start, final_fold))


def main() -> None:
    arguments = build_parser().parse_args()
    split_directory = arguments.split_directory or arguments.model_directory
    required_paths = [arguments.dataset_csv, arguments.feature_root, arguments.model_directory, split_directory]
    missing_paths = [path for path in required_paths if not path.exists()]
    if missing_paths:
        raise FileNotFoundError(f"Required input path does not exist: {missing_paths[0]}")

    arguments.n_classes = len(set(arguments.label_map.values()))
    arguments.drop_out = arguments.dropout
    evaluation_directory = arguments.output_directory / arguments.evaluation_name
    evaluation_directory.mkdir(parents=True, exist_ok=True)
    (evaluation_directory / "settings.json").write_text(
        json.dumps({key: str(value) if isinstance(value, Path) else value for key, value in vars(arguments).items()}, indent=2),
        encoding="utf-8",
    )

    pathology_dataset = build_dataset(arguments)
    split_column = {"train": 0, "val": 1, "test": 2}
    fold_metrics: list[dict[str, float | int]] = []
    for fold_index in evaluation_folds(arguments):
        checkpoint_path = arguments.model_directory / f"s_{fold_index}_checkpoint.pt"
        if not checkpoint_path.is_file():
            raise FileNotFoundError(f"Missing checkpoint: {checkpoint_path}")

        if arguments.split == "all":
            evaluation_dataset = pathology_dataset
        else:
            split_csv = split_directory / f"splits_{fold_index}.csv"
            if not split_csv.is_file():
                raise FileNotFoundError(f"Missing split file: {split_csv}")
            evaluation_dataset = pathology_dataset.return_splits(from_id=False, csv_path=str(split_csv))[split_column[arguments.split]]

        _, _, error_rate, auc_score, prediction_table = evaluate_dataset(evaluation_dataset, arguments, str(checkpoint_path))
        prediction_table.to_csv(evaluation_directory / f"fold_{fold_index}_predictions.csv", index=False)
        fold_metrics.append({"fold": fold_index, "auc": auc_score, "accuracy": 1.0 - error_rate})

    pd.DataFrame(fold_metrics).to_csv(evaluation_directory / "evaluation_summary.csv", index=False)
    print(f"Evaluation complete. Results written to {evaluation_directory}")


if __name__ == "__main__":
    main()
