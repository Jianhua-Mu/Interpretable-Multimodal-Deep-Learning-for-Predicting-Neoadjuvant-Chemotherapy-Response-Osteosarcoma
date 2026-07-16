"""Command-line training workflow for whole-slide pathology classification."""

from __future__ import annotations

import argparse
import json
import os
import random
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

from pathomics.data.dataset_generic import Generic_MIL_Dataset
from pathomics.support.core_utils import train_fold
from pathomics.support.file_utils import save_pkl


def _label_key(value: str) -> str | int | float:
    """Restore numeric JSON keys so they match labels read from a CSV file."""
    try:
        numeric_value = float(value)
    except ValueError:
        return value
    return int(numeric_value) if numeric_value.is_integer() else numeric_value


def parse_label_map(raw_mapping: str) -> dict[str | int | float, int]:
    try:
        parsed_mapping = json.loads(raw_mapping)
    except json.JSONDecodeError as exc:
        raise argparse.ArgumentTypeError("--label-map must be valid JSON") from exc

    if not isinstance(parsed_mapping, dict) or not parsed_mapping:
        raise argparse.ArgumentTypeError("--label-map must contain at least one label")

    normalized_mapping = {_label_key(str(label)): int(class_id) for label, class_id in parsed_mapping.items()}
    class_ids = sorted(set(normalized_mapping.values()))
    if class_ids != list(range(len(class_ids))):
        raise argparse.ArgumentTypeError("label-map class IDs must be consecutive integers starting at 0")
    return normalized_mapping


def configure_reproducibility(random_seed: int) -> None:
    random.seed(random_seed)
    os.environ["PYTHONHASHSEED"] = str(random_seed)
    np.random.seed(random_seed)
    torch.manual_seed(random_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(random_seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def build_dataset(arguments: argparse.Namespace) -> Generic_MIL_Dataset:
    return Generic_MIL_Dataset(
        csv_path=str(arguments.dataset_csv),
        data_dir=str(arguments.feature_root),
        shuffle=False,
        seed=arguments.seed,
        print_info=True,
        label_dict=arguments.label_map,
        patient_strat=False,
        label_col=arguments.label_column,
        ignore=[],
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train a pathology-only CLAM or MIL classifier.")
    parser.add_argument("--dataset-csv", type=Path, required=True, help="CSV with case_id, slide_id, and label columns")
    parser.add_argument("--feature-root", type=Path, required=True, help="Folder containing pt_files/<slide_id>.pt")
    parser.add_argument("--split-directory", type=Path, required=True, help="Folder containing splits_<fold>.csv files")
    parser.add_argument("--output-directory", type=Path, default=Path("outputs"), help="Root directory for checkpoints and metrics")
    parser.add_argument("--experiment-name", default="pathology_mil", help="Name used for the run directory")
    parser.add_argument("--label-map", type=parse_label_map, default={0: 0, 1: 1}, help='JSON label map, e.g. "{\\"negative\\": 0, \\"positive\\": 1}"')
    parser.add_argument("--label-column", default="label", help="Source label column in dataset-csv")
    parser.add_argument("--fold-count", type=int, default=10, help="Total number of cross-validation folds")
    parser.add_argument("--fold-start", type=int, default=0, help="First fold index to train")
    parser.add_argument("--fold-end", type=int, default=None, help="Exclusive final fold index; defaults to fold-count")
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=82)
    parser.add_argument("--optimizer", choices=["adam", "sgd"], default="adam")
    parser.add_argument("--model-type", choices=["clam_sb", "clam_mb", "mil"], default="clam_sb")
    parser.add_argument("--model-size", choices=["small", "big"], default="small")
    parser.add_argument("--dropout", type=float, default=0.25, help="Dropout probability; set 0 to disable")
    parser.add_argument("--bag-loss", choices=["ce", "svm"], default="ce")
    parser.add_argument("--instance-loss", choices=["ce", "svm"], default="svm")
    parser.add_argument("--instance-samples", type=int, default=8, help="Positive/negative patches sampled by CLAM")
    parser.add_argument("--bag-loss-weight", type=float, default=0.7, help="Relative weight of the bag-level CLAM loss")
    parser.add_argument("--disable-instance-clustering", action="store_true")
    parser.add_argument("--subtyping", action="store_true")
    parser.add_argument("--use-weighted-sampling", action="store_true")
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--tensorboard", action="store_true", help="Write TensorBoard logs when tensorboardX is installed")
    parser.add_argument("--testing", action="store_true", help="Use a small random subset for a smoke test")
    return parser


def serialize_settings(arguments: argparse.Namespace) -> dict[str, Any]:
    return {name: str(value) if isinstance(value, Path) else value for name, value in vars(arguments).items()}


def main() -> None:
    arguments = build_parser().parse_args()
    if arguments.fold_count < 1:
        raise ValueError("--fold-count must be positive")
    if arguments.fold_start < 0:
        raise ValueError("--fold-start must be non-negative")

    final_fold = arguments.fold_end if arguments.fold_end is not None else arguments.fold_count
    if final_fold <= arguments.fold_start or final_fold > arguments.fold_count:
        raise ValueError("fold range must satisfy 0 <= fold-start < fold-end <= fold-count")

    missing_inputs = [path for path in (arguments.dataset_csv, arguments.feature_root, arguments.split_directory) if not path.exists()]
    if missing_inputs:
        raise FileNotFoundError(f"Required input path does not exist: {missing_inputs[0]}")

    configure_reproducibility(arguments.seed)
    arguments.n_classes = len(set(arguments.label_map.values()))
    arguments.lr = arguments.learning_rate
    arguments.reg = arguments.weight_decay
    arguments.opt = arguments.optimizer
    arguments.drop_out = arguments.dropout
    arguments.inst_loss = arguments.instance_loss
    arguments.B = arguments.instance_samples
    arguments.bag_weight = arguments.bag_loss_weight
    arguments.no_inst_cluster = arguments.disable_instance_clustering
    arguments.weighted_sample = arguments.use_weighted_sampling
    arguments.log_data = arguments.tensorboard

    run_directory = arguments.output_directory / f"{arguments.experiment_name}_seed{arguments.seed}"
    run_directory.mkdir(parents=True, exist_ok=True)
    arguments.results_dir = str(run_directory)
    (run_directory / "settings.json").write_text(json.dumps(serialize_settings(arguments), indent=2), encoding="utf-8")

    pathology_dataset = build_dataset(arguments)
    fold_summaries: list[dict[str, float | int]] = []
    for fold_index in range(arguments.fold_start, final_fold):
        split_csv = arguments.split_directory / f"splits_{fold_index}.csv"
        if not split_csv.is_file():
            raise FileNotFoundError(f"Missing split file: {split_csv}")

        configure_reproducibility(arguments.seed)
        fold_datasets = pathology_dataset.return_splits(from_id=False, csv_path=str(split_csv))
        patient_predictions, test_auc, validation_auc, test_accuracy, validation_accuracy = train_fold(
            fold_datasets, fold_index, arguments
        )
        save_pkl(str(run_directory / f"fold_{fold_index}_predictions.pkl"), patient_predictions)
        fold_summaries.append(
            {
                "fold": fold_index,
                "test_auc": test_auc,
                "validation_auc": validation_auc,
                "test_accuracy": test_accuracy,
                "validation_accuracy": validation_accuracy,
            }
        )

    pd.DataFrame(fold_summaries).to_csv(run_directory / "cross_validation_summary.csv", index=False)
    print(f"Training complete. Results written to {run_directory}")


if __name__ == "__main__":
    main()
