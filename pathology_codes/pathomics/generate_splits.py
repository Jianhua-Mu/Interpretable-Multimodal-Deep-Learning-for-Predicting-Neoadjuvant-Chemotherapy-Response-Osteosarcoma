"""Create reproducible patient-level train, validation, and test folds for pathology data."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pathomics.data.dataset_generic import Generic_WSI_Classification_Dataset, save_splits
from pathomics.train import parse_label_map


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate cross-validation splits for pathology slide data.")
    parser.add_argument("--dataset-csv", type=Path, required=True, help="CSV with case_id, slide_id, and label columns")
    parser.add_argument("--output-directory", type=Path, required=True, help="Destination for splits_<fold>.csv files")
    parser.add_argument("--label-map", type=parse_label_map, default={0: 0, 1: 1})
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--fold-count", type=int, default=10)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--label-fraction", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=82)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    if not arguments.dataset_csv.is_file():
        raise FileNotFoundError(f"Dataset CSV not found: {arguments.dataset_csv}")
    for fraction_name in ("validation_fraction", "test_fraction", "label_fraction"):
        fraction = getattr(arguments, fraction_name)
        if not 0 < fraction <= 1:
            raise ValueError(f"--{fraction_name.replace('_', '-')} must be in (0, 1]")
    if arguments.validation_fraction + arguments.test_fraction >= 1:
        raise ValueError("validation-fraction plus test-fraction must be below 1")

    pathology_dataset = Generic_WSI_Classification_Dataset(
        csv_path=str(arguments.dataset_csv),
        shuffle=False,
        seed=arguments.seed,
        print_info=True,
        label_dict=arguments.label_map,
        patient_strat=True,
        label_col=arguments.label_column,
        ignore=[],
    )
    class_patient_counts = np.array([len(patient_ids) for patient_ids in pathology_dataset.patient_cls_ids])
    validation_counts = np.floor(class_patient_counts * arguments.validation_fraction).astype(int)
    test_counts = np.floor(class_patient_counts * arguments.test_fraction).astype(int)
    if np.any(validation_counts == 0) or np.any(test_counts == 0):
        raise ValueError("At least one class is too small for the selected validation/test fractions")

    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    pathology_dataset.create_splits(
        k=arguments.fold_count,
        val_num=validation_counts,
        test_num=test_counts,
        label_frac=arguments.label_fraction,
    )
    for fold_index in range(arguments.fold_count):
        pathology_dataset.set_splits()
        split_datasets = pathology_dataset.return_splits(from_id=True)
        descriptor_table = pathology_dataset.test_split_gen(return_descriptor=True)
        save_splits(split_datasets, ["train", "val", "test"], arguments.output_directory / f"splits_{fold_index}.csv")
        save_splits(
            split_datasets,
            ["train", "val", "test"],
            arguments.output_directory / f"splits_{fold_index}_bool.csv",
            boolean_style=True,
        )
        descriptor_table.to_csv(arguments.output_directory / f"splits_{fold_index}_descriptor.csv")

    print(f"Created {arguments.fold_count} folds in {arguments.output_directory}")


if __name__ == "__main__":
    main()
