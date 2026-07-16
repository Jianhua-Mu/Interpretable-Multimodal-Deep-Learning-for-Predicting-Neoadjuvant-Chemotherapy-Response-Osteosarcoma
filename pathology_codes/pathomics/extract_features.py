"""Extract ResNet features from pathology WSI patch bags."""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import h5py
import openslide
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from pathomics.data.dataset_h5 import Dataset_All_Bags, Whole_Slide_Bag_FP
from pathomics.models.resnet_custom import resnet50_baseline
from pathomics.support.file_utils import save_hdf5
from pathomics.support.utils import collate_features


COMPUTE_DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def extract_bag_features(
    patch_bag_path: Path,
    feature_h5_path: Path,
    slide_handle: openslide.OpenSlide,
    feature_encoder: nn.Module,
    batch_size: int,
    custom_downsample: int,
    target_patch_size: int,
) -> Path:
    patch_dataset = Whole_Slide_Bag_FP(
        file_path=str(patch_bag_path),
        wsi=slide_handle,
        pretrained=True,
        custom_downsample=custom_downsample,
        target_patch_size=target_patch_size,
    )
    loader_options = {"num_workers": 4, "pin_memory": True} if COMPUTE_DEVICE.type == "cuda" else {}
    patch_loader = DataLoader(patch_dataset, batch_size=batch_size, collate_fn=collate_features, **loader_options)
    print(f"Extracting {patch_bag_path.name}: {len(patch_loader)} batches")

    write_mode = "w"
    for batch_index, (image_batch, coordinates) in enumerate(patch_loader):
        if batch_index % 20 == 0:
            print(f"  batch {batch_index}/{len(patch_loader)}")
        with torch.no_grad():
            image_batch = image_batch.to(COMPUTE_DEVICE, non_blocking=True)
            feature_batch = feature_encoder(image_batch).cpu().numpy()
        save_hdf5(str(feature_h5_path), {"features": feature_batch, "coords": coordinates}, attr_dict=None, mode=write_mode)
        write_mode = "a"
    return feature_h5_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract pathology patch embeddings into H5 and PT files.")
    parser.add_argument("--patch-root", type=Path, required=True, help="Directory containing patches/<slide_id>.h5")
    parser.add_argument("--slide-directory", type=Path, required=True, help="Directory containing original WSI files")
    parser.add_argument("--process-list", type=Path, required=True, help="CSV listing the slide files to process")
    parser.add_argument("--feature-directory", type=Path, required=True, help="Destination for h5_files and pt_files")
    parser.add_argument("--slide-extension", default=".svs")
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--custom-downsample", type=int, default=1)
    parser.add_argument("--target-patch-size", type=int, default=-1)
    parser.add_argument("--overwrite", action="store_true", help="Recompute embeddings that already have a PT file")
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    required_paths = [arguments.patch_root, arguments.slide_directory, arguments.process_list]
    missing_paths = [path for path in required_paths if not path.exists()]
    if missing_paths:
        raise FileNotFoundError(f"Required input path does not exist: {missing_paths[0]}")

    h5_feature_directory = arguments.feature_directory / "h5_files"
    pt_feature_directory = arguments.feature_directory / "pt_files"
    h5_feature_directory.mkdir(parents=True, exist_ok=True)
    pt_feature_directory.mkdir(parents=True, exist_ok=True)

    feature_encoder = resnet50_baseline(pretrained=True).to(COMPUTE_DEVICE)
    if torch.cuda.device_count() > 1:
        feature_encoder = nn.DataParallel(feature_encoder)
    feature_encoder.eval()

    slide_bags = Dataset_All_Bags(str(arguments.process_list))
    for bag_index in range(len(slide_bags)):
        listed_slide = str(slide_bags[bag_index])
        slide_filename = os.path.basename(listed_slide)
        slide_id = Path(slide_filename).stem
        output_pt_path = pt_feature_directory / f"{slide_id}.pt"
        if output_pt_path.is_file() and not arguments.overwrite:
            print(f"Skipping existing feature file: {output_pt_path.name}")
            continue

        patch_bag_path = arguments.patch_root / "patches" / f"{slide_id}.h5"
        slide_path = arguments.slide_directory / f"{slide_id}{arguments.slide_extension}"
        if not patch_bag_path.is_file():
            raise FileNotFoundError(f"Patch bag not found: {patch_bag_path}")
        if not slide_path.is_file():
            raise FileNotFoundError(f"Slide not found: {slide_path}")

        feature_h5_path = h5_feature_directory / f"{slide_id}.h5"
        print(f"[{bag_index + 1}/{len(slide_bags)}] {slide_id}")
        started_at = time.time()
        with openslide.open_slide(str(slide_path)) as slide_handle:
            extract_bag_features(
                patch_bag_path,
                feature_h5_path,
                slide_handle,
                feature_encoder,
                arguments.batch_size,
                arguments.custom_downsample,
                arguments.target_patch_size,
            )

        with h5py.File(feature_h5_path, "r") as feature_file:
            slide_features = torch.from_numpy(feature_file["features"][:])
        torch.save(slide_features, output_pt_path)
        print(f"Finished {slide_id} in {time.time() - started_at:.1f}s")


if __name__ == "__main__":
    main()
