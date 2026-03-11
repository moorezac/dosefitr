"""
Run Trackmate on a folder of images with associated masks.

Usage:
    python incur_trackmate.py --input_dir --fiji_path [...]
"""

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import tifffile
from tqdm import tqdm


@dataclass
class TrackMateConfig:
    """Configuration for TrackMate tracking parameters."""

    # Detector
    target_channel: int = field(
        default=1,
        metadata={"help": "Target channel for detection", "group": "Detector"},
    )
    simplify_contours: bool = field(
        default=True, metadata={"help": "Simplify contours", "group": "Detector"}
    )

    # Spot
    quality_threshold: float = field(
        default=0.0,
        metadata={"help": "Minimum quality threshold for spots", "group": "Spot"},
    )

    # Tracker - Frame Linking
    allow_frame_linking: bool = field(
        default=True,
        metadata={"help": "Link spots between frames", "group": "Tracker Linking"},
    )
    frame_linking_max_distance: float = field(
        default=200.0,
        metadata={"help": "Max distance for frame linking", "group": "Tracker Linking"},
    )

    # Tracker - Gap Closing
    allow_gap_closing: bool = field(
        default=True,
        metadata={"help": "Allow gap closing", "group": "Tracker Gap Closing"},
    )
    gap_closing_max_distance: float = field(
        default=200.0,
        metadata={
            "help": "Max distance for gap closing",
            "group": "Tracker Gap Closing",
        },
    )
    gap_closing_max_frame_gap: int = field(
        default=2,
        metadata={"help": "Max frame gap allowed", "group": "Tracker Gap Closing"},
    )

    # Tracker - Splitting
    allow_track_splitting: bool = field(
        default=True,
        metadata={"help": "Allow track splitting", "group": "Tracker Splitting"},
    )
    track_splitting_max_distance: float = field(
        default=15.0,
        metadata={
            "help": "Max distance for track splitting",
            "group": "Tracker Splitting",
        },
    )

    # Tracker - Merging
    allow_track_merging: bool = field(
        default=False,
        metadata={"help": "Allow track merging", "group": "Tracker Merging"},
    )
    track_merging_max_distance: float = field(
        default=15.0,
        metadata={"help": "Max distance for track merging", "group": "Tracker Merging"},
    )

    @classmethod
    def add_arguments(cls, parser: argparse.ArgumentParser) -> None:
        """Add dataclass fields as arguments to argparse parser."""
        groups = {}
        for field_name, field_def in cls.__dataclass_fields__.items():
            group_name = field_def.metadata.get("group", "TrackMate Settings")
            groups.setdefault(group_name, parser.add_argument_group(group_name))
            arg_name = f"--{field_name}"
            default_value = field_def.default
            help_text = (
                f"{field_def.metadata.get('help','')} (default: {default_value})"
            )
            if field_def.type == bool:
                groups[group_name].add_argument(
                    arg_name, action="store_true", default=default_value, help=help_text
                )
                groups[group_name].add_argument(
                    f"--no-{field_name.replace('_','-')}",
                    dest=field_name,
                    action="store_false",
                    help=f"Disable {field_name.replace('_',' ')}",
                )
            else:
                groups[group_name].add_argument(
                    arg_name, type=field_def.type, default=default_value, help=help_text
                )

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "TrackMateConfig":
        return cls(
            **{
                f: getattr(args, f)
                for f in cls.__dataclass_fields__
                if hasattr(args, f)
            }
        )

    def __str__(self) -> str:
        grouped = {}
        for f, fd in self.__dataclass_fields__.items():
            group = fd.metadata.get("group", "Other")
            grouped.setdefault(group, []).append((f, getattr(self, f)))
        lines = ["TrackMate Configuration:"]
        for group, fields in grouped.items():
            lines.append(f"\n{group}:")
            lines.extend([f"  {f}: {v}" for f, v in fields])
        return "\n".join(lines)


def generate_trackmate_macro(config: TrackMateConfig, input_dir: str) -> str:
    """Generate a Python macro for TrackMate."""
    link_dist = config.frame_linking_max_distance if config.allow_frame_linking else 0.0
    gap_dist = config.gap_closing_max_distance if config.allow_gap_closing else 0.0
    max_gap = config.gap_closing_max_frame_gap if config.allow_gap_closing else 0
    split_dist = (
        config.track_splitting_max_distance if config.allow_track_splitting else 0.0
    )
    merge_dist = (
        config.track_merging_max_distance if config.allow_track_merging else 0.0
    )

    return f"""# @String input_dir
import os
import glob

from ij import IJ
from fiji.plugin.trackmate import Logger, Model, SelectionModel, Settings, TrackMate
from fiji.plugin.trackmate.detection import LabelImageDetectorFactory
from fiji.plugin.trackmate.features import FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.io import CSVExporter, TmXmlWriter
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from java.io import File


def run_trackmate(image_path):
    imp = IJ.openImage(image_path)
    if imp is None: return
    logger = Logger.IJ_LOGGER
    model = Model(); model.setLogger(Logger.IJ_LOGGER)
    settings = Settings(imp)
    settings.detectorFactory = LabelImageDetectorFactory()
    settings.detectorSettings["TARGET_CHANNEL"] = {config.target_channel}
    settings.detectorSettings["SIMPLIFY_CONTOURS"] = {str(config.simplify_contours)}
    settings.addSpotFilter(FeatureFilter('QUALITY', {config.quality_threshold}, True))
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings.update({{
        "LINKING_MAX_DISTANCE": {link_dist},
        "GAP_CLOSING_MAX_DISTANCE": {gap_dist},
        "MAX_FRAME_GAP": {max_gap},
        "ALLOW_TRACK_SPLITTING": {str(config.allow_track_splitting)},
        "SPLITTING_MAX_DISTANCE": {split_dist},
        "ALLOW_TRACK_MERGING": {str(config.allow_track_merging)},
        "MERGING_MAX_DISTANCE": {merge_dist}
    }})
    settings.addAllAnalyzers()
    trackmate = TrackMate(model, settings)
    if not trackmate.checkInput() or not trackmate.process(): imp.close(); return
    base_name = os.path.splitext(os.path.basename(image_path))[0]
    folder = os.path.dirname(image_path)
    ds = DisplaySettingsIO.readUserDefault()
    
    out_xml = File(os.path.join(folder, base_name+".xml"))
    writer = TmXmlWriter(out_xml)
    writer.appendSettings(settings)
    writer.appendModel(model)		
    writer.appendDisplaySettings(ds)	
    writer.writeToFile()
    CSVExporter.exportSpots(os.path.join(folder, base_name+"_spots.csv"), model, True)
    imp.close()


input_folder = "{input_dir}"
if os.path.exists(input_folder):
    for tiff_path in sorted(glob.glob(os.path.join(input_folder, "*merge.tif"))):
        name = os.path.splitext(os.path.basename(tiff_path))[0]
        print("NAME " + name)
        run_trackmate(tiff_path)
        print("SUCCESS")
else:
    print("Error: folder not found")
"""


def merge_mask_and_image(mask_path: Path, image_path: Path, output_dir: Path) -> Path:
    """Merge mask and image into a multi-channel TIFF."""
    mask = tifffile.imread(mask_path)
    image = tifffile.imread(image_path)

    if mask.ndim == 3:
        mask = np.expand_dims(mask, 1)
    if image.ndim == 3:
        image = np.expand_dims(image, 1)

    name = image_path.stem
    output_path = output_dir / f"{name}_merge.tif"

    merged = np.concatenate([mask, image], axis=1)
    tifffile.imwrite(output_path, merged, imagej=True, metadata={"axes": "TCYX"})
    return output_path


def batch_merge_masks_and_images(input_dir: Path) -> list[Path]:
    """Merge all mask files with their corresponding images into trackmate folder."""
    output_dir = input_dir / "trackmate"
    output_dir.mkdir(exist_ok=True, parents=True)

    mask_paths = sorted(input_dir.glob("*_cp_masks.tif"))
    if not mask_paths:
        print("No mask files found.")
        return

    total = len(mask_paths)
    print(f"\nFound {total} pairs of mask(s) to merge")

    merged_files = []
    progress_bar = tqdm(total=total, unit="image", leave=True)

    for mask_path in mask_paths:
        image_name = mask_path.name.replace("_cp_masks.tif", ".tif")
        image_path = input_dir / image_name

        if not image_path.exists():
            print(f"No matching image found for {mask_path.name}")
            continue

        progress_bar.set_description(f"{image_name}")
        try:
            merged_file = merge_mask_and_image(mask_path, image_path, output_dir)
            merged_files.append(merged_file)
        except Exception as e:
            print(f"Error merging {mask_path.name}: {e}")

        progress_bar.update()

    return merged_files


def process_with_trackmate(
    input_dir: Path, fiji_path: Path, config: TrackMateConfig, verbose: bool = False
) -> bool:
    """Process images with TrackMate using the given configuration."""
    merged_files = batch_merge_masks_and_images(input_dir)

    output_dir = input_dir / "trackmate"

    if not merged_files:
        print("No merged files to process.")
        return False

    macro_file = output_dir / "macro.py"
    macro_file.write_text(generate_trackmate_macro(config, output_dir))

    command = [
        fiji_path,
        "-Djava.awt.headless=true",
        "--run",
        str(macro_file),
        f"input_dir='{output_dir}'",
    ]
    if verbose:
        print("Running command:", " ".join(command))

    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )

    total = len(merged_files)
    print(f"\nFound {total} file(s) to process with Trackmate")

    progress_bar = tqdm(total=total, unit="image", leave=False)
    pattern_name = re.compile(r"NAME (.+)")
    pattern_success = re.compile(r"SUCCESS")
    current = 0

    log_path = output_dir / "trackmate.log"
    if log_path.exists():
        log_path.unlink()

    with log_path.open("a") as log:
        for line in process.stdout:
            log.write(line)

            # Update current file name
            match_name = pattern_name.match(line)
            if match_name:
                name = match_name.group(1)
                progress_bar.set_description_str(name)

            # Update progress
            match_success = pattern_success.match(line)
            if match_success:
                progress_bar.update()

    process.wait

    return process.returncode == 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="TrackMate pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing images and masks"
    )
    parser.add_argument(
        "--fiji_path", required=True, help="Path to Fiji/ImageJ executable"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    TrackMateConfig.add_arguments(parser)
    args = parser.parse_args()

    config = TrackMateConfig.from_args(args)

    if args.verbose:
        print(config)

    args.input_dir, args.fiji_path = Path(args.input_dir), Path(args.fiji_path)
    success = process_with_trackmate(
        args.input_dir, args.fiji_path, config, verbose=args.verbose
    )
    if not success:
        sys.exit("TrackMate processing failed.")


if __name__ == "__main__":
    main()
