"""
Run Cellpose on all multi-series .tif images in a folder.

Usage:
    python incur_cellpose.py --input_dir --channels [-- --cellpose_args...]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import numpy as np
import tifffile
from tqdm import tqdm


def detect_gpu() -> bool:
    """
    Check for GPU
    """
    cuda_visible = os.getenv("CUDA_VISIBLE_DEVICES")
    if cuda_visible not in (None, "", "None"):
        return True

    try:
        result = subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def extract_frames_to_dir(
    tiff_path: Path,
    out_dir: Path,
    channels: Optional[List[int]] = None,
) -> List[Path]:
    """
    Extract individual frames from a time-series .tif file
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: List[Path] = []

    with tifffile.TiffFile(tiff_path) as tiff:
        series = tiff.series[0]
        axes = series.axes
        data = series.asarray()
        axis_map = {ax: i for i, ax in enumerate(axes)}

        # Frame/stack and channel dimensions
        frame_dim = axis_map.get("T", axis_map.get("Z", 0))
        channel_dim = axis_map.get("C")

        # Move frame dimension to the front
        if frame_dim != 0:
            data = np.moveaxis(data, frame_dim, 0)

        # Handle single-frame multi-channel case
        if data.ndim == 3 and "C" in axes and data.shape[channel_dim] < 10:
            data = data[np.newaxis, ...]

        # Filter channels if specified
        if channels is not None and channel_dim is not None:
            channels_idx = [c - 1 for c in channels]
            slicer = [slice(None)] * data.ndim
            slicer[channel_dim] = channels_idx
            data = data[tuple(slicer)]

        # Write each frame to disk
        for i, frame in enumerate(data):
            frame_path = out_dir / f"frame_{i:04d}.tif"
            tifffile.imwrite(frame_path, frame, imagej=True)
            frame_paths.append(frame_path)

    return frame_paths


def run_cellpose_dir(input_dir: Path, cellpose_args: list[str] | None = None) -> bool:
    """Run Cellpose CLI on all TIFFs in a folder with progress bar and logs"""
    use_gpu = detect_gpu()
    tiff_paths = sorted(input_dir.glob("*.tif")) + sorted(input_dir.glob("*.tiff"))
    if not tiff_paths:
        print(f"No .tif*(s) found in {input_dir}")
        return False

    command = ["cellpose", "--dir", str(input_dir), "--save_tif"]
    if use_gpu and "--use_gpu" not in command:
        command.append("--use_gpu")
    if "--verbose" not in command:
        command.append("--verbose")

    cellpose_args = cellpose_args or []
    command.extend(cellpose_args)

    log_path = input_dir / "cellpose.log"
    if log_path.exists():
        log_path.unlink()
    with open(log_path, "a") as log:
        log.write(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] {' '.join(command)}\n")

    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )

    total = len(tiff_paths)
    progress_bar = tqdm(total=total, unit="frame", leave=False)
    pattern = re.compile(r"(\d+)[/]+(\d+)")
    current = 0

    with open(log_path, "a") as log:
        for line in process.stdout:
            log.write(line)

            match = pattern.search(line)
            if not match:
                continue

            if match.group(0):
                new = current + 1
            else:
                continue

            if new > current and new <= total:
                progress_bar.update()
                current = new

    process.wait()
    progress_bar.close()

    ok = process.returncode == 0
    with open(log_path, "a") as log:
        log.write(f"[{'OK' if ok else 'ERROR'}] exit code {process.returncode}\n")

    if not ok:
        print(f"\nFailed in {input_dir}", file=sys.stderr)

    return ok


def stack_masks_to_timeseries(mask_files: List[Path], out_path: Path):
    """
    Stack multiple single-frame Cellpose mask TIFFs into one multi-frame .tif.
    """
    if not mask_files:
        return
    mask_files = sorted(mask_files)
    stack = np.stack([tifffile.imread(f) for f in mask_files], axis=0)
    tifffile.imwrite(out_path, stack, imagej=True)


def process_all(
    input_dir: Path,
    channels: Optional[List[int]],
    cellpose_args: Optional[List[str]] = None,
):
    """Extract all frames, run Cellpose, and clean up temporary frames."""
    tiff_paths = sorted(input_dir.glob("*.tif")) + sorted(input_dir.glob("*.tiff"))
    if not tiff_paths:
        print(f"No .tif files found in {input_dir}")
        return

    total = len(tiff_paths)
    print(f"\nFound {total} file(s) to process with Cellpose")

    progress_bar = tqdm(total=len(tiff_paths), unit="image", leave=False)

    for tiff_path in tiff_paths:

        progress_bar.set_description(f"{tiff_path.name}")

        with tempfile.TemporaryDirectory(dir=input_dir) as tmp_dir:
            tmp_path = Path(tmp_dir)

            extract_frames_to_dir(tiff_path, tmp_path, channels)

            run_cellpose_dir(tmp_path, cellpose_args)

            shutil.copy(
                tmp_path / "cellpose.log", input_dir / f"{tiff_path.stem}_cellpose.log"
            )

            mask_files = list(tmp_path.glob("*_cp_masks.tif"))
            if mask_files:
                stacked_path = input_dir / f"{tiff_path.stem}_cp_masks.tif"
                stack_masks_to_timeseries(mask_files, stacked_path)

        progress_bar.update()


def concat_logs_and_clean(input_dir: Path) -> None:
    """Concatenate all .log files in a folder, prefixing each with its filename"""
    input_dir = Path(input_dir)
    output_path = input_dir / "cellpose.log"
    
    log_files = sorted(p for p in input_dir.glob("*.log") if p.name != "cellpose.log")
    if not log_files:
        return
    
    with open(output_path, "w") as outfile:
        for log_file in log_files:
            outfile.write(f"\n{log_file.name}\n")
            with open(log_file, "r") as infile:
                outfile.write(infile.read().strip())
                outfile.write("\n")
    
    # Delete original log files
    for log_file in log_files:
        log_file.unlink()
    

def main():
    parser = argparse.ArgumentParser(
        description="Run Cellpose on all .tif files in a folder."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="Path to folder containing .tif files",
    )
    parser.add_argument(
        "--channels",
        nargs="*",
        type=int,
        default=1,
        help="Optional list of channels to include (1-based indexing)",
    )
    # Everything after -- is passed to Cellpose
    parser.add_argument(
        "cellpose_args",
        nargs=argparse.REMAINDER,
        help="Additional arguments to pass to Cellpose (use -- to separate from script arguments)",
    )

    args = parser.parse_args()

    # Remove leading '--'
    if args.cellpose_args and args.cellpose_args[0] == "--":
        args.cellpose_args = args.cellpose_args[1:]

    process_all(args.input_dir, args.channels, args.cellpose_args)
    concat_logs_and_clean(args.input_dir)

if __name__ == "__main__":
    main()
