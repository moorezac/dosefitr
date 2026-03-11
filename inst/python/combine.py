import os
import sys
import xml.etree.ElementTree as ET
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Union

import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm


def check_folders_contain_same_names(directory_paths):
    """
    Check if all folders contain exactly the same set of names (files and subfolders).

    Args:
    folder_paths (list): A list of folder paths to compare

    Returns:
    tuple: (boolean indicating if all folders have the same names, set of common names)
    """
    if not directory_paths:
        return False, set()

    # Get the list of names (files and subfolders) for the first folder
    first_directory_names = set(os.listdir(directory_paths[0]))

    # Compare names with each subsequent folder
    for directory in directory_paths[1:]:
        # Get names for the current folder
        current_directory_names = set(os.listdir(directory))

        # If the sets of names are different, return False
        if current_directory_names != first_directory_names:
            return False, first_directory_names

    # If we've made it through all folders, they have the same names
    return True, first_directory_names


def extract_metadata_from_multipage_tiff(tiff):
    """
    Extract metadata from each page of a multi-page TIFF file
    """
    try:
        # List to store metadata for each page
        page_metadata_list = []

        # Iterate through each page in the TIFF
        for page_index, page in enumerate(tiff.pages):
            # Base metadata dictionary for this page
            metadata_dict = {
                "page_index": page_index,
                "image_width": page.shape[1] if len(page.shape) > 1 else None,
                "image_height": page.shape[0] if len(page.shape) > 0 else None,
                "dtype": str(page.dtype),
            }

            # Try to extract XML metadata from page tags
            for tag in page.tags.values():
                if tag.name == "ImageDescription":
                    try:
                        # Try to parse XML metadata
                        xml_str = tag.value
                        root = ET.fromstring(xml_str)

                        # Extract props from MetaData section
                        for prop in root.findall(".//prop"):
                            prop_name = prop.get("id", "Unknown")
                            prop_value = prop.get("value", "Unknown")
                            prop_type = prop.get("type", "string")

                            # Convert value based on type
                            if prop_type == "float":
                                try:
                                    prop_value = float(prop_value)
                                except (ValueError, TypeError):
                                    pass
                            elif prop_type == "int":
                                try:
                                    prop_value = int(prop_value)
                                except (ValueError, TypeError):
                                    pass
                            elif prop_type == "bool":
                                prop_value = prop_value.lower() == "true"

                            metadata_dict[prop_name] = prop_value

                        # Extract custom props
                        for prop in root.findall(".//custom-prop"):
                            prop_name = prop.get("id", "Unknown")
                            prop_value = prop.get("value", "Unknown")
                            metadata_dict[prop_name] = prop_value

                    except ET.ParseError:
                        # If XML parsing fails, store raw description
                        metadata_dict["image_description"] = xml_str
                    except Exception as e:
                        metadata_dict["metadata_extraction_error"] = str(e)

            # Add timing information if possible
            if hasattr(page, "time"):
                metadata_dict["page_time"] = page.time

            # Try to extract acquisition time from metadata
            if "acquisition-time-local" in metadata_dict:
                metadata_dict["acquisition_time"] = metadata_dict[
                    "acquisition-time-local"
                ]

            page_metadata_list.append(metadata_dict)

        return page_metadata_list

    except Exception as e:
        return [{"error": str(e)}]


def convert_16bit_to_8bit(
    input_data: Union[str, np.ndarray, List[np.ndarray]],
    output_path: Union[str, List[str]] = None,
) -> Union[np.ndarray, List[np.ndarray]]:
    """
    Convert 16-bit image(s) to 8-bit.

    Args:
        input_data: Can be:
            - A string path to a single TIFF file
            - A numpy array representing a single image
            - A list/series of numpy arrays
        output_path: Optional. Can be:
            - A string path for single image output
            - A list of paths for multiple images
            - None to skip saving to disk

    Returns:
        Converted 8-bit image(s) as numpy array(s)
    """

    def process_single_image(img_array: np.ndarray) -> np.ndarray:

        min_val = np.min(img_array)
        max_val = np.max(img_array)

        # avoid division by zero
        if max_val == min_val:
            return np.zeros_like(img_array, dtype=np.uint8)

        img_8bit = ((img_array - min_val) / (max_val - min_val) * 255).astype(np.uint8)
        return img_8bit

    # handle different input types
    if isinstance(input_data, str):
        # Single image path - use tifffile instead of PIL
        img_16bit = tifffile.imread(input_data)
        result = process_single_image(img_16bit)

        if output_path:
            tifffile.imwrite(output_path, result)

        return result

    elif isinstance(input_data, np.ndarray):
        # Single numpy array
        result = process_single_image(input_data)

        if output_path:
            tifffile.imwrite(output_path, result)

        return result

    elif isinstance(input_data, (list, tuple, np.ndarray)) and len(input_data) > 0:
        # list/series of numpy arrays
        results = []

        for i, img_array in enumerate(input_data):
            result = process_single_image(img_array)
            results.append(result)

            if output_path and isinstance(output_path, (list, tuple)):
                if i < len(output_path):
                    tifffile.imwrite(output_path[i], result)

        return results

    else:
        raise ValueError("invalid input type")


def process_individual_well(file_name, directory_list, output_directory):
    """
    For each well, extract metadata from each channel (as multi-page TIFF file), filter on common times, merge, and save

    Args:
    file_name: Name of the file across all directories to merge. Must be consistent across directories
    directory_list (list): A list of folder paths to compare
    output_directory (str): Directory to save output files
    """
    try:
        # Read in metadata information
        tiff_list = [
            tifffile.TiffFile(os.path.join(directory, file_name)) for directory in directory_list
        ]

        # Extract metadata for all images across all channels
        metadata_list = [
            extract_metadata_from_multipage_tiff(tiff=tiff) for tiff in tiff_list
        ]

        # Extract the datetime
        datetime_list = [
            pd.DataFrame(metadata)["acquisition-time-local"]
            for metadata in metadata_list
        ]

        # Find the common datetimes
        common_datetime = set(datetime_list[0]).intersection(
            *[set(dates) for dates in datetime_list[1:]]
        )

        # Find which indexes to keep for each channel
        index_list = [
            datetime[datetime.isin(common_datetime)].index for datetime in datetime_list
        ]

        # Read in the actual data
        tiff_np_list = [
            tifffile.imread(os.path.join(directory, file_name))
            for directory in directory_list
        ]

        # Filter down to the common datetime
        tiff_np_list_filtered = [
            np.take(a=tiff, indices=index, axis=0)
            for tiff, index in zip(tiff_np_list, index_list)
        ]

        # Make 8-bit
        tiff_np_list_filtered = [
            convert_16bit_to_8bit(tiff_np) for tiff_np in tiff_np_list_filtered
        ]

        # Merge them into a multi-channel time series
        time_series = np.stack(tiff_np_list_filtered, axis=1)

        # Save the multi-channel time series as a TIFF
        output_path = os.path.join(output_directory, file_name)
        tifffile.imwrite(
            output_path, time_series, metadata={"axes": "TCYX"}, imagej=True
        )

        # Save both full metadata and common datetime
        file_name_sans_ext = os.path.splitext(file_name)[0]

        metadata_df = [pd.DataFrame(metadata) for metadata in metadata_list]
        metadata_df = pd.concat(metadata_df)
        metadata_df.to_csv(
            os.path.join(output_directory, file_name_sans_ext + "_metadata.csv")
        )

        common_datetime_df = pd.DataFrame({"datetime": sorted(common_datetime)})
        common_datetime_df.to_csv(
            os.path.join(output_directory, file_name_sans_ext + "_datetime.csv"),
            header=False,
        )

        # print(f"Processed: {file_name}")
        return True

    except Exception as e:
        print(f"Error processing {file_name}: {str(e)}")
        return False


def process_files_parallel(file_names, directory_list, output_directory, n_cores=None):
    """
    Process multiple files in parallel using multiprocessing

    Args:
        file_names: List of file names to process
        directory_list: List of input directories
        output_directory: Output directory path
        n_cores: Number of cores to use (default: all available cores)
    """
    if n_cores is None:
        n_cores = cpu_count()

    print(f"Processing {len(file_names)} files using {n_cores} cores:")

    # Create a partial function with fixed arguments
    process_func = partial(
        process_individual_well,
        directory_list=directory_list,
        output_directory=output_directory,
    )

    # Use multiprocessing Pool with tqdm progress bar
    with Pool(processes=n_cores) as pool:
        results = list(
            tqdm(
                pool.imap(process_func, file_names),
                total=len(file_names),
                # desc="Processing:",
                unit=" file",
            )
        )

    # Summary
    successful = sum(results)
    failed = len(results) - successful
    print(f"\nProcessing complete: {successful} successful, {failed} failed")

    return results


if __name__ == "__main__":
    # Set multiprocessing start method to avoid cleanup issues
    import multiprocessing

    multiprocessing.set_start_method("spawn", force=True)

    if len(sys.argv) < 3:
        print(
            "Usage: python script.py <input_dir1> <input_dir2> ... <output_dir> [--cores N]"
        )
        print("  --cores N: Number of CPU cores to use (default: all available)")
        print(
            "  Need at least 2 arguments: one input directory and one output directory"
        )
        sys.exit(1)

    # Parse arguments
    args = sys.argv[1:]
    n_cores = None

    # Check for --cores argument
    if "--cores" in args:
        cores_idx = args.index("--cores")
        if cores_idx + 1 < len(args):
            try:
                n_cores = int(args[cores_idx + 1])
                # Remove --cores and the number from args
                args = args[:cores_idx] + args[cores_idx + 2 :]
            except ValueError:
                print("Error: --cores must be followed by a number")
                sys.exit(1)
        else:
            print("Error: --cores must be followed by a number")
            sys.exit(1)

    # Now check if we have enough arguments after removing --cores
    if len(args) < 2:
        print("Error: Need at least one input directory and one output directory")
        sys.exit(1)

    directory_list = args[:-1]
    output_directory = args[-1]

    print(f"Input: {directory_list}")
    print(f"Output: {output_directory}")
    # if n_cores:
        # print(f"Using {n_cores} cores")

    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Verify all input directories exist
    for directory in directory_list:
        if not os.path.exists(directory):
            print(f"Error: Input directory does not exist: {directory}")
            sys.exit(1)
        if not os.path.isdir(directory):
            print(f"Error: Path is not a directory: {directory}")
            sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Check that directories all contain the same files
    all_equal, file_names = check_folders_contain_same_names(
        directory_paths=directory_list
    )

    if not all_equal:
        print("Error: Mismatch in names of files between directories")
        print("Make sure all input directories contain the same set of files")
        sys.exit(1)

    if not file_names:
        print("Error: No files found in directories")
        sys.exit(1)

    # Convert set to list and sort for consistent processing order
    file_list = sorted(list(file_names))

    # Filter for TIFF files only (you can remove this filter if you want to process all files)
    tiff_files = [f for f in file_list if f.lower().endswith((".tif", ".tiff"))]

    if not tiff_files:
        print("Warning: No TIFF files found to process")
        print(f"Available files: {file_list[:10]}...")  # Show first 10 files
        # If you want to process all files regardless of extension, use:
        # tiff_files = file_list
        sys.exit(1)

    # print(f"Found {len(tiff_files)} files to process")

    # Process files in parallel
    process_files_parallel(tiff_files, directory_list, output_directory, n_cores)
