#' Run Cellpose segmentation on a directory of images
#'
#' Calls the Cellpose segmentation pipeline on all \code{.tif} files in the
#' specified directory. Each image is split into individual frames, segmented,
#' and the resulting masks are stacked back into a time-series TIFF.
#'
#' @param input_dir Character. Path to the directory containing \code{.tif}
#'   files to segment.
#' @param channels Integer or integer vector. Channel(s) to use for
#'   segmentation (1-based indexing). Defaults to \code{1}.
#' @param cellpose_args Character vector of additional arguments passed directly
#'   to the Cellpose CLI. These are appended after \code{--} and passed
#'   through verbatim. Defaults to \code{NULL}.
#' @param python Character. Path to the Python executable. Defaults to
#'   \code{"python3"}.
#'
#' @return A character vector of stdout/stderr output from the Python process,
#'   returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' run_cellpose("/data/experiment1")
#' }
#'
#' @export
run_cellpose <- function(
  input_dir,
  channels = 1,
  cellpose_args = NULL,
  python = "python3"
) {
  script <- system.file("python", "cellpose.py", package = "incur")

  args <- c(
    script,
    "--input_dir",
    input_dir,
    "--channels",
    channels
  )

  if (!is.null(cellpose_args)) {
    args <- c(args, "--", cellpose_args)
  }

  system2(python, args = args, stdout = TRUE, stderr = TRUE)
}


#' Merge multi-channel time-series images
#'
#' Combines corresponding \code{.tif} files from multiple input directories
#' (one per channel) into a single multi-channel time-series TIFF. Frames are
#' aligned by acquisition timestamp, converted to 8-bit, and stacked along the
#' channel axis.
#'
#' @param input_dirs Character vector of paths to input directories, one per
#'   channel. All directories must contain identically named \code{.tif} files.
#' @param output_dir Character. Path to the directory where merged \code{.tif}
#'   files will be saved.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to \code{NULL} (all available cores).
#' @param python Character. Path to the Python executable. Defaults to
#'   \code{"python3"}.
#'
#' @return The stdout/stderr output from the Python process, returned
#'   invisibly. Called for its side effect of writing merged TIFFs to
#'   \code{output_dir}.
#'
#' @examples
#' \dontrun{
#' merge_channels(
#'   input_dirs = c("/data/channel1", "/data/channel2"),
#'   output_dir = "/data/merged"
#' )
#'
#' # Limit CPU usage
#' run_merge_channels(
#'   input_dirs = c("/data/ch1", "/data/ch2", "/data/ch3"),
#'   output_dir = "/data/merged",
#'   cores      = 4
#' )
#' }
#'
#' @export
run_merge_channels <- function(
  input_dirs,
  output_dir,
  cores = NULL,
  python = "python3"
) {
  script <- system.file("python", "merge_channels.py", package = "incur")

  args <- c(script, input_dirs, output_dir)

  if (!is.null(cores)) {
    args <- c(args, "--cores", cores)
  }

  result <- system2(python, args = args, stdout = TRUE, stderr = TRUE)

  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop("Channel merging failed:\n", paste(result, collapse = "\n"))
  }

  invisible(result)
}


#' Run TrackMate cell tracking on segmented images
#'
#' Merges Cellpose mask TIFFs with their corresponding images, generates a
#' TrackMate macro, and runs it via the Fiji/ImageJ CLI. Outputs per-image
#' \code{.xml} track files and \code{_spots.csv} tables to a \code{trackmate/}
#' subdirectory within \code{input_dir}.
#'
#' @param input_dir Character. Path to the directory containing \code{.tif}
#'   images and their corresponding \code{_cp_masks.tif} Cellpose mask files.
#' @param fiji_path Character. Path to the Fiji/ImageJ executable (e.g.
#'   \code{"/opt/fiji/Fiji.app/ImageJ-linux64"}).
#' @param verbose Logical. If \code{TRUE}, prints the TrackMate configuration
#'   and the Fiji command before execution. Defaults to \code{FALSE}.
#' @param target_channel Integer. Channel index (1-based) used for spot
#'   detection. Defaults to \code{1}.
#' @param simplify_contours Logical. Whether to simplify detected contours.
#'   Defaults to \code{TRUE}.
#' @param quality_threshold Numeric. Minimum spot quality score; spots below
#'   this value are filtered out. Defaults to \code{0.0}.
#' @param allow_frame_linking Logical. Whether to link spots between consecutive
#'   frames. Defaults to \code{TRUE}.
#' @param frame_linking_max_distance Numeric. Maximum distance (in pixels) for
#'   linking spots between frames. Defaults to \code{200.0}.
#' @param allow_gap_closing Logical. Whether to close gaps in tracks caused by
#'   missed detections. Defaults to \code{TRUE}.
#' @param gap_closing_max_distance Numeric. Maximum distance (in pixels) for
#'   gap closing. Defaults to \code{200.0}.
#' @param gap_closing_max_frame_gap Integer. Maximum number of frames over
#'   which a gap can be closed. Defaults to \code{2}.
#' @param allow_track_splitting Logical. Whether to allow tracks to split (e.g.
#'   cell division). Defaults to \code{TRUE}.
#' @param track_splitting_max_distance Numeric. Maximum distance (in pixels)
#'   for a track split event. Defaults to \code{15.0}.
#' @param allow_track_merging Logical. Whether to allow tracks to merge.
#'   Defaults to \code{FALSE}.
#' @param track_merging_max_distance Numeric. Maximum distance (in pixels) for
#'   a track merge event. Defaults to \code{15.0}.
#' @param python Character. Path to the Python executable. Defaults to
#'   \code{"python3"}.
#'
#' @return The stdout/stderr output from the Python process, returned
#'   invisibly. Called for its side effect of writing TrackMate \code{.xml}
#'   and \code{_spots.csv} files.
#'
#' @examples
#' \dontrun{
#' # Basic usage with defaults
#' run_trackmate(
#'   input_dir = "/data/experiment1",
#'   fiji_path = "/opt/fiji/Fiji.app/ImageJ-linux64"
#' )
#'
#' # Tighter linking distances, no splitting or merging
#' run_trackmate(
#'   input_dir                  = "/data/experiment1",
#'   fiji_path                  = "/opt/fiji/Fiji.app/ImageJ-linux64",
#'   frame_linking_max_distance = 100,
#'   gap_closing_max_frame_gap  = 3L,
#'   allow_track_splitting      = FALSE,
#'   allow_track_merging        = FALSE,
#'   verbose                    = TRUE
#' )
#' }
#'
#' @export
run_trackmate <- function(
  input_dir,
  fiji_path,
  verbose = FALSE,
  target_channel = 1,
  simplify_contours = TRUE,
  quality_threshold = 0.0,
  allow_frame_linking = TRUE,
  frame_linking_max_distance = 200.0,
  allow_gap_closing = TRUE,
  gap_closing_max_distance = 200.0,
  gap_closing_max_frame_gap = 2L,
  allow_track_splitting = TRUE,
  track_splitting_max_distance = 15.0,
  allow_track_merging = FALSE,
  track_merging_max_distance = 15.0,
  python = "python3"
) {
  script <- system.file("python", "trackmate.py", package = "incur")

  bool_arg <- function(flag, value) {
    if (isTRUE(value)) {
      paste0("--", flag)
    } else {
      paste0("--no-", gsub("_", "-", flag))
    }
  }

  args <- c(
    script,
    "--input_dir",
    input_dir,
    "--fiji_path",
    fiji_path,
    "--target_channel",
    target_channel,
    bool_arg("simplify_contours", simplify_contours),
    "--quality_threshold",
    quality_threshold,
    bool_arg("allow_frame_linking", allow_frame_linking),
    "--frame_linking_max_distance",
    frame_linking_max_distance,
    bool_arg("allow_gap_closing", allow_gap_closing),
    "--gap_closing_max_distance",
    gap_closing_max_distance,
    "--gap_closing_max_frame_gap",
    gap_closing_max_frame_gap,
    bool_arg("allow_track_splitting", allow_track_splitting),
    "--track_splitting_max_distance",
    track_splitting_max_distance,
    bool_arg("allow_track_merging", allow_track_merging),
    "--track_merging_max_distance",
    track_merging_max_distance
  )

  if (verbose) {
    args <- c(args, "--verbose")
  }

  result <- system2(python, args = args, stdout = TRUE, stderr = TRUE)

  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop("TrackMate failed:\n", paste(result, collapse = "\n"))
  }

  invisible(result)
}
