#' Import TrackMate Spot Data
#' Reads and processes TrackMate spot CSV files, optionally joining datetime
#' information from companion files. Returns a list of data frames organised
#' by well.
#' @param spot_path Character string specifying the directory containing
#'   TrackMate spot CSV files.
#' @param datetime_path Character string specifying the directory containing
#'   datetime CSV files. If NA (default), datetime information is not joined.
#' @param spot_suffix Character string specifying the filename suffix used to
#'   identify spot files. Default is \code{"_merge_spots.csv"}.
#' @param datetime_suffix Character string specifying the filename suffix used
#'   to identify datetime files. Default is \code{"_datetime.csv"}.
#' @return
#' A named list of data frames, one per well. Each data frame contains:
#'  \itemize{
#'    \item `datetime`: POSIXct timestamp (only if \code{datetime_path} provided).
#'    \item `well`: Well identifier extracted from filename (e.g., "B10"). If no
#'      well can be extracted, the filename with its suffix removed is used.
#'    \item `frame`: Frame number (1-indexed).
#'    \item `...`: Additional TrackMate spot statistics with cleaned column names.
#'  }
#' The following columns are dropped during import: \code{TRACK_VISIBLE},
#'  \code{POSITION_T}, \code{VISIBILITY}, \code{MANUAL_SPOT_COLOR}.
#' @note
#' \itemize{
#'   \item Empty spot files (no data rows) are silently excluded from the output.
#'   \item Frame numbers are converted from 0-indexed (TrackMate) to 1-indexed.
#'   \item Files are matched between spot and datetime directories by their
#'     basename prefix (before the suffix). Spot files without a matching
#'     datetime file are dropped when \code{datetime_path} is provided.
#'   \item If no well ID is found in a filename, the filename (minus suffix)
#'     is used as the identifier instead.
#' }
#' @export
import_trackmate <- function(
  spot_path,
  datetime_path = NA,
  spot_suffix = "_merge_spots.csv",
  datetime_suffix = "_datetime.csv"
) {
  # Helper: use extracted well ID, or fall back to filename minus suffix
  well_or_filename <- function(path, suffix, common_prefix_str) {
    w <- extract_well_id(path, common_prefix_str)
    if (length(w) == 0 || is.na(w) || !nzchar(w)) {
      w <- sub(suffix, "", basename(path), fixed = TRUE)
    }
    w
  }

  # Get spot files
  spot_paths <- list.files(spot_path, pattern = spot_suffix, full.names = TRUE)

  if (length(spot_paths) == 0) {
    stop("No spot files found in: ", spot_path)
  }

  # Get datetime files if path provided
  if (!is.na(datetime_path)) {
    time_paths <- list.files(
      datetime_path,
      pattern = datetime_suffix,
      full.names = TRUE
    )

    if (length(time_paths) == 0) {
      stop("No datetime files found in: ", datetime_path)
    }

    # Match files by basename
    spot_basenames <- sub(spot_suffix, "", basename(spot_paths), fixed = TRUE)
    time_basenames <- sub(
      datetime_suffix,
      "",
      basename(time_paths),
      fixed = TRUE
    )

    # Pair each spot file with its datetime file; drop unmatched spot files
    sort_order <- match(spot_basenames, time_basenames)
    keep <- !is.na(sort_order)
    spot_paths <- spot_paths[keep]
    time_paths <- time_paths[sort_order[keep]]

    if (length(spot_paths) == 0) {
      stop("No spot files have a matching datetime file.")
    }
  }

  common <- common_prefix(spot_paths)

  # Read spot files
  spot_list <- lapply(spot_paths, function(x) {
    # Extract well ID (e.g., "B10" from "VID207_B10_1-spots.csv"),
    # falling back to the filename if no well is found
    well <- well_or_filename(x, spot_suffix, common)

    # Read CSV
    df <- read.csv(x, row.names = NULL)
    df <- df[-c(1:3), ]

    # Return NA if empty
    if (nrow(df) == 0) {
      return(NA)
    }
    rownames(df) <- NULL

    # Convert columns
    df <- convert_cols(df)

    # Add well column at start
    df <- cbind(well = well, df, stringsAsFactors = FALSE)

    to_drop <- c(
      "TRACK_VISIBLE",
      "POSITION_T",
      "VISIBILITY",
      "MANUAL_SPOT_COLOR"
    )
    df <- df[, !names(df) %in% to_drop]

    return(df)
  })

  # Track which files were non-empty (needed to keep datetime files aligned)
  non_empty <- !vapply(
    spot_list,
    function(x) length(x) == 1 && is.na(x),
    logical(1)
  )
  spot_list <- spot_list[non_empty]

  if (length(spot_list) == 0) {
    stop("All spot files were empty.")
  }

  # Name list by wells
  names(spot_list) <- vapply(
    spot_list,
    function(x) as.character(unique(x$well)[1]),
    character(1)
  )

  # If datetime path provided, read and join
  if (!is.na(datetime_path)) {
    time_paths <- time_paths[non_empty]

    datetime_list <- lapply(time_paths, function(x) {
      df <- read.csv(x, header = FALSE, col.names = c("FRAME", "datetime"))

      df$datetime <- as.POSIXct(
        df$datetime,
        format = "%Y%m%d %H:%M:%OS",
        tz = "UTC"
      )

      return(df)
    })

    # Spot and datetime files are already paired positionally, so name the
    # datetime list to match the spot list directly
    names(datetime_list) <- names(spot_list)

    # Join datetime to spots
    spot_list <- mapply(
      function(spots, times) {
        # Merge by FRAME column
        merged <- merge(spots, times, by = "FRAME", all.x = TRUE)

        # Convert FRAME from 0-based to 1-based
        merged$FRAME <- merged$FRAME + 1

        # Reorder columns to put datetime first
        col_order <- c("datetime", setdiff(names(merged), "datetime"))
        merged <- merged[, col_order]

        return(merged)
      },
      spot_list,
      datetime_list,
      SIMPLIFY = FALSE
    )
  }

  # Clean column names
  spot_list <- lapply(spot_list, function(x) {
    colnames(x) <- tolower(colnames(x))
    colnames(x) <- gsub("[^a-z0-9]+", "_", colnames(x))
    colnames(x) <- gsub("^_+|_+$", "", colnames(x))
    colnames(x) <- gsub("_+", "_", colnames(x))

    return(x)
  })

  # Sort by datetime if available
  if (!is.na(datetime_path)) {
    spot_list <- lapply(spot_list, function(x) {
      x <- x[order(x$datetime), ]
      rownames(x) <- NULL

      return(x)
    })
  }

  return(spot_list)
}


#' Summarise TrackMate Spot Data by Timepoint
#'
#' @description
#' Aggregates spot-level TrackMate data to per-timepoint summaries, applying
#' user-specified functions to numeric columns and preserving non-numeric
#' columns.
#' @param spot_list A named list of data frames as returned by
#'   \code{\link{import_trackmate}}.
#' @param datetime_column Character string specifying the column to group by.
#'   Default is \code{"datetime"}.
#' @param multi_spot_fun A named list of summary functions to apply to numeric
#'   columns. Each function should accept a numeric vector and return a scalar.
#'   Default is \code{list(mean = \\(x) mean(x, na.rm = TRUE))}.
#' @param drop_cols Character vector of column names to exclude from the output.
#'   Default is \code{c("id", "label", "track_id", "track_name")}.
#' @param na.rm TODO: description.
#' @return
#' A named list of data frames (same structure as input), where each
#'  data frame contains one row per unique timepoint with:
#'   \itemize{
#'     \item \code{<datetime_column>}: The grouping timepoint.
#'     \item `num_spots`: Count of spots at that timepoint.
#'     \item \code{<fun>_<col>}: Summary statistics for each numeric column, named as
#'       `function_column`: (e.g., `mean_area`).
#'     \item \code{...}: Non-numeric columns preserved as first unique value per group.
#'   }
#' @seealso \code{\link{import_trackmate}} for importing raw TrackMate data
#' @examples
#' \dontrun{
#' summary_data <- summarise_trackmate(spot_list)
#' # Multiple summary functions
#' summary_data <- summarise_trackmate(
#'   spot_list,
#'   multi_spot_fun = list(
#'     mean = \(x) mean(x, na.rm = TRUE),
#'     sd = \(x) sd(x, na.rm = TRUE),
#'     median = \(x) median(x, na.rm = TRUE)
#'   )
#' )
#' # Custom grouping column and dropped columns
#' summary_data <- summarise_trackmate(
#'   spot_list,
#'   datetime_column = "frame",
#'   drop_cols = c("id", "label")
#' )
#' }
#' @export
summarise_trackmate <- function(
  spot_list,
  datetime_column = "datetime",
  multi_spot_fun = list(mean = \(x) mean(x, na.rm = TRUE)),
  drop_cols = c("id", "label", "track_id", "track_name")
) {
  result <- lapply(spot_list, function(df) {
    # Columns to summarise (excluding time and dropped)
    cols <- setdiff(
      names(df),
      c(
        datetime_column,
        if (any(is.na(drop_cols))) character() else drop_cols
      )
    )

    # Split data by time column
    split_data <- split(df, df[[datetime_column]])

    # Summarise per time group
    modified_list <- lapply(split_data, function(subdf) {
      out <- list()
      out[[datetime_column]] <- unique(subdf[[datetime_column]])
      out[["num_spots"]] <- nrow(subdf)

      for (col in cols) {
        col_data <- subdf[[col]]

        if (is.numeric(col_data)) {
          # Apply all numeric summary functions
          for (fname in names(multi_spot_fun)) {
            fun <- multi_spot_fun[[fname]]
            new_col <- paste(fname, col, sep = "_")
            out[[new_col]] <- fun(col_data)
          }
        } else {
          # Non-numeric: keep the first unique value, preserving its type
          val <- unique(col_data)[1]

          # If factor, keep factor value with same levels
          if (is.factor(col_data)) {
            out[[col]] <- factor(as.character(val), levels = levels(col_data))
          } else {
            # Date, POSIXct, character, etc. keep as-is (typed)
            out[[col]] <- val
          }
        }
      }

      out <- as.data.frame(out, stringsAsFactors = FALSE)

      return(out)
    })

    modified_df <- do.call(rbind, modified_list)
    rownames(modified_df) <- NULL

    return(modified_df)
  })

  return(result)
}


#' Longest common leading substring across a character vector
#'
#' @param x A character vector.
#' @return The shared prefix, or `""` when there is none.
#' @keywords internal
common_prefix <- function(x) {
  if (length(x) < 2) {
    return("")
  }

  chars <- strsplit(x, "")
  n <- min(lengths(chars))

  i <- 0L
  while (
    i < n && length(unique(vapply(chars, `[`, character(1), i + 1L))) == 1L
  ) {
    i <- i + 1L
  }

  substr(x[1], 1L, i)
}

#' Extract a well identifier from a filename
#' @keywords internal
extract_well_id <- function(path, common_prefix_str) {
  tail <- substring(path, nchar(common_prefix_str))
  m <- gregexpr(
    "(?<![A-Z0-9])[A-Z]{1,2}\\d{1,3}(?![A-Z0-9])",
    tail,
    perl = TRUE
  )
  if (m[[1]][1] == -1) {
    return(NA_character_)
  }
  toupper(regmatches(tail, m)[[1]][1])
}

#' Convert character columns to numeric where possible
#'
#' @param df A data frame.
#' @return `df` with each column converted to numeric where the conversion
#'   produces at least one non-`NA` value. Genuinely all-`NA` columns are left
#'   as they are.
#' @keywords internal
convert_cols <- function(df) {
  for (col in names(df)) {
    suppressWarnings(converted <- as.numeric(df[[col]]))
    if (!all(is.na(converted)) || all(is.na(df[[col]]))) {
      df[[col]] <- converted
    }
  }
  df
}
