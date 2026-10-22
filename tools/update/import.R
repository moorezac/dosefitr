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
#'    \item `well`: Well identifier extracted from filename (e.g., "B10").
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
#'     basename prefix (before the suffix).
#' }
#' @export
import_trackmate <- function(
  spot_path,
  datetime_path = NA,
  spot_suffix = "_merge_spots.csv",
  datetime_suffix = "_datetime.csv"
) {
  # message(paste("Processing spot data from:\n", spot_path))

  # Get spot files
  spot_paths <- list.files(spot_path, pattern = spot_suffix, full.names = TRUE)

  if (length(spot_paths) == 0) {
    stop("No spot files found in: ", spot_path)
  }

  # Get datetime files if path provided
  if (!is.na(datetime_path)) {
    # message(paste0("Processing datetime data from:\n", datetime_path))
    time_paths <- list.files(
      datetime_path,
      pattern = datetime_suffix,
      full.names = TRUE
    )

    # if (length(spot_paths) != length(time_paths)) {
    # message("Warning: No. of spot files does not equal no. of time files")
    # }

    # Match files by basename
    spot_basenames <- sub(spot_suffix, "", basename(spot_paths))
    time_basenames <- sub(datetime_suffix, "", basename(time_paths))

    # Sort to match
    sort_order <- match(spot_basenames, time_basenames)
    spot_paths <- spot_paths[sort_order]
    time_paths <- time_paths[sort_order]

    # Remove NA matches
    spot_paths <- spot_paths[!is.na(spot_paths)]
  }

  # Helper function to convert character columns to numeric where possible
  convert_cols <- function(df) {
    for (col in names(df)) {
      # Try converting to numeric
      suppressWarnings({
        converted <- as.numeric(df[[col]])
      })
      # If conversion successful (not all NA), use it
      if (!all(is.na(converted)) || all(is.na(df[[col]]))) {
        df[[col]] <- converted
      }
    }
    df
  }

  common_prefix <- function(x) {
    s <- strsplit(x, "")
    n <- min(lengths(s))
    i <- 1
    while (i <= n && length(unique(sapply(s, `[`, i))) == 1) {
      i <- i + 1
    }
    return(paste(s[[1]][1:(i - 1)], collapse = ""))
  }

  common <- common_prefix(spot_paths)

  # Read spot files
  spot_list <- lapply(spot_paths, function(x) {
    # Extract well ID (e.g., "B10" from "VID207_B10_1-spots.csv")
    sub_x <- substring(x, nchar(common))
    m <- gregexpr(
      "(?<![A-Z0-9])[A-Z]{1,2}\\d{1,3}(?![A-Z0-9])",
      sub_x,
      perl = TRUE
    )
    well <- if (m[1] != -1) {
      toupper(regmatches(sub_x, m))
    } else {
      NA
    }

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
    if (!is.na(well)) {
      df <- cbind(well = well, df, stringsAsFactors = FALSE)
    }

    to_drop <- c(
      "TRACK_VISIBLE",
      "POSITION_T",
      "VISIBILITY",
      "MANUAL_SPOT_COLOR"
    )
    df <- df[, !names(df) %in% to_drop]

    return(df)
  })

  # Remove NA entries
  spot_list <- spot_list[
    !sapply(spot_list, function(x) length(x) == 1 && is.na(x))
  ]

  # Name list by wells
  names(spot_list) <- sapply(spot_list, function(x) unique(x$well)[1])

  # If datetime path provided, read and join
  if (!is.na(datetime_path)) {
    datetime_list <- lapply(time_paths, function(x) {
      df <- read.csv(x, header = FALSE, col.names = c("FRAME", "datetime"))

      df$datetime <- as.POSIXct(
        df$datetime,
        format = "%Y%m%d %H:%M:%OS",
        tz = "UTC"
      )
      # This converts to character
      # df$datetime <- format(df$datetime, "%Y-%m-%d %H:%M:%S")

      return(df)
    })

    # Extract well names from datetime files
    time_well_names <- sapply(basename(time_paths), function(x) {
      parts <- strsplit(x, "_")[[1]]
      if (length(parts) >= 2) parts[2] else NA
    })
    names(datetime_list) <- time_well_names

    # Filter datetime list to match spot data
    datetime_list <- datetime_list[names(datetime_list) %in% names(spot_list)]

    # Join datetime to spots
    spot_list <- mapply(
      function(spots, times) {
        # Merge by T column
        merged <- merge(spots, times, by = "FRAME", all.x = TRUE)

        # Convert T from 0-based to 1-based
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

      return(x)
    })
  }

  return(spot_list)
}


#' Summarise TrackMate Spot Data by Timepoint
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
