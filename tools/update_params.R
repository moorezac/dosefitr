library(tidyverse)

template_path <- "tools/params_base.R"
template <- readLines(template_path)

# Helper: extract multi-line @param blocks from template
extract_template_params <- function(lines) {
  param_starts <- which(str_detect(lines, "^#'\\s*@param"))
  if (length(param_starts) == 0) {
    return(list())
  }

  params <- list()
  for (i in seq_along(param_starts)) {
    start_idx <- param_starts[i]
    end_idx <- start_idx

    # extend while next line is roxygen continuation (not a new tag)
    while (
      end_idx + 1 <= length(lines) &&
        str_detect(lines[end_idx + 1], "^#'") &&
        !str_detect(lines[end_idx + 1], "^#'\\s*@")
    ) {
      end_idx <- end_idx + 1
    }

    param_name <- str_match(lines[start_idx], "@param\\s+([^\\s]+)")[, 2]
    param_lines <- lines[start_idx:end_idx]

    # remove trailing empty roxygen lines from this param block
    while (
      length(param_lines) > 0 &&
        str_detect(param_lines[length(param_lines)], "^#'\\s*$")
    ) {
      param_lines <- param_lines[-length(param_lines)]
    }

    params[[param_name]] <- param_lines
  }
  params
}

template_params <- extract_template_params(template)

# helper: get next function line after an index
find_next_function_line <- function(lines, start_idx, look_ahead = 200) {
  n <- length(lines)
  end_search <- min(n, start_idx + look_ahead)
  hits <- which(str_detect(lines[(start_idx + 1):end_search], "<-\\s*function"))
  if (length(hits) == 0) {
    return(NA_integer_)
  }
  return(start_idx + hits[1])
}

# helper: remove empty roxygen lines (lines that are just "#'" or "#'  ")
remove_empty_roxy_lines <- function(lines) {
  empty_roxy <- str_detect(lines, "^#'\\s*$")
  lines[!empty_roxy]
}

update_file <- function(path) {
  lines <- readLines(path)
  n <- length(lines)
  if (n == 0) {
    return()
  }

  # first pass: remove empty roxygen lines throughout the file
  lines <- remove_empty_roxy_lines(lines)
  n <- length(lines)

  is_roxy <- str_detect(lines, "^#'")
  if (!any(is_roxy)) {
    return()
  }

  # find contiguous roxygen blocks
  rle_roxy <- rle(is_roxy)
  ends <- cumsum(rle_roxy$lengths)
  starts <- ends - rle_roxy$lengths + 1
  roxy_blocks <- data.frame(
    start = starts[rle_roxy$values],
    end = ends[rle_roxy$values]
  )

  for (b in seq_len(nrow(roxy_blocks))) {
    start <- roxy_blocks$start[b]
    end <- roxy_blocks$end[b]
    block <- lines[start:end]

    fun_line_idx <- find_next_function_line(lines, end)
    if (is.na(fun_line_idx)) {
      next
    }

    # extract formal args from function signature
    sig_lines <- lines[fun_line_idx]
    sig_i <- fun_line_idx
    while (
      !str_detect(sig_lines, "function\\(.*\\)") &&
        sig_i < min(n, fun_line_idx + 50)
    ) {
      sig_i <- sig_i + 1
      sig_lines <- paste0(sig_lines, collapse = " ")
      sig_lines <- paste0(sig_lines, " ", lines[sig_i])
    }

    arg_str <- str_match(sig_lines, "function\\((.*)\\)")[, 2]
    if (is.na(arg_str)) {
      arg_way <- str_match(sig_lines, "function\\((.*)")[, 2]
      if (is.na(arg_way) || arg_way == "") {
        next
      }
      arg_str <- arg_way
    }
    args <- trimws(unlist(str_split(arg_str, ",")))
    args <- args[args != ""]
    arg_names <- sub("=.*", "", args)
    arg_names <- trimws(arg_names)

    # find existing @param positions in block
    param_starts <- which(str_detect(block, "^#'\\s*@param\\b"))
    if (length(param_starts) > 0) {
      param_ends <- integer(length(param_starts))
      for (i in seq_along(param_starts)) {
        s_idx <- param_starts[i]
        e_idx <- s_idx
        while (
          e_idx + 1 <= length(block) &&
            str_detect(block[e_idx + 1], "^#'") &&
            !str_detect(block[e_idx + 1], "^#'\\s*@")
        ) {
          e_idx <- e_idx + 1
        }
        param_ends[i] <- e_idx
      }
      last_param_idx_in_block <- max(param_ends)
      insert_after_in_block <- last_param_idx_in_block
    } else {
      first_tag <- which(str_detect(block, "^#'\\s*@"))
      if (length(first_tag) > 0) {
        insert_after_in_block <- first_tag[1] - 1
      } else {
        insert_after_in_block <- length(block)
      }
    }

    # existing param names
    existing_param_lines <- block[str_detect(block, "^#'\\s*@param\\b")]
    existing_names <- if (length(existing_param_lines) > 0) {
      str_match(existing_param_lines, "@param\\s+([^\\s]+)")[, 2]
    } else {
      character(0)
    }

    missing <- setdiff(arg_names, existing_names)
    if (length(missing) == 0) {
      next
    }

    # build insert lines - handling multi-line templates
    insert_lines <- character(0)
    for (param_name in missing) {
      if (param_name %in% names(template_params)) {
        insert_lines <- c(insert_lines, template_params[[param_name]])
      } else {
        insert_lines <- c(
          insert_lines,
          paste0("#' @param ", param_name, " TODO: description.")
        )
      }
    }

    # rebuild block with insertion
    new_block <- append(block, insert_lines, after = insert_after_in_block)

    # write back into file lines
    lines <- c(
      if (start > 1) lines[1:(start - 1)] else character(0),
      new_block,
      if (end < length(lines)) lines[(end + 1):length(lines)] else character(0)
    )

    # recompute roxygen blocks for next iteration
    n <- length(lines)
    is_roxy <- str_detect(lines, "^#'")
    rle_roxy <- rle(is_roxy)
    ends <- cumsum(rle_roxy$lengths)
    starts <- ends - rle_roxy$lengths + 1
    roxy_blocks <- data.frame(
      start = starts[rle_roxy$values],
      end = ends[rle_roxy$values]
    )
  }

  path_new <- paste0("tools/update/", basename(path))
  writeLines(lines, path_new)
  message("Processed: ", basename(path))
}

files <- dir("R", full.names = TRUE, pattern = "\\.R$")
walk(files, update_file)
