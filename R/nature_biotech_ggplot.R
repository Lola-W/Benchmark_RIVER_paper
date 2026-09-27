# Nature Biotechnology ggplot styling helpers
#
# Palette: qualitative, colour-blind-aware (Okabe-Ito based), with one stable
# hue per family: Reference #4D4D4D; ppmSeq #0072B2; UDSeq #009E73;
# NanoSeq #E69F00; HiDEF-seq #D55E00; IlluminaSeq #CC79A7;
# PTA-Seq #56B4E9.
# Subtypes use family shades: lower coverage/singleton/single-stranded are
# lighter; multiread/double-stranded are darker. A lone subtype uses the base.
# Names are normalized before matching (illumina -> IlluminaSeq;
# ppmSeq_singleton_HC -> ppmSeq_singleton; leading Duplex_ is removed;
# existing_any_PTA1_16 -> PTA-Seq; nanoseq_raw -> NanoSeq (raw)).
#
# Source this file, then add the theme and method scales to a ggplot:
#
#   source("R/nature_biotech_ggplot.R")
#
#   dat$method <- normalize_seq_method(dat$method)
#   ggplot2::ggplot(dat, ggplot2::aes(x, y, colour = method)) +
#     ggplot2::geom_line() +
#     scale_colour_seq_method(dat$method) +
#     theme_nature_biotech()
#
# Or add everything in one call:
#
#   p + style_nature_biotech(dat$method, aesthetics = "colour")
#
# This is a standalone script. ggplot2 is the only non-base dependency.


# Stable, colour-blind-aware base colours for the sequencing technologies.
# These values are public so the same palette can also be used outside ggplot.
NBT_METHOD_BASE_COLOURS <- c(
  "ppmSeq"    = "#0072B2",
  "UDSeq"     = "#009E73",
  "NanoSeq"   = "#E69F00",
  "HiDEF-seq" = "#D55E00",
  "IlluminaSeq" = "#CC79A7",
  "PTA-Seq"     = "#56B4E9",
  "Reference" = "#4D4D4D"
)

# American-spelling alias.
NBT_METHOD_BASE_COLORS <- NBT_METHOD_BASE_COLOURS


.nbt_assert_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "The 'ggplot2' package is required. Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
}


.nbt_clean_method_key <- function(x) {
  x <- trimws(x)
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x)
  x <- gsub("[[:space:]-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}


.nbt_clean_suffix <- function(x) {
  if (!nzchar(x)) {
    return("")
  }

  # Collapse known spelling variants to a stable vocabulary.
  x <- sub("^\\(?raw\\)?$", "raw", x, perl = TRUE)
  x <- gsub("singletons?_hc", "singleton", x, perl = TRUE)
  x <- gsub("singletons?", "singleton", x, perl = TRUE)
  x <- gsub("multi_?reads?", "multiread", x, perl = TRUE)
  x <- gsub("double_?strand(?:ed)?", "double_stranded", x, perl = TRUE)
  x <- gsub("single_?strand(?:ed)?", "single_stranded", x, perl = TRUE)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  # Keep common technical abbreviations readable without altering sample IDs.
  parts <- strsplit(x, "_", fixed = TRUE)[[1L]]
  parts[parts == "hc"] <- "HC"
  paste(parts, collapse = "_")
}


#' Normalize sequencing-method labels
#'
#' Canonicalizes method-family case and punctuation while retaining meaningful
#' subtype/sample suffixes. In particular, ppmSeq_singleton_HC and plural
#' singleton variants become ppmSeq_singleton, every capitalization of illumina
#' becomes IlluminaSeq, raw variants use the display suffix " (raw)", leading
#' Duplex_ is removed, and existing_any_PTA1_16 becomes PTA-Seq.
#'
#' @param x A character or factor vector.
#' @return A character vector with names preserved.
normalize_seq_method <- function(x) {
  original_names <- names(x)
  input <- as.character(x)
  output <- input
  present <- !is.na(input)

  for (i in which(present)) {
    key <- .nbt_clean_method_key(input[[i]])
    key <- sub("^duplex_", "", key)
    if (key == "existing_any_pta1_16") key <- "pta"

    specs <- list(
      list(pattern = "^ppm_?seq(?:_|$)", strip = "^ppm_?seq_?", family = "ppmSeq"),
      list(pattern = "^ud_?seq(?:_|$)", strip = "^ud_?seq_?", family = "UDSeq"),
      list(pattern = "^nano_?seq(?:_|$)", strip = "^nano_?seq_?", family = "NanoSeq"),
      list(pattern = "^hidef(?:_?seq)?(?:_|$)", strip = "^hidef(?:_?seq)?_?", family = "HiDEF-seq"),
      list(pattern = "^illumina(?:_?seq)?(?:_|$)", strip = "^illumina(?:_?seq)?_?", family = "IlluminaSeq"),
      list(pattern = "^pta(?:_?seq)?(?:_|$)", strip = "^pta(?:_?seq)?_?", family = "PTA-Seq")
    )

    matched <- FALSE
    for (spec in specs) {
      if (grepl(spec$pattern, key, perl = TRUE)) {
        suffix <- sub(spec$strip, "", key, perl = TRUE)
        suffix <- .nbt_clean_suffix(suffix)
        output[[i]] <- if (identical(suffix, "raw")) {
          paste0(spec$family, " (raw)")
        } else if (nzchar(suffix)) {
          paste0(spec$family, "_", suffix)
        } else {
          spec$family
        }
        matched <- TRUE
        break
      }
    }

    if (!matched && key %in% c("ref", "reference", "ref_hg38", "hg38")) {
      output[[i]] <- "Reference"
      matched <- TRUE
    }

    if (!matched) {
      output[[i]] <- trimws(input[[i]])
    }
  }

  names(output) <- original_names
  output
}


#' Normalize one method column in a data frame
#'
#' @param data A data.frame or tibble.
#' @param column Character name of the column to normalize.
#' @return A copy of data with the selected column normalized.
normalize_seq_method_column <- function(data, column) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or tibble.", call. = FALSE)
  }
  if (!is.character(column) || length(column) != 1L || !column %in% names(data)) {
    stop("'column' must be the name of a column in 'data'.", call. = FALSE)
  }

  data[[column]] <- normalize_seq_method(data[[column]])
  data
}


#' Identify the technology family for normalized or raw labels
#'
#' @param x A character or factor vector of method labels.
#' @return A character vector. Unknown labels are returned as "Other".
seq_method_family <- function(x) {
  x <- normalize_seq_method(x)
  family <- rep("Other", length(x))
  family[is.na(x)] <- NA_character_

  known <- names(NBT_METHOD_BASE_COLOURS)
  for (candidate in known) {
    family[
      !is.na(x) & (
        x == candidate |
          startsWith(x, paste0(candidate, "_")) |
          startsWith(x, paste0(candidate, " ("))
      )
    ] <- candidate
  }
  family
}


.nbt_mix_colour <- function(colour, amount) {
  amount <- max(-1, min(1, amount))
  rgb <- grDevices::col2rgb(colour) / 255
  target <- if (amount >= 0) 1 else 0
  mixed <- rgb + (target - rgb) * abs(amount)
  grDevices::rgb(mixed[1L], mixed[2L], mixed[3L])
}


.nbt_stable_hash <- function(x) {
  code <- utf8ToInt(enc2utf8(x))
  if (!length(code)) {
    return(0L)
  }
  as.integer(sum((code + 17) * (seq_along(code) + 31)) %% 2147483647)
}


.nbt_semantic_shade <- function(label, family) {
  suffix <- sub(paste0("^", family, "(?:_| \\()?"), "", label)
  suffix <- sub("\\)$", "", suffix)

  if (!nzchar(suffix)) {
    return(0)
  }
  if (suffix == "raw") {
    return(0.34)
  }
  if (grepl("multiread|double_stranded", suffix)) {
    return(-0.30)
  }
  if (grepl("singleton|single_stranded", suffix)) {
    return(0.34)
  }
  NA_real_
}


.nbt_coverage_suffix <- function(label) {
  # Treat a delimited value such as _10x, _30.5x, or _30x_sample as coverage.
  # Plain numeric sample IDs (for example, _7614) deliberately do not match.
  match <- regexec("(?:^|_)([0-9]+(?:\\.[0-9]+)?)x(?:_|$)", label, perl = TRUE)
  value <- regmatches(label, match)[[1L]]
  if (length(value) != 2L) {
    return(NA_real_)
  }
  as.numeric(value[[2L]])
}


.nbt_assign_shades <- function(labels, family) {
  # Zero is the canonical family colour. Remaining values span darker and
  # lighter variants while retaining the family's hue.
  fallback <- c(-0.18, 0.20, -0.36, 0.38, -0.50, 0.52, 0.65)
  assigned <- setNames(rep(NA_real_, length(labels)), labels)

  desired <- vapply(labels, .nbt_semantic_shade, numeric(1L), family = family)
  coverage <- vapply(labels, .nbt_coverage_suffix, numeric(1L))
  coverage_index <- which(!is.na(coverage))

  if (length(coverage_index)) {
    coverage_order <- coverage_index[order(coverage[coverage_index], labels[coverage_index])]

    if (family %in% labels) {
      # With the unsuffixed method present, all explicit coverage variants are
      # lighter than the base, and lower coverage is lightest.
      coverage_shades <- seq(0.48, 0.18, length.out = length(coverage_order))
    } else if (length(coverage_order) == 1L) {
      coverage_shades <- 0.34
    } else {
      # Without an unsuffixed method, span both sides of the base colour while
      # preserving the same lower-is-lighter ordering.
      coverage_shades <- seq(0.46, -0.18, length.out = length(coverage_order))
    }

    desired[coverage_order] <- coverage_shades
  }

  known_order <- order(is.na(desired), labels)
  used <- numeric()

  for (index in known_order) {
    label <- labels[[index]]
    candidate <- desired[[index]]

    if (is.na(candidate)) {
      start <- (.nbt_stable_hash(label) %% length(fallback)) + 1L
      candidates <- fallback[c(
        seq.int(start, length(fallback)),
        if (start > 1L) seq_len(start - 1L) else integer()
      )]
    } else {
      candidates <- c(
        candidate,
        fallback[order(abs(fallback - candidate))]
      )
    }

    available <- candidates[!vapply(
      candidates,
      function(value) any(abs(used - value) < 1e-12),
      logical(1L)
    )]
    if (!length(available)) {
      # More than eight variants in one family is unusual; continue producing
      # deterministic colours rather than silently recycling a shade.
      available <- seq(-0.58, 0.72, length.out = length(labels) + 2L)
      available <- available[!vapply(
        available,
        function(value) any(abs(used - value) < 1e-12),
        logical(1L)
      )]
    }

    assigned[[label]] <- available[[1L]]
    used <- c(used, assigned[[label]])
  }

  assigned
}


#' Return stable colours for sequencing methods and their subtypes
#'
#' When only one label from a technology family is present, it receives the
#' family's base colour. When two or more labels from one family are present,
#' each receives a deterministic shade of that same colour. Known contrasts
#' have semantic direction: multiread/double-stranded are darker, while
#' singleton/single-stranded are lighter. Coverage suffixes are ordered so a
#' smaller value such as _10x is lighter than _30x or _60x. A coverage-suffixed
#' label shown by itself still receives the family's unmodified base colour.
#'
#' @param methods Method labels present in the plot.
#' @param alpha Opacity from 0 to 1.
#' @param other_colour Colour for methods outside the registered panel.
#' @return A named character vector keyed by normalized method label.
seq_method_colours <- function(
  methods,
  alpha = 1,
  other_colour = "#7F7F7F"
) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be one number between 0 and 1.", call. = FALSE)
  }

  labels <- unique(normalize_seq_method(methods))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) {
    return(setNames(character(), character()))
  }

  families <- seq_method_family(labels)
  result <- setNames(rep(NA_character_, length(labels)), labels)

  for (family in unique(families)) {
    members <- labels[families == family]
    base <- if (family %in% names(NBT_METHOD_BASE_COLOURS)) {
      NBT_METHOD_BASE_COLOURS[[family]]
    } else {
      other_colour
    }

    if (length(members) == 1L) {
      result[[members]] <- base
    } else {
      shades <- .nbt_assign_shades(sort(members), family)
      result[names(shades)] <- vapply(
        shades,
        function(amount) .nbt_mix_colour(base, amount),
        character(1L)
      )
    }
  }

  if (alpha < 1) {
    result <- grDevices::adjustcolor(result, alpha.f = alpha)
  }
  result
}

# American-spelling alias.
seq_method_colors <- seq_method_colours


.nbt_scale_spec <- function(methods, alpha, other_colour) {
  raw <- as.character(methods)
  raw <- raw[!is.na(raw)]
  raw <- unique(raw)
  normalized <- normalize_seq_method(raw)
  normalized_colours <- seq_method_colours(
    normalized,
    alpha = alpha,
    other_colour = other_colour
  )

  keep <- !duplicated(normalized)
  list(
    values = stats::setNames(unname(normalized_colours[normalized]), raw),
    breaks = raw[keep],
    labels = normalized[keep]
  )
}


#' Manual colour scale for sequencing methods
#'
#' @param methods The vector mapped to colour. Passing the observed vector lets
#'   the scale decide whether a family needs its base colour or subtype shades.
#' @param alpha Opacity from 0 to 1.
#' @param other_colour Colour used for unregistered method families.
#' @param breaks,labels Optional legend overrides.
#' @param ... Passed to ggplot2::scale_colour_manual().
scale_colour_seq_method <- function(
  methods,
  ...,
  alpha = 1,
  other_colour = "#7F7F7F",
  breaks = NULL,
  labels = NULL
) {
  .nbt_assert_ggplot2()
  spec <- .nbt_scale_spec(methods, alpha, other_colour)
  if (is.null(breaks)) breaks <- spec$breaks
  if (is.null(labels)) labels <- normalize_seq_method(breaks)

  ggplot2::scale_colour_manual(
    ...,
    values = spec$values,
    breaks = breaks,
    labels = labels
  )
}

# American-spelling alias.
scale_color_seq_method <- scale_colour_seq_method


#' Manual fill scale for sequencing methods
#'
#' @inheritParams scale_colour_seq_method
#' @param ... Passed to ggplot2::scale_fill_manual().
scale_fill_seq_method <- function(
  methods,
  ...,
  alpha = 1,
  other_colour = "#7F7F7F",
  breaks = NULL,
  labels = NULL
) {
  .nbt_assert_ggplot2()
  spec <- .nbt_scale_spec(methods, alpha, other_colour)
  if (is.null(breaks)) breaks <- spec$breaks
  if (is.null(labels)) labels <- normalize_seq_method(breaks)

  ggplot2::scale_fill_manual(
    ...,
    values = spec$values,
    breaks = breaks,
    labels = labels
  )
}


#' Minimal Nature Biotechnology font theme
#'
#' By default this changes only the font family. It is deliberately a partial
#' theme: axes, ticks, gridlines, panel backgrounds, spacing, text sizes, title
#' styling, and legend placement are inherited unchanged from the plot.
#'
#' @param base_size Optional font size. NULL preserves the plot's current size.
#' @param base_family Font family.
#' @param line_width Optional root line width. NULL preserves existing lines.
#' @param transparent Optional logical. NULL preserves existing backgrounds;
#'   TRUE makes plot/panel/legend backgrounds transparent; FALSE makes them
#'   white.
#' @param legend_position Optional ggplot2 legend position. NULL preserves the
#'   plot's existing position.
#' @return A partial ggplot2 theme.
theme_nature_biotech <- function(
  base_size = NULL,
  base_family = "Arial",
  line_width = NULL,
  transparent = NULL,
  legend_position = NULL
) {
  .nbt_assert_ggplot2()
  if (!is.null(base_size) &&
      (!is.numeric(base_size) || length(base_size) != 1L || is.na(base_size) || base_size <= 0)) {
    stop("'base_size' must be NULL or one positive number.", call. = FALSE)
  }
  if (!is.null(line_width) &&
      (!is.numeric(line_width) || length(line_width) != 1L || is.na(line_width) || line_width <= 0)) {
    stop("'line_width' must be NULL or one positive number.", call. = FALSE)
  }

  elements <- list(
    text = ggplot2::element_text(family = base_family, size = base_size)
  )
  if (!is.null(line_width)) {
    elements$line <- ggplot2::element_line(linewidth = line_width)
  }
  if (!is.null(transparent)) {
    fill <- if (isTRUE(transparent)) NA else "white"
    # Set only the fill so existing borders and other rectangle properties are
    # inherited unchanged from the plot's current theme.
    elements$plot.background <- ggplot2::element_rect(fill = fill)
    elements$panel.background <- ggplot2::element_rect(fill = fill)
    elements$legend.background <- ggplot2::element_rect(fill = fill)
    elements$legend.key <- ggplot2::element_rect(fill = fill)
  }
  if (!is.null(legend_position)) {
    elements$legend.position <- legend_position
  }

  do.call(ggplot2::theme, elements)
}

# Short aliases for interactive work.
theme_nbt <- theme_nature_biotech
theme_natbiotech <- theme_nature_biotech


#' Add the publication theme and method scales together
#'
#' @param methods Optional vector mapped to a method aesthetic. If NULL, only
#'   the theme is returned.
#' @param aesthetics Any combination of "colour"/"color" and "fill".
#' @inheritParams theme_nature_biotech
#' @param alpha Opacity for method colours.
#' @param other_colour Colour used for unregistered method families.
#' @return A list that can be added directly to a ggplot with `+`.
style_nature_biotech <- function(
  methods = NULL,
  aesthetics = if (is.null(methods)) character() else "colour",
  base_size = NULL,
  base_family = "Arial",
  line_width = NULL,
  transparent = NULL,
  legend_position = NULL,
  alpha = 1,
  other_colour = "#7F7F7F"
) {
  result <- list(theme_nature_biotech(
    base_size = base_size,
    base_family = base_family,
    line_width = line_width,
    transparent = transparent,
    legend_position = legend_position
  ))

  aesthetics <- unique(sub("^color$", "colour", aesthetics))
  unsupported <- setdiff(aesthetics, c("colour", "fill"))
  if (length(unsupported)) {
    stop(
      "Unsupported aesthetic(s): ", paste(unsupported, collapse = ", "),
      ". Use 'colour', 'color', and/or 'fill'.",
      call. = FALSE
    )
  }
  if (length(aesthetics) && is.null(methods)) {
    stop("Supply 'methods' when requesting method colour/fill scales.", call. = FALSE)
  }

  if ("colour" %in% aesthetics) {
    result <- c(result, list(scale_colour_seq_method(
      methods,
      alpha = alpha,
      other_colour = other_colour
    )))
  }
  if ("fill" %in% aesthetics) {
    result <- c(result, list(scale_fill_seq_method(
      methods,
      alpha = alpha,
      other_colour = other_colour
    )))
  }
  result
}

# Short alias.
style_nbt <- style_nature_biotech


#' Save a plot at a Nature-compatible figure width
#'
#' Preserves the plot background by default and uses a Cairo PDF device when
#' available so Arial text remains editable. Nature standard widths are 89 mm
#' (single column) and 183 mm (double column).
#'
#' @param plot A ggplot object.
#' @param filename Output filename; its extension selects the format.
#' @param width "single", "double", or a numeric width in mm.
#' @param height Numeric height in mm. Nature's stated maximum is 170 mm.
#' @param dpi Raster resolution. Nature requests at least 450 dpi for images in
#'   the final artwork guide; vector PDF output does not use this value.
#' @param bg Export background.
#' @param device Optional ggsave device override.
#' @param ... Additional arguments passed to ggplot2::ggsave().
#' @return The normalized output path, invisibly.
save_nature_plot <- function(
  plot,
  filename,
  width = c("single", "double"),
  height = 60,
  dpi = 450,
  bg = NULL,
  device = NULL,
  ...
) {
  .nbt_assert_ggplot2()

  if (is.character(width)) {
    width <- match.arg(width)
    width <- unname(c(single = 89, double = 183)[[width]])
  }
  if (!is.numeric(width) || length(width) != 1L || is.na(width) || width <= 0) {
    stop("'width' must be 'single', 'double', or one positive number in mm.", call. = FALSE)
  }
  if (!is.numeric(height) || length(height) != 1L || is.na(height) || height <= 0) {
    stop("'height' must be one positive number in mm.", call. = FALSE)
  }
  if (height > 170) {
    warning("Nature's stated maximum figure height is 170 mm.", call. = FALSE)
  }

  extension <- tolower(tools::file_ext(filename))
  if (!nzchar(extension)) {
    stop("'filename' must include an output extension such as .pdf or .png.", call. = FALSE)
  }
  if (is.null(device) && extension == "pdf" && capabilities("cairo")) {
    device <- grDevices::cairo_pdf
  }

  args <- list(
    filename = filename,
    plot = plot,
    device = device,
    width = width,
    height = height,
    units = "mm",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE,
    ...
  )
  if (is.null(device)) {
    args$device <- NULL
  }
  do.call(ggplot2::ggsave, args)
  invisible(normalizePath(filename, mustWork = FALSE))
}
