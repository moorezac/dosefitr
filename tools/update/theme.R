#' ggplot2 Theme for dosefitr Package
#' @param font_size TODO: description.
#' @param font_small TODO: description.
#' @param font_large TODO: description.
#' @param linewidth TODO: description.
#' @description
#' A clean ggplot2 theme designed for the dosefitr package. Provides consistent
#'    styling across plotting functions.
#' @examples
#' \dontrun{
#' # Apply dosefitr theme to any ggplot
#' library(ggplot2)
#' ggplot(data, aes(x, y)) +
#'   geom_point() +
#'   theme_dosefitr()
#' }
#' @export
theme_dosefitr <- function(
  font_size = 10,
  font_small = 8,
  font_large = 12,
  linewidth = 1
) {
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = linewidth / ggplot2::.pt),
    axis.text = ggplot2::element_text(size = font_small),
    axis.title.x = ggplot2::element_text(size = font_size),
    axis.title.y = ggplot2::element_text(size = font_size),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key = ggplot2::element_rect(
      colour = "black",
      linewidth = linewidth / ggplot2::.pt
    ),
    legend.key.size = ggplot2::unit(1.25 * font_size, "pt"),
    legend.title = ggplot2::element_text(size = font_size),
    panel.grid = ggplot2::element_line(
      colour = "grey90",
      linewidth = linewidth / ggplot2::.pt
    ),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    title = ggplot2::element_text(size = font_size),
    plot.title = ggplot2::element_text(size = font_size),
    text = ggplot2::element_text(size = font_size)
  )
}
