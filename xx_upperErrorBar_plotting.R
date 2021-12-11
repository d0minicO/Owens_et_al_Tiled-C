geom_uperrorbar <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}



GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
                          default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
                                            alpha = NA),
                          
                          draw_key = draw_key_path,
                          required_aes = c("x", "y", "ymax"),
                          
                          setup_data = function(data, params) {
                            data$width <- data$width %||%
                              params$width %||% (resolution(data$x, FALSE) * 0.9)
                            
                            transform(data,
                                      xmin = x - width / 2, xmax = x + width / 2, width = NULL
                            )
                          },
                          draw_panel = function(data, panel_scales, coord, width = NULL) {
                            GeomPath$draw_panel(data.frame(
                              x = as.vector(rbind(data$xmin, data$xmax, NA, data$x,   data$x)),
                              y = as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y)),
                              colour = rep(data$colour, each = 5),
                              alpha = rep(data$alpha, each = 5),
                              size = rep(data$size, each = 5),
                              linetype = rep(data$linetype, each = 5),
                              group = rep(1:(nrow(data)), each = 5),
                              stringsAsFactors = FALSE,
                              row.names = 1:(nrow(data) * 5)
                            ), panel_scales, coord)
                          }
)
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}
                          