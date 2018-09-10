require(ggthemes)

source("../../colors/color_map.R")

all_categories = names(category_colors)
all_categories = factor(all_categories, levels=all_categories, ordered=TRUE)

print(all_categories)

category_abs_colors <- primary_colors[unlist(category_colors)]
names(category_abs_colors) <- names(category_colors) 

scale_colour_skp <- function(...) {
  scale_colour_manual(values=category_abs_colors, labels=category_names, ...)
}
scale_color_skp <- scale_colour_skp

scale_fill_skp <- function(...) {
  scale_fill_manual(values=category_abs_colors, labels=category_names, ...)
}

