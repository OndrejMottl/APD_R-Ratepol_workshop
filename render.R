#----------------------------------------------------------#
#
#
#                 APD R-Ratepol workshop
#
#                       Render
#
#
#                      O. Mottl
#                        2023
#
#----------------------------------------------------------#
# Render all files to both `.dm` and `.html` files`

library(quarto)
library(rmarkdown)

# helper function o render fie in both formats
render_md_and_html <- function(file_name) {
  quarto::quarto_render(
    input = paste0(file_name, ".qmd"),
    output_format = "gfm",
    output_file = paste0(file_name, ".md")
  )

  rmarkdown::render(
    input = paste0(file_name, ".qmd"),
    output_format = rmarkdown::html_document(
      theme = "readable"
    ),
    output_file = paste0("docs/", file_name, ".html")
  )
}

# README -----
render_md_and_html(file_name = "README")

# workshop_info -----
render_md_and_html(file_name = "workshop_info")

# pre_workshop -----
render_md_and_html(file_name = "pre_workshop")

# step_by_step_guide -----
render_md_and_html(file_name = "step_by_step_guide")