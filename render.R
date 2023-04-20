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


#----------------------------------------------------------#
# 1. README -----
#----------------------------------------------------------#
quarto::quarto_render(
  input = "README.qmd",
  output_format = "gfm",
  output_file = "README.md"
)

rmarkdown::render(
  input = "README.qmd",
  output_format = rmarkdown::html_document(
    theme = "readable"
  ),
  output_file = "docs/index.html"
)


#----------------------------------------------------------#
# 2. workshop_info -----
#----------------------------------------------------------#
quarto::quarto_render(
  input = "workshop_info.qmd",
  output_format = "gfm",
  output_file = "workshop_info.md"
)

rmarkdown::render(
  input = "workshop_info.qmd",
  output_format = rmarkdown::html_document(
    theme = "readable",
    self_contained = TRUE
  ),
  output_file = "docs/workshop_info.html"
)


#----------------------------------------------------------#
# 3. step_by_step_guide -----
#----------------------------------------------------------#
quarto::quarto_render(
  input = "step_by_step_guide.qmd",
  output_format = "gfm",
  output_file = "step_by_step_guide.md"
)

rmarkdown::render(
  input = "step_by_step_guide.qmd",
  output_format = rmarkdown::html_document(
    theme = "readable"
  ),
  output_file = "docs/step_by_step_guide.html"
)
