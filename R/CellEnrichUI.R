#' @import shinymaterial
#' @import highcharter
#' @export

CellEnrichUI <- function() {
  require(shinymaterial)
  require(highcharter)
  material_page(
    shinyjs::useShinyjs(),

    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}")),

    # waitress declare
    use_waitress(color = "#697682", percent_color = "#333333"),

    title = paste0(
      "CellEnrich ",
      "<a href = 'https://github.com/jhk0530/cellenrich' target = '_blank'> ", # github link
      "<i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>" # icon tag
    ),
    nav_bar_fixed = FALSE,
    nav_bar_color = "blue darken-2",
    font_color = "#ffffff",
    include_fonts = FALSE,
    include_nav_bar = TRUE,
    include_icons = FALSE,

    # CellEnrich options
    material_row(
      material_column(
        material_card(
          title = "Options",
          divider = TRUE,
          style = 'border : solid 0.5em #1976d2',
          material_row(
            material_column(
              material_card(
                material_radio_button(
                  input_id = "FCoption",
                  label = "Select FoldChange Option",
                  choices = c("median", "mean", "zero", "GSVA"),
                  selected = "median",
                  color = '#1976d2'
                ),
                material_radio_button(
                  input_id = "plotOption",
                  label = "Select Plot Option",
                  choices = c("t-SNE", "U-MAP"),
                  selected = "t-SNE",
                  color = '#1976d2'
                )
              ),
              width = 4
            ),
            material_column(
              material_card(
                material_number_box(
                  input_id = "minGenesetSize",
                  label = "Minimum Gene-set Size",
                  min_value = 10,
                  max_value = 30,
                  initial_value = 15,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "maxGenesetSize",
                  label = "Maximum Gene-set Size",
                  min_value = 250,
                  max_value = 750,
                  initial_value = 500,
                  step_size = 5
                ),
                material_number_box(
                  input_id = "qvalueCutoff",
                  label = "Q-value threshold",
                  min_value = 0,
                  max_value = 0.25,
                  initial_value = 0.05,
                  step_size = 0.01
                )
              ),
              width = 4
            ),
            material_column(
              material_card(
                material_radio_button(
                  input_id = "genesetOption",
                  label = "Select Gene-set",
                  color = '#1976d2',
                  choices = c(
                    "User-Geneset",
                    "Human-Curated", # c2
                    "Human-KEGG", # KEGG
                    "Human-GO",
                    "Human-GO-BP",
                    "Human-GO-CC",
                    "Human-GO-MF",
                    "Mouse-KEGG", # Mouse
                    "Mouse-GO",
                    "Mouse-GO-BP",
                    "Mouse-GO-CC",
                    "Mouse-GO-MF"
                  ),
                  selected = NULL
                )
              ),
              width = 4
            )
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "Start CellEnrich",
            style = "margin-left:45%; background-color: #1976d2",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 6,
        offset = 3 # center half layout
      )
    ),

    # tSNE/UMAP plot
    material_row(
      material_column(
        material_card(
          title = "Group plot / Distribution",
          depth = 3,
          material_column(
            plotOutput("CellPlot", height = "700px"),
            width = 6
          ),
          material_column(
            highchartOutput("CellBar", height = "700px"), # cell distribution
            width = 6
          ),
          material_button("colorbtn", "toColor", icon = "color_lens", color = "blue darken-2"),
          material_button("graybtn", "toGray", icon = "clear", color = "blue darken-2"),
          material_button("freqbtn", "Frequent", icon = "grain", color = "blue darken-2"),
          material_button("sigbtn", "Significant", icon = "grade", color = "blue darken-2")
        ),
        width = 12
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # marker table
    material_row(
      material_card(
        title = "MarkerGenes",
        material_row(
          material_column(
            material_card(
              title = 'DE from each Cell specific',
              DT::dataTableOutput("markerL1")
            ),
            width = 6
          )
          ,
          material_column(
            material_card(
              title = 'DE - Pathway from each Cell specific',
              DT::dataTableOutput("markerL2")
            ),
            width = 6
          )
        )
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    ),

    # emphasize tables
    material_row(
      material_card(
        title = "",
        material_card(
          title = "Pathway Emphasize", divider = TRUE,
          tags$h3("To be recognized by application, Please move element's position"),
          rank_list(text = "Pathways", labels = "Please Clear First", input_id = "sortList", css_id = "mysortableCell"),
          material_button("OrderEmphasize", "Emphasize with Order", icon = "timeline", color = 'blue darken-2'),
          material_button("Emphasize", "Emphasize without Order", icon = "bubble_chart", color = 'blue darken-2'),
          material_button("ClearList", "Clear List", icon = "clear_all", color = 'blue darken-2'),
        ),
        uiOutput("dynamicTable"),
        depth = 3
      ),
      style = "margin : 1em; border : solid 0.5em #1976d2"
    )
  )
}
