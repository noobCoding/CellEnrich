#' @import shinymaterial
#' @export

CellEnrichUI <- function() {
  require(shinymaterial)
  material_page(
    shinyjs::useShinyjs(),

    # dynamic datatable full width

    tags$head(tags$style(type = "text/css", ".display.dataTable.no-footer{width : 100% !important;}")),
    tags$head(tags$style(type = "text/css", ".markerP {display : none;}")),

    # waitress declare
    use_waitress(color = "#697682", percent_color = "#333333"),

    title = paste0(
      "CellEnrich ",
      "<a href = 'https://github.com/jhk0530/cellenrich' target = '_blank'> ", # github link
      "<i class='material-icons' style = 'font-size:1.3em;'>info</i> </a>" # icon tag
    ),
    nav_bar_fixed = FALSE,
    nav_bar_color = "light-blue darken-1",
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
          material_row(
            material_column(
              material_radio_button(
                input_id = "FCoption",
                label = "Select FoldChange Option",
                choices = c("median", "mean", "zero", "GSVA"),
                selected = "median"
              ),
              width = 4
            ),
            material_column(
              material_radio_button(
                input_id = "plotOption",
                label = "Select Plot Option",
                choices = c("t-SNE", "U-MAP"),
                selected = "t-SNE"
              ),
              width = 4
            ),
            material_column(
              material_radio_button(
                input_id = "genesetOption",
                label = "Select Gene-set",
                choices = c(
                  "Curated", # c2
                  "GeneOntology", # c5 = GO
                  "KEGG", # KEGG
                  "Mouse-KEGG", # Mouse
                  "Mouse-GO"
                ),
                selected = "Curated"
              ),
              width = 4
            )
          ),

          material_row(
            material_column(
              material_number_box(
                input_id = "minGenesetSize",
                label = "Minimum Gene-set Size",
                min_value = 10,
                max_value = 30,
                initial_value = 15,
                step_size = 5
              ),
              width = 4
            ),
            material_column(
              material_number_box(
                input_id = "maxGenesetSize",
                label = "Maximum Gene-set Size",
                min_value = 250,
                max_value = 750,
                initial_value = 500,
                step_size = 5
              ),
              width = 4
            ),
            material_column(
              material_number_box(
                input_id = "qvalueCutoff",
                label = "Q-value threshold",
                min_value = 0,
                max_value = 0.25,
                initial_value = 0.05,
                step_size = 0.01
              ),
              width = 4
            )
          ),
          solvedButton(
            inputId = "StartCellEnrich",
            label = "Start CellEnrich",
            style = "margin-left:45%;",
            onClick = 'console.log("CellEnrich");'
          ),
          depth = 3
        ),
        width = 6, # Centered Layout
        offset = 3
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
            plotOutput("CellBar", height = "700px"), # cell distribution
            width = 6
          ),
          material_button("colorbtn", "toColor", icon = "color_lens", color = "blue lighten-1"),
          material_button("graybtn", "toGray", icon = "clear", color = "blue lighten-1"),
          material_button("freqbtn", "Frequent", icon = "grain", color = "blue lighten-1"),
          material_button("sigbtn", "Significant", icon = "grade", color = "blue lighten-1")
        ),
        width = 12
      ),
      style = "margin : 1em"
    ),

    # marker table
    material_row(
      material_card(
        title = "MarkerGenes",
        p("DE from each Cell specific", class = 'markerP'),
        DT::dataTableOutput("markerL1"),
        p("DE - Pathway from each Cell specific", class = 'markerP'),
        DT::dataTableOutput("markerL2")
      ),
      style = "margin : 1em"
    ),

    # emphasize tables
    material_row(
      material_card(
        title = "",
        material_card(
          title = "Pathway Emphasize", divider = TRUE,
          tags$h3("To be recognized by application, Please move element's position"),
          rank_list(text = "Pathways", labels = "Please Clear First", input_id = "sortList", css_id = "mysortableCell"),
          material_button("OrderEmphasize", "Emphasize with Order", icon = "timeline"),
          material_button("Emphasize", "Emphasize without Order", icon = "bubble_chart"),
          material_button("ClearList", "Clear List", icon = "clear_all"),
        ),
        uiOutput("dynamicTable"),
        depth = 3
      ),
      style = "margin : 1em"
    )
  )
}
