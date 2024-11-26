# this script render the UI of home page
# connected to ./server/home.R

# read in home page introduction
home.intro <- read_lines("./db/home_introduction.txt")
home.intro.html <- paste(home.intro, "</p>  <p>", collapse = " ")
home.intro.html <- paste("<p>", home.intro.html)
home.intro.html <- str_replace(home.intro.html, "</p>  <p>$", "</p>")
home.intro.html <- str_remove_all(home.intro.html, "<p>\\s+</p>")
#print(home.intro.html)

# read in news and updates
home.news.and.updates <- read_lines("./db/home_news_and_updates.txt")
home.news.and.updates.html <- paste(home.news.and.updates, "</p>  <p>", collapse = " ")
home.news.and.updates.html <- paste("<p>", home.news.and.updates.html)
home.news.and.updates.html <- str_replace(home.news.and.updates.html, "</p>  <p>$", "</p>")
home.news.and.updates.html <- str_remove_all(home.news.and.updates.html, "<p>\\s+</p>")


# read in citation
# home.citation <- read_lines("./db/home_citation.txt")
# home.citation.html <- paste(home.citation, "</p>  <p>", collapse = " ")
# home.citation.html <- paste("<p>", home.citation.html)
# home.citation.html <- str_replace(home.citation.html, "</p>  <p>$", "</p>")
# home.citation.html <- str_remove_all(home.citation.html, "<p>\\s+</p>")


# home page carousel images
home.carousel.images = list.files("./www/home_carousel/") # all imgs contained in this path will be explored
home.carousel.images.url = paste("./www/home_carousel/", home.carousel.images, sep = "")


# set up with UI
ui.page_home <- function() {
  tabPanel(
    title = "Home",
    value = "Home",
    
    # create icon http://shiny.rstudio.com/reference/shiny/latest/icon.html
    fluidPage(
      style = "width:80%;",
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "./static/css/mystyle.css")
      ),
      tags$br(),
      div(
        img(src = "logo_single_small.png", height = 80, width = 85, style="padding-top:0px;padding-bottom:0px;"),
        p("M O F S", style = "font-weight:bold;font-size:220%;color:#276e7f;padding-bottom:0px;margin-bottom:0px;font-family:\"Comic Sans MS\", \"Comic Sans\", cursive;"),
        tags$p("Multi-omics fusion subtype", 
               style = "font-size:110%;color:#999;font-style:italic;padding-top:0px;margin-top:0px;"),
        
        style = "text-align:center;"
      ),
      
      #----section1-----
      fluidRow(
        # column(
        #   6,
        #   tags$div(
        #     div(
        #       HTML(home.intro.html), 
        #       class = "scrollbar-hidden",
        #       style = "height:400px; overflow-y:auto;color:white; font-size:115%;text-align:justify;padding:30px 50px 30px 50px;border-radius:10px;background-color:#296D7F;"
        #     ),
        #     
        #     id = "HomePageIntroduction",
        #     style = "display:flex;align-items:center;height:450px;box-shadow:5px 5px 10px 5px #ccc;border-radius:10px;background-color:#296D7F;"
        #   )
        # ),
        column(
          12,
          div(
            div(
              slickR(obj = home.carousel.images.url,
                     height = 390, 
                     width = "95%") + 
                settings(dots = FALSE, autoplay = FALSE),
            ),
            style = "height:450px;box-shadow:5px 5px 10px 5px #ccc;padding:20px 0px 20px 0px;border-radius:10px;"
          )
        )
      ),
      
      
      #tags$br(),
      #tags$hr(),
      tags$br(),
      
      #### Section 2 info #####
      tags$hr(),
      tags$p("Please input your multimodal MRI feature file.",
             style = "font-size:120%;font-weight:bold;"),
      tags$p("Note: Only one sample at a time. You can use the InputExample data
             to make adjustments",
             style = "font-size:90%; color:#777777"),
      tags$br(),
      
      div(
        fluidRow(
          column(
            4,
            fluidRow(
              column(
                7,
                shiny::fileInput(
                  inputId = "inputFileBtn",
                  multiple = FALSE,
                  label = NULL,
                  accept = ".csv",
                  buttonLabel = "Select file"
                )
              ),
              column(
                5,
                shiny::downloadButton(
                  outputId = "demoDownloadBtn",
                  label = "InputExample",
                  icon = NULL,
                  style = "border-width:0px;background-color:#fff;font-style:italic;text-decoration:underline;"
                )
              )
            )
            
          ),
          column(
            2,
            # shiny::actionButton(
            #   inputId = "runBtn",
            #   label = "Run"
            # )
            shinyWidgets::actionBttn(
              inputId = "runBtn",
              label = "Run",
              icon = icon("rocket"),
              style = "pill",
              color = "success",
              block = TRUE
            )
          ),
          column(
            1
          ),
          column(
            4,
            shiny::htmlOutput(
              outputId = "resultOutput"
            )
          )
        )
      ),
      
      tags$br(),
      tags$hr(),
      tags$br(),
      
      #### section 4 news and update ####
      div(
        bs4Dash::box(
          
          div(
            HTML(home.news.and.updates.html),
            style="height:170px;overflow-y:scroll;border:1px solid #cecece;padding:10px 20px 10px 20px;text-align:justify;"
          ),
          
          
          # box configure
          id = "HomePageLastInfoCard",
          title = "News and updates",
          solidHeader = FALSE,
          height = "200px",
          closable = FALSE,
          maximizable = FALSE,
          width = 12,
          collapsible = FALSE,
          icon = icon("dove")
        ),
        
        style = "box-shadow:2px 2px 5px 2px #ccc;"
      ),
      
      
      tags$br(),
      #### section 5 citation ####
      # div(
      #   bs4Dash::box(
      #     
      #     div(
      #       
      #       HTML(home.citation.html),
      #       
      #       style="height:170px;overflow-y:scroll;border:1px solid #cecece;padding:10px 20px 10px 20px;text-align:justify;"
      #     ),
      #     
      #     # box configure
      #     id = "HomePageLastInfoCard",
      #     title = "Citation",
      #     solidHeader = FALSE,
      #     height = "200px",
      #     closable = FALSE,
      #     maximizable = FALSE,
      #     width = 12,
      #     collapsible = FALSE,
      #     icon = icon("book")
      #   ),
      #   
      #   style = "box-shadow:2px 2px 5px 2px #ccc; display:none;"
      # ),
      
      tags$br(),
      tags$br()
      
      
    )
  )
}
