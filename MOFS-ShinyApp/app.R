# -----------------------------Declaration---------------------------------
# 
# MOFS
#
#
# -----------------------------Author Information--------------------------
# 
# Author: Mingjie Wang
# Email: huzai920621@126.com
# Affiliation: Shanghai Tengyun Biotech. co., Ltd. 
#
# ----------------------------Software updates-----------------------------
# 
# Version: 0.1.1
# Date: 2023-02-08
# Log: initial setup
# 
#

# -----Packages-----
message("[+] Checking depedencies...")
source("./lib/load_packages.R")

# -----Setting up-----
options(shiny.maxRequestSize=1024*1024^2)
message("[+] Starting...")



# -----UI part-----
## ------UI files-----
# read main page
pages_file <- dir("./ui", pattern = "\\.R$", full.names = TRUE)
sapply(pages_file, function(x, y) source(x, local = y), y = environment())

## -----Set up UI-----
ui <- shinyUI(
  fluidPage(
    tags$head(
      tags$title("MOFS"),
      tags$link(rel = "stylesheet", type = "text/css", href = "./static/css/mystyle.css"),
      tags$link(rel = "shortcut icon", href = "logo.ico")
	),
    
    shinyjs::useShinyjs(),
    
    #navbar
    bslib::page_navbar( # use bslib to place tabs to the right, see https://github.com/rstudio/bslib/pull/319
      id = "navbar",
      bg = "#00b4d8",
      # title = div(
      #   tags$a(img(src = "logo_short_small_white.png", height = 40, width = 35, style="padding-top:0px;padding-bottom:0px;"),href="")
      # ),
      
      ##### web pages ####
      ## 自定义css控制nav字体大小
      !!!list(
        
        # nav spacer moves all tabs to the right
        nav_spacer(),
        
        # home page
        nav("Home",ui.page_home(),icon = icon("home")),
        
        # help page, four sub pages
        # nav_menu("Help",
        #          icon = icon("question-circle"),
        #          nav("Usage",ui.page_help_usage()),
        #          nav("Citation",ui.page_help_citation()),
        #          # nav("FAQ",ui.page_help_faq())
        #          ),
        # 
        # about page
        # nav("About", ui.page_about(), icon = icon("address-book")),
        # nav(NULL, ui.page_about()),
        
        # contact page
        nav("Contact",ui.page_contact(),icon = icon("paper-plane")),
        
        nav(title = NULL,value="PrivacyPolicy",ui.page_help_privacy_policy()),
        nav(title = NULL,value = "TermsAndConditions",ui.page_help_terms_and_conditions())
      ),
      footer = ui.footer(),
      collapsible = TRUE,
      theme = bs_theme(bootswatch = "yeti")
    )
  )
)
  

# -----Server part-----
server <- function(input, output, session) {
  # pages in nav bar
  source("./server/home.R", local = TRUE)
  source("./server/contact.R", local = TRUE)
  # source("./server/help.R", local = TRUE)
  
  # success info
  message("[+] Shiny app run successfully! Enjoy it!\n") 
}

# ----- Main wrapper -----
shiny::shinyApp(
  ui = ui,
  server = server
)
