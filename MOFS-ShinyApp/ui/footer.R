ui.footer <- function() {
  tags$footer(tags$hr(),
              HTML("Copyright &copy; 2022"),
              HTML("<span><a href=\"https://rookieutopia.com/\" target=\"_blank\">MOFS</a>, All Rights Reserved -</span>"),
              span("An open platform for enabling clinicians obtain MOFS"),
              tags$br(),
              span("Covered by "),
              HTML("<span><a href=\"https://creativecommons.org/licenses/by-nc/4.0/\" target=\"_blank\">CC BY-NC License</a></span>"),
              span(" | "),
              HTML("<span><a href=\"https://www.beian.miit.gov.cn/\" target=\"_blank\">沪ICP备18048749号-4</a></span>"),
              align = "center", style = "
                           position:relative;
                           bottom:0;
                           width:100%;
                           height:50px;   /* Height of the footer */
                           padding: 10px;
                           z-index: 1000;"
  )
}
