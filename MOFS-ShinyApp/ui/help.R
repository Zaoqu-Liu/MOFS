# this script render the UI of help page
set_md_path <- function(filename) {
  paste("./doc/help/", filename, sep = "")
}

ui.page_help_usage <- function() {
  tabPanel(
    title = "Usage",
    value = "Usage",
    
    fluidPage(
      id = "HelpPageUsageSubpage",
      shiny::includeMarkdown(set_md_path("usage.md")),
      style = "width:80%;")
    )
}


ui.page_help_citation <- function() {
  tabPanel(
    title = "Citation",
    value = "Citation",
    
    fluidPage(
      shiny::includeMarkdown(set_md_path("citation.md")),
      style = "width:80%;")
    )
}


ui.page_help_privacy_policy <- function() {
  tabPanel(
    title = "Privacy Policy",
    value = "PrivacyPolicy",
      
    fluidPage(
      shiny::includeMarkdown(set_md_path("privacy_policy.md")),
      style = "width:80%;")
    )
}

ui.page_help_terms_and_conditions <- function(){
  tabPanel(
    title = "Terms and Conditions",
    value = "TermsAndConditions",
    
    fluidPage(
      shiny::includeMarkdown(set_md_path("terms_and_conditions.md")),
      style = "width:80%;")
  )
}
