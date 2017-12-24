# load all packages
library(shiny)
library(stringr)
library(stringi)
library(purrr)
library(rebus)
library(dplyr)
library(ggplot2)
library(tm)
library(qdap)
library(ggrepel)
source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")
library(graph)
library(igraph)
library(plotrix)
library(googleVis)

# load datasets
load(file = "./data/all_cardiovascular_case_reports.Rdata") # all cardiovascular case reports
load(file = "./data/cv_MH.Rdata") # all cardiovascular MeSH terms

# use a theme for webpage
ui <- fluidPage(theme = "bootstrap.css",

    titlePanel("Associations in MeSH and RN (drug) terms between two population groups"), br(), 
      tabsetPanel(
        tabPanel("Home", sidebarLayout(sidebarPanel(
        
          textInput(inputId = "phrases",
                    label = "Please enter all the MeSH terms for your first search, each separated by a comma:",
                    value = ""),
          
          helpText("Example: female, infant, stroke volume"),
          
          br(),
          
          textInput(inputId = "compare_with",
                    label = "Please enter all the MeSH terms for your second search, each separated by a comma:",
                    value = ""),
          
          helpText("Example: male, middle-aged, aneurysm"),
          
          br(),
          
          sliderInput(inputId = "CR_num",
                      label = "How many case reports should I analyze?",
                      value = "20000",
                      step = 5000,
                      min = 1,
                      max = 428036),
          
          helpText("The maximum number is 428,036. The more you choose, the more accurate the results will be; however, 
                          it will take up to 6 minutes to load. 20,000 is a feasible starting point."), br(), br(), br()),
                 
          mainPanel( 
                 h2("Welcome to this innovative tool!"), br(),
                 p("Instructions:"), 
                 p(("Please search for words/phrases that are in the MeSH browser, which are found at this site:"), a("MeSH Browser", href = "https://meshb.nlm.nih.gov/search?searchInField=allTerms&sort=&size=20&searchType=exactMatch&searchMethod=FullWord", target = "_blank")),
                 p("Your searches are not case-sensitive. At this time, you must search 'middle aged' as 'middle-aged.' 
                   Also, for simplicity, we will refer to the search you type in first as 'first search' and the search to which you compare as 'second search.' 
                   You can navigate between tabs to see different outputs. You can also customize your outputs within each tab."),
                 p("If you run into any errors, then either one of your MeSH phrases is not recognized 
                   or that there is not enough significance between the two groups. Also, please be patient and allow the graphics to load; big data science is at work!"),
                 p("This technology randomly samples the number of case reports you choose each time you run a new search. This means that, unless you choose to load 428,036 case reports,
                   the outputs will vary slightly from search to search, even when using the same number of case reports. This is to prevent any bias to occur from a certain selection of case reports."),
                 br(),
                 p("Background:"),
                 p("There are hundreds of thousands of cardiovascular case reports from PubMed. Each case report contains MeSH terms (MEdical Subject Headings) that 
                   describe demographics, symptoms, diseases, and drugs of a patient, and RN terms (Registry Number) which describe
                   drugs, chemicals, or other substances mentioned in the case report. This tool analyzes which MeSH and RN terms are most unique to certain demographics via the two groups you choose to compare.
                   This tool aids in improveing clinical diagnostics by customizing treatments."),
                 br(), 
                 p("Author:"),
                 p(("Kitu Komya, an undergraduate student at UCLA, working in Dr. Peipei Ping's lab at the NIH - HeartBD2K Lab at UCLA, 
                   under Dr. Harry Caufield's supervision. Peek at my"), a("LinkedIn", href = "https://tinyurl.com/linkedin-kitu", target = "_blank"), ("or e-mail me at kitu.komya@gmail.com! I don't bite! :^)"))))), 
        
        tabPanel("Heat Maps", sidebarLayout(sidebarPanel(
              sliderInput(inputId = "MH_num",
                          label = "How many top, distinct MeSH terms should I plot?",
                          value = "12",
                          min = 1,
                          max = 30),

              helpText("10 or 12 is standard. Any more, and the graphic will become too cluttered."),

              br(),

              sliderInput(inputId = "RN_num",
                          label = "How many top, distinct RN terms should I plot?",
                          value = "12",
                          min = 1,
                          max = 30),

              helpText("10 or 12 is standard. Any more, and the graphic will become too cluttered.")),
          
        mainPanel(
          textOutput("text1"), textOutput("text2"), br(), 
          p("These heat maps are useful in determining the most distinct MeSH or RN terms between each group you searched. In order to interpret the results,
                   you may simply dictate: 'The MeSH term (a) is (b) times more prevalent in case reports concerning (d) than in case reports concerning (d)', 
            where (a) is the MeSH term the line corresponds to, (b) is the numeric factor value, (c) is the search in which it's more prevalent, and (d) is the opposite search.
            Note that numbers in black represent distinctness in the first search, whereas the white numbers correspond to the second search."), br(),
          plotOutput("MH"), plotOutput("RN")))),
        
        tabPanel("Tree Maps", sidebarLayout(sidebarPanel(
          sliderInput(inputId = "treemap1",
                      label = "How many top MeSH terms should I plot from your first search?",
                      value = "20",
                      min = 2,
                      max = 100),
          
          helpText("20 is standard. Any more, and the graphic will become too cluttered."),
          
          br(),
          
          sliderInput(inputId = "treemap2",
                      label = "How many top MeSH terms should I plot from your second search?",
                      value = "20",
                      min = 2,
                      max = 100),
          
          helpText("20 is standard. Any more, and the graphic will become too cluttered.")),
          
          mainPanel(
                 textOutput("text3"), textOutput("text4"), br(), 
                 p("Currently, the treemaps will open in two different browser tabs. These treemaps represent the top MeSH terms within each search. The more
                   pink the value, the more distinct the MeSH term is to that search. The more blue it is, the less distinct it is, and in fact, more similar to the other search.
                   You can hover over a box to view the entire MeSH term name, if it's trimmed."), br(),
                 textOutput("jaccard"), plotOutput("colorbox"), plotOutput("treemap_1"), plotOutput("treemap_2")))),
        
        tabPanel("PCA Plots", sidebarLayout(sidebarPanel(
              sliderInput(inputId = "pca_num",
                          label = "How many top, distinct MeSH terms should I plot from your first search?",
                          value = "20",
                          min = 1,
                          max = 100),

              helpText("20 is standard. Any more, and the graphic will become too cluttered."),
              
              br(),
              
              sliderInput(inputId = "pca_num_against",
                          label = "Please choose how many top, distinct MeSH terms should I plot from your second search?",
                          value = "20",
                          min = 1,
                          max = 100),
              
              helpText("20 is standard. Any more, and the graphic will become too cluttered.")),
                 
              mainPanel(
                 textOutput("text5"), textOutput("text6"), br(), 
                 p("These PCA plots show the spatial relations among the top MeSH terms within each search. The MeSH terms are clustered by k means, meaning that MeSH terms most 
                   similar to each other (and hence closer to each other in the plot) are of the same color."), br(),
                 plotOutput("pca_MH_search"), br(), br(), plotOutput("pca_MH_against")))),
        
        tabPanel("Correlation Plots", sidebarLayout(sidebarPanel(
              sliderInput(inputId = "corr_search",
                          label = "How many top, distinct MeSH terms should I plot from your first search?",
                          value = 20,
                          min = 2,
                          max = 100),

              helpText("20 is standard. Any more, and the graphic will become too cluttered."),

              br(),

              sliderInput(inputId = "corr_against",
                          label = "How many top, distinct MeSH terms should I plot from your second search?",
                          value = 20,
                          min = 2,
                          max = 100),

              helpText("20 is standard. Any more, and the graphic will become too cluttered.")),

              mainPanel(
                 textOutput("text7"), textOutput("text8"), br(), 
                 p("These correlation plots show the relations of various clusters to each other. However, this method does not read in phrases, and only singular words,
                   yet this plot is still useful in understanding how some phrases might be related to each other in one search versus in the other."), br(),
                 plotOutput("corr_MH_search"), br(), plotOutput("corr_MH_against"))))
  )
)
