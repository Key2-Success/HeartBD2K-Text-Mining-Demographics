# load packages
library(shiny)
library(purrr)
library(stringi)
library(stringr)
library(dplyr)
library(rebus)
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
load(file = "./data/named_RN.Rdata") # all RN terms
load(file = "./data/final_RN.Rdata") # all cardiovascular case reports for RN

# back-end code goes here
server <- function(input, output)
{
  output$text1 <- renderText({
  paste("You have chosen: ", input$phrases)})
  
  output$text2 <- renderText({
    paste("You have chosen to compare with: ", input$compare_with)})
  
  output$text3 <- renderText({
    paste("You have chosen: ", input$phrases)})
  
  output$text4 <- renderText({
    paste("You have chosen to compare with: ", input$compare_with)})
  
  output$text5 <- renderText({
    paste("You have chosen: ", input$phrases)})
  
  output$text6 <- renderText({
    paste("You have chosen to compare with: ", input$compare_with)})
  
  output$text7 <- renderText({
    paste("You have chosen: ", input$phrases)})
  
  output$text8 <- renderText({
    paste("You have chosen to compare with: ", input$compare_with)})
  
  
    # dataframe of # of CR as decided by user
    df_cr <- reactive({
      final2 <- final[sample(nrow(final), input$CR_num, replace = F), ]
      return(final2)
    })
  
    
    # dataframe of user's search
    df <- reactive({
      
      final2 <- df_cr()
      # cleans user's search
      a <- paste0("\\b", sapply(strsplit(input$phrases, ", ")[[1]],  function(x) tolower(noquote(x))), "\\b")
      
      # returns user's dataframe
      indices <- colSums(do.call(rbind,lapply(a, function(x) grepl(pattern = x, x = final2$search, ignore.case = TRUE))))==length(a)
      return(final2[indices, ])
    })
    
    # dataframe of user's compare with search
    df2 <- reactive({
      final3 <- df_cr()
      
      # repeat for compare_with
      b <- paste0("\\b", sapply(strsplit(input$compare_with, ", ")[[1]],  function(x) tolower(noquote(x))), "\\b")
      
      # returns user's dataframe
      indices2 = colSums(do.call(rbind,lapply(b, function(x) grepl(pattern = x, x = final3$search, ignore.case = TRUE))))==length(b)
      return(final3[indices2, ])
      
    })
    
    
    # calculate jaccard index similarity in R
    jaccard <- reactive({
      # obtain dataframes of user and bind them
      df <- df()
      df2 <- df2()
      test <- rbind(df, df2)
      
      # find number of unique intersections between both searches
      intersect <- nrow(unique(test[test$PMID %in% intersect(df$PMID, df2$PMID),]))
      total <- nrow(unique(test))
      
      # find jaccard index, and round to 3 sigfigs
      jaccard <- intersect/total
      jaccard <- lapply(jaccard, round, 3)
      
      # return value
      text <- paste("Your Jaccard index is ", jaccard, ". This value ranges from 0 to 1 and shows how similar the two searches are (1 being both searches are identical to each other).")
      return(text)
      
    })
    
    # return just the jaccard index number
    jaccard_num <- reactive({
      # obtain dataframes of user and bind them
      df <- df()
      df2 <- df2()
      test <- rbind(df, df2)
      
      # find number of unique intersections between both searches
      intersect <- nrow(unique(test[test$PMID %in% intersect(df$PMID, df2$PMID),]))
      total <- nrow(unique(test))
      
      # find jaccard index, and round to 3 sigfigs
      jaccard <- intersect/total
      jaccard <- lapply(jaccard, round, 3)
      
      # return value
      return(jaccard)
    })
    
    # returns occurrences of all cv_MH in both searches
    cv_occur <- reactive({
      
      # uses user dataframe
      result <- df()
      result2 <- df2()
      
      # size of both dataframes: find p-value for 2 prop z test
      size <- dim(result)[1]
      size2 <- dim(result2)[1]
      
      # cleans into a vector
      result <- as.data.frame(result[, 1:3])
      result_MH <- result[ , 3]
      names(result_MH)[1] <- "MH"
      
      # cleans into a vector
      result2 <- as.data.frame(result2[ , 1:3])
      result_MH2 <- result2[ , 3]
      names(result_MH2)[1] <- "MH"
      
      # find occurrences of initial dataframe
      cv_MH$search <- map_int(cv_MH$V1, function(x){sum(str_detect(result_MH, pattern = x))})
      
      # find occurrences of initial dataframe
      cv_MH$against <- map_int(cv_MH$V1, function(x){sum(str_detect(result_MH2, pattern = x))})
      
      # delete all duplicate entries
      cv_MH <- subset(cv_MH, !duplicated(cv_MH$V1))
      
      # proportions
      cv_MH$search_prop <- cv_MH$search/size
      cv_MH$against_prop <- cv_MH$against/size2
      
      return(cv_MH)
    })
    
    # treemap for first search
    treemap_1 <- reactive({
      cv_MH <- cv_occur()
      
      # take and store as character variables the searches
      input_phrases <- input$phrases
      input_phrases <- as.character(input_phrases)
      
      # make user input numeric
      # treemap1 <- as.numeric(input$treemap1)
      
      # input_comparewith <- input$compare_with
      # input_comparewith <- as.character(input_comparewith)
      
      # create new df for search 1
      search1 <- cv_MH[order(-cv_MH$search_prop), ] # top occurrences
      search1 <- search1[1:input$treemap1, ] # treemap df for first search
      search1$usersearch <- input_phrases
      search1 <- search1[ , c(1, 6, 4, 5)]
      # search1 <- search1[ , c(1, 4, 2, 3)]
      # add1 <- data.frame(v1 = c(input_phrases), search = c(80), against = c(90), search_prop = c(20), against_prop = c(20), search = c(NA))
      # names(add1) <- names(search1)
      # search1 <- rbind(search1, add1)
      
      # # create new df for search 2
      # search2 <- cv_MH[order(-cv_MH$against_prop), ] # top occurrences
      # search2 <- search2[1:20, ] # treemap df for second search
      # search2$usersearch <- input$compare_with
      # 
      # # reverse variables for the second dataframe
      # search2$temp <- search2$search_prop
      # search2$search_prop <- search2$against_prop
      # search2$against_prop <- search2$temp
      # search2 <- search2[ , -c(7)]
      # 
      # names(search2) <- names(search1)
      # search2 <- search2[ , c(1, 4, 2, 3)]
      # add2 <- data.frame(v1 = c(input_comparewith), search = c(10), against = c(10), search_prop = c(20), against_prop = c(20), search = c(NA))
      # names(add2) <- names(search2)
      # search2 <- rbind(search2, add2)

      # # combine both dataframes for treemap
      # search_all <- rbind(search1, search2)
      # search_all <- search_all[ , c(1, 6, 4, 5)]
      
      # add "parental" nodes
      add1 <- data.frame(input$phrases, NA, 20, 20)
      # add2 <- data.frame(input$compare_with, NA, 20, 20)
      names(add1) <- names(search1)
      # names(add2) <- names(search_all)
      
      # confirm NA is there
      add1$usersearch <- as.character(add1$usersearch)
      add1$usersearch[1] <- NA
      
      # add2$usersearch <- as.character(add2$usersearch)
      # add2$usersearch[1] <- NA
      
      # bind dataframes with the "parental" nodes dataframes
      search_all <- rbind(search1, add1)
      # search_all <- rbind(search_all, add2)
      
      # change class of first two variables to factor
      search_all$V1 <- as.factor(search_all$V1)
      search_all$usersearch <- as.factor(search_all$usersearch)

      # create treemap
      allMap <- gvisTreeMap(data = search_all, idvar = "V1", parentvar = "usersearch",
                            sizevar = "search_prop", colorvar = "against_prop",
                            options = list(showScale = TRUE, minColor = "#F8766D", maxColor = "#619CFF"))
      plot(allMap)

      # return(search_all)
    })
    

    # treemap for second search
    treemap_2 <- reactive({
      cv_MH <- cv_occur()
      
      # take and store as character variables the searches
      input_phrases <- input$compare_with
      input_phrases <- as.character(input_phrases)
      
      # create new df for search 1
      search1 <- cv_MH[order(-cv_MH$against_prop), ] # top occurrences
      search1 <- search1[1:input$treemap2, ] # treemap df for second search
      search1$usersearch <- input_phrases
      search1 <- search1[ , c(1, 6, 4, 5)]

      # add "parental" nodes
      add1 <- data.frame(input$compare_with, NA, 20, 20)
      names(add1) <- names(search1)
      
      # confirm NA is there
      add1$usersearch <- as.character(add1$usersearch)
      add1$usersearch[1] <- NA
      
      # bind dataframes with the "parental" nodes dataframes
      search_all <- rbind(search1, add1)
      
      # change class of first two variables to factor
      search_all$V1 <- as.factor(search_all$V1)
      search_all$usersearch <- as.factor(search_all$usersearch)
      
      # create treemap
      allMap <- gvisTreeMap(data = search_all, idvar = "V1", parentvar = "usersearch",
                            sizevar = "against_prop", colorvar = "search_prop",
                            options = list(showScale = TRUE, minColor = "#F8766D", maxColor = "#619CFF"))
      plot(allMap)

    })
    
    # MH terms heatmap
    df_MH <- reactive({
      # uses user dataframe
      result <- df()
      result2 <- df2()

      # cleans into a vector
      result <- as.data.frame(result[, 1:3])
      result_MH <- result[ , 3]
      names(result_MH)[1] <- "MH"
      
      # returns cv dataframe
      cv_MH <- cv_occur()
      
      # cleans into a vector
      result2 <- as.data.frame(result2[ , 1:3])
      result_MH2 <- result2[ , 3]
      names(result_MH2)[1] <- "MH"
      
      # size of both dataframes: find p-value for 2 prop z test
      size <- dim(result)[1]
      size2 <- dim(result2)[1]
      
      # remove chi squared warnings by finding and removing those that don't satisfy the conditions
      count <- numeric(dim(cv_MH)[1]) # dim of test

      for (i in 1:nrow(cv_MH))
      {
        a <- c(cv_MH$search[i], cv_MH$against[i])
        b <- c(size, size2) # dim of each df
        m <- matrix(c(a, b-a), ncol = 2)
        if (sum(chisq.test(m)$expected > 5) != 4)
        {
          count <- append(count, i)
        }
      }

      cv_MH <- cv_MH[-count, ]
      
      
      # calculate p-value
      cv_MH$pval <- prop.test(x = c(cv_MH$search, cv_MH$against), 
                              n = c(rep(size, length(cv_MH$search)),
                                    rep(size2, length(cv_MH$against))), correct = FALSE)$p.value
      # 
      # # proportions
      # cv_MH$search_prop <- cv_MH$search/size
      # cv_MH$against_prop <- cv_MH$against/size2
      
      
      # keep only significant MH terms
      cv_MH <- cv_MH[cv_MH$pval <= 0.05, ]

      # differences
      cv_MH$diffs <- cv_MH$search_prop - cv_MH$against_prop
      cv_MH$diffs <- ifelse(cv_MH$diffs == 0, 99999, cv_MH$diffs)
      
      # factors
      cv_MH$factor <- ifelse(cv_MH$search_prop > cv_MH$against_prop, cv_MH$search_prop/cv_MH$against_prop, -cv_MH$against_prop/cv_MH$search_prop)
      cv_MH$factor <- ifelse(cv_MH$factor == Inf, cv_MH$against, cv_MH$factor)
      cv_MH$factor <- ifelse(cv_MH$factor == -Inf, cv_MH$search, cv_MH$factor)
      
      # order by factor
      cv_MH <- cv_MH[order(-abs(cv_MH$factor)), ]
      
      # only keep terms requested as by user
      cv_MH <- cv_MH[1:input$MH_num, ]
      
      # plot graph by scaling diffs to range from -1 to 1
      cv_MH$scaled <- (((2)*(cv_MH$factor - min(cv_MH$factor)))/(max(cv_MH$factor) - min(cv_MH$factor))) - 1

      # heat map aesthetic
      cv_MH$V1 <- factor(cv_MH$V1, levels = (cv_MH$V1)[order(cv_MH$scaled)]) # in order
      
      # round factor to 2 digits
      cv_MH$factor <- lapply(cv_MH$factor, round, 2)
      cv_MH$factor <- as.numeric(cv_MH$factor)
      cv_MH$type <- ifelse(cv_MH$factor < 0, "against", "search")
      cv_MH$type <- as.factor(cv_MH$type)

      # plot heat map
      title_plot <- paste("Prevalence Factors of Top", input$MH_num, "Distinct MeSH terms in\n", input$phrases, "vs", input$compare_with)
      
      a <- ggplot(cv_MH, aes(x = " ", y = V1, fill = scaled)) + geom_tile() +
        xlab("") + ylab("MeSH Terms") + labs(title = title_plot) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1)) + labs(fill = "") +
        scale_fill_gradient(low = "#619CFF", high = "#F8766D",
                            breaks = c(-1, 0, 1), labels = c(input$compare_with, "0", input$phrases)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 12)) + 
        geom_text(aes(label = cv_MH$factor, color = cv_MH$type), show.legend = FALSE) +
        scale_colour_manual(values=c("white", "black"))

        
      return(a)
    })
    
    
    
    # do for RN terms
    df_RN <- reactive({
      # uses user's search
      result <- df()
      result2 <- df2()
      
      # cleans into a vector
      result <- as.data.frame(result[ , ])
      result_MH <- result[ , 3]
      names(result_MH)[1] <- "MH"
      
      # find occurrences of initial dataframe
      named_RN$search <- map_int(named_RN$V1, function(x){sum(str_detect(result_MH, fixed(x)))})
      
      
      # cleans into a vector
      result2 <- as.data.frame(result2[ , ])
      result_MH2 <- result2[ , 3]
      names(result_MH2)[1] <- "MH"
      
      # find occurrences of initial dataframe
      named_RN$against <- map_int(named_RN$V1, function(x){sum(str_detect(result_MH2, fixed(x)))})
      
      # delete all duplicate entries
      named_RN <- subset(named_RN, !duplicated(named_RN$V1))
      
      
      
      # size of both dataframes: find p-value for 2 prop z test
      size <- dim(result)[1]
      size2 <- dim(result2)[1]
      
      # remove chi squared warnings by finding and removing those that don't satisfy the conditions
      count <- numeric(dim(named_RN)[1]) # dim of test
      
      for (i in 1:nrow(named_RN))
      {
        a <- c(named_RN$search[i], named_RN$against[i])
        b <- c(size, size2) # dim of each df
        m <- matrix(c(a, b-a), ncol = 2)
        if (sum(chisq.test(m)$expected > 5) != 4)
        {
          count <- append(count, i)
        }
      }
      
      named_RN <- named_RN[-count, ]
      
      
      # calculate p-value
      named_RN$pval <- prop.test(x = c(named_RN$search, named_RN$against), 
                              n = c(rep(size, length(named_RN$search)),
                                    rep(size2, length(named_RN$against))), correct = FALSE)$p.value
      
      # proportions
      named_RN$search_prop <- named_RN$search/size
      named_RN$against_prop <- named_RN$against/size2

      # keep only significant MH terms
      named_RN <- named_RN[named_RN$pval <= 0.05, ]
      
      # differences
      named_RN$diffs <- named_RN$search_prop - named_RN$against_prop
      named_RN$diffs <- ifelse(named_RN$diffs == 0, 99999, named_RN$diffs)
      
      # factors
      named_RN$factor <- ifelse(named_RN$search_prop > named_RN$against_prop, named_RN$search_prop/named_RN$against_prop, -named_RN$against_prop/named_RN$search_prop)
      named_RN$factor <- ifelse(named_RN$factor == Inf, named_RN$against, named_RN$factor)
      named_RN$factor <- ifelse(named_RN$factor == -Inf, named_RN$search, named_RN$factor)
      
      
      # order by factor
      named_RN <- named_RN[order(-abs(named_RN$factor)), ]
      
      # only keep terms requested as by user
      named_RN <- named_RN[1:input$RN_num, ]
      
      # plot graph by scaling diffs to range from -1 to 1
      named_RN$scaled <- (((2)*(named_RN$factor - min(named_RN$factor)))/(max(named_RN$factor) - min(named_RN$factor))) - 1
      
      # heat map aesthetic
      named_RN$V1 <- factor(named_RN$V1, levels = (named_RN$V1)[order(named_RN$scaled)]) # in order
      
      # round factor to 2 digits
      named_RN$factor <- lapply(named_RN$factor, round, 2)
      named_RN$factor <- as.numeric(named_RN$factor)
      named_RN$type <- ifelse(named_RN$factor < 0, "against", "search")
      named_RN$type <- as.factor(named_RN$type)
      
      # plot heat map
      title_plot <- paste("Prevalence Factors of Top", input$RN_num, "Distinct RN terms in\n", input$phrases, "vs", input$compare_with)
      
      a <- ggplot(named_RN, aes(x = " ", y = V1, fill = scaled)) + geom_tile() +
        xlab("") + ylab("RN Terms") + labs(title = title_plot) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1)) + labs(fill = "") +
        scale_fill_gradient(low = "#619CFF", high = "#F8766D",
                            breaks = c(-1, 0, 1), labels = c(input$compare_with, "0", input$phrases)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 12)) + 
        geom_text(aes(label = named_RN$factor, color = named_RN$type), show.legend = FALSE) +
        scale_colour_manual(values=c("white", "black"))
      
      
      return(a)
    })
    
    dtm_MH_search <- reactive({
      test_sub <- df()
      
      # create corpus with phrases kept together based off https://stackoverflow.com/questions/24038498/corpus-build-with-phrases
      dat <- test_sub[ , 3]
      colnames(dat) <- c("text")
      
      # create 2 variables to combine into 1
      dat$docs <- "doc"
      dat$num <- ""
      dat$num <- 1:nrow(dat)
      
      # combine both variables
      dat$docs <- paste(dat$docs, dat$num, sep = "")
      dat <- dat[ , -c(3)]
      
      # dat <- dat[1:1000, ]
      
      x <- sub_holder(", ", dat$text)
      
      MH_parsed <- apply_as_tm(t(wfm(x$unhold(gsub(" ", "~~", x$output)), dat$docs)), 
                               weightTfIdf, to.qdap = FALSE)
      
      return(MH_parsed)
    })
    
    
    # pca plot on MH that user searched
    pca_MH_search <- reactive({
      MH_parsed <- dtm_MH_search()
      
      # PCA plot on MH and RN
      dtm_tf2 <- weightTfIdf(MH_parsed) # trying out RN_parsed
      m2 <- as.matrix(dtm_tf2)
      rownames(m2) <- 1:nrow(m2)
      norm_eucl <- function(m2) m2/apply(m2, MARGIN = 1, FUN = function(x) sum(x^2)^0.5) # normalize vectors so Euclidean distance makes sense
      m_norm2 <- norm_eucl(m2)
      m_norm2 <- m_norm2[, order(colSums(-m_norm2))]
      m_norm2 <- t(m_norm2)
      m_norm2[is.na(m_norm2)] <- 0
      m_norm2 <- m_norm2[order(-rowSums(m_norm2)), ] # orders them 
      m_norm2_sub <- m_norm2[1:100, ] # take only top 100
      pca <- prcomp((m_norm2_sub)) # temp changed from m_norm2_sub to see all values
      dat.loadings <- pca$x[ , 1:2]
      c <- kmeans(dat.loadings, centers = 10)
      pca1 <- pca$x[ , 1]
      pca2 <- pca$x[ , 2]
      mydf <- data.frame(ID = names(pca1), PCA1 = pca1, PCA2 = pca2, Cluster = factor(c$cluster))
      a <- which(mydf$ID == "doc")
      mydf <- mydf[-c(a), ]
      mydf2 <- mydf[1:input$pca_num, ]
      
      # plot pca plot on top 100 MH
      pca_title <- paste("The Top", input$pca_num, "MeSH terms within", input$phrases)
      
      d <- ggplot(mydf2, aes(x = PCA1, y = PCA2, label = ID, fill = Cluster)) + geom_point(aes(color = Cluster)) + 
        geom_label_repel(aes(label = ID)) + theme(legend.position = "none") + ggtitle(pca_title) +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) #+ geom_text()
    
      return(d)
      
    })
    
    
    dtm_MH_against <- reactive({
      test_sub <- df2()
      
      # create corpus with phrases kept together based off https://stackoverflow.com/questions/24038498/corpus-build-with-phrases
      dat <- test_sub[ , 3]
      colnames(dat) <- c("text")
      
      # create 2 variables to combine into 1
      dat$docs <- "doc"
      dat$num <- ""
      dat$num <- 1:nrow(dat)
      
      # combine both variables
      dat$docs <- paste(dat$docs, dat$num, sep = "")
      dat <- dat[ , -c(3)]
      
      x <- sub_holder(", ", dat$text)
      
      MH_parsed <- apply_as_tm(t(wfm(x$unhold(gsub(" ", "~~", x$output)), dat$docs)), 
                               weightTfIdf, to.qdap = FALSE)
      
      return(MH_parsed)
    })
    
    # pca plot on MH terms user did compare with
    pca_MH_against <- reactive({
      MH_parsed <- dtm_MH_against()
      
      # PCA plot on MH and RN
      dtm_tf2 <- weightTfIdf(MH_parsed) # trying out RN_parsed
      m2 <- as.matrix(dtm_tf2)
      rownames(m2) <- 1:nrow(m2)
      norm_eucl <- function(m2) m2/apply(m2, MARGIN = 1, FUN = function(x) sum(x^2)^0.5) # normalize vectors so Euclidean distance makes sense
      m_norm2 <- norm_eucl(m2)
      m_norm2 <- m_norm2[, order(colSums(-m_norm2))]
      m_norm2 <- t(m_norm2)
      m_norm2[is.na(m_norm2)] <- 0
      m_norm2 <- m_norm2[order(-rowSums(m_norm2)), ] # orders them 
      m_norm2_sub <- m_norm2[1:100, ] # take only top 100
      pca <- prcomp((m_norm2_sub)) # temp changed from m_norm2_sub to see all values
      dat.loadings <- pca$x[ , 1:2]
      c <- kmeans(dat.loadings, centers = 10)
      pca1 <- pca$x[ , 1]
      pca2 <- pca$x[ , 2]
      mydf <- data.frame(ID = names(pca1), PCA1 = pca1, PCA2 = pca2, Cluster = factor(c$cluster))
      a <- which(mydf$ID == "doc")
      mydf <- mydf[-c(a), ]
      mydf2 <- mydf[1:input$pca_num_against, ]
      
      # plot pca plot on top 100 MH
      pca_title <- paste("The Top", input$pca_num, "MeSH terms within", input$compare_with)
      
      d <- ggplot(mydf2, aes(x = PCA1, y = PCA2, label = ID, fill = Cluster)) + geom_point(aes(color = Cluster)) + 
        geom_label_repel(aes(label = ID)) + theme(legend.position = "none") + ggtitle(pca_title) +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) #+ geom_text()
      
      return(d)
      
    })
    
    
    # make correlation plot on MH search by user
    corr_MH_search <- reactive({
      final <- df()
      
      # creating corpus on variable that I want to create plot on
      myCorpus <- Corpus(VectorSource(final$MH2))  
      dtm2 <- DocumentTermMatrix(myCorpus)
      
      # correlation of terms plot
      title_search <- paste("Correlation Plot among Top", input$corr_search, "MeSH terms within", input$phrases)

      freq.terms <- findFreqTerms(dtm2)[1:input$corr_search] # choose top 25 terms
      f <- plot(dtm2, term = freq.terms, corThreshold = 0.1, weighting = T, main = title_search,
                edgeAttrs = eAttrs) # choose terms with correlation of at least 0.1
      
      return(f)
    })
    
    
    # make correlation plot on MH against by user
    corr_MH_against <- reactive({
      final <- df2()
      
      # creating corpus on variable that I want to create plot on
      myCorpus <- Corpus(VectorSource(final$MH2))  
      dtm2 <- DocumentTermMatrix(myCorpus)
      
      # correlation of terms plot
      title_against <- paste("Correlation Plot among Top", input$corr_against, "MeSH terms within", input$compare_with)
      
      freq.terms <- findFreqTerms(dtm2)[1:input$corr_against] # choose top 25 terms
      g <- plot(dtm2, term = freq.terms, corThreshold = 0.1, weighting = T, main = title_against,
                edgeAttrs = eAttrs) # choose terms with correlation of at least 0.1
      
      return(g)
    })
    
    output$MH <- renderPlot({df_MH()})
    output$RN <- renderPlot({df_RN()})
    output$pca_MH_search <- renderPlot({pca_MH_search()})
    output$pca_MH_against <- renderPlot({pca_MH_against()})
    output$corr_MH_search <- renderPlot({corr_MH_search()})
    output$corr_MH_against <- renderPlot({corr_MH_against()})
    output$jaccard <- renderText({jaccard()})
    output$colorbox <- renderPlot({
      jaccard_num <- jaccard_num()
      jaccard_num <- as.numeric(jaccard_num)
      # get an empty box
      plot(0, 0, xlim=c(0,1), ylim=c(0,1))
      plot(0:10, type = "n", axes = FALSE, xlab = NA, ylab = NA)
      # rectangle filled with a gradient
      gradient.rect(0, 0, 10, 3, col = smoothColors("#F8766D", 99, "#619CFF"), border = NA)
      # vertical bar 
      segments(jaccard_num*10, 0, jaccard_num*10, 3, lwd = 3)
      text(jaccard_num*10, 0, jaccard_num, pos = 1, xpd = TRUE)
      text(1, 4, "dissimilar search", pos = 1, xpd = TRUE)
      text(10, 4, "identical search", pos = 1, xpd = TRUE)
    })
  output$treemap_1 <- renderPlot({treemap_1()})
  output$treemap_2 <- renderPlot({treemap_2()})
}
