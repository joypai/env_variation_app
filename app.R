library(shiny)
library(shinythemes)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(plotly)
library(ggseqlogo)
library(ggsci)
library(stringr)
library(reshape2)


edit_binding_sites <- function(ab_binding_pos_orig, hxb2) {
  ab_binding_pos <- c()
  gaps <- 0
  seq <- as.character(hxb2)
  
  for (pos in ab_binding_pos_orig) {
    gaps <- str_count(substr(seq, 1, pos+gaps), "-")
    ab_binding_pos <- c(ab_binding_pos, pos+gaps)
  }
  return(ab_binding_pos)
}

calc_variation <- function(seqs_aa, consensus_aa) {
  mismatches <- sum(sapply(seqs_aa, function(x) x!=consensus_aa))
  return(mismatches/length(seqs_aa))
}

AA_breakdown <- function(seqs, timepoint) {
  all_pos_df <- data.frame(pos=numeric(0), timepoint=character(0), AA=character(0), count=numeric(0), frequency=numeric(0))
  for (i in 1:width(seqs)[1]) {
    AA_column <- subseq(seqs, i, i)
    AA_df <- as.data.frame(table(as.character(AA_column)))
    colnames(AA_df) <- c("AA", "count")
    AA_df$frequency <- AA_df$count/length(seqs)
    all_pos_df <- rbind(all_pos_df, cbind(pos=i, timepoint=timepoint, AA_df))
  }
  
  return(all_pos_df)
}

find_new_mutations <- function(AA_by_pos) {
  new_mut_pos <- data.frame(pos=numeric(0), type=character(0))
  
  for (i in 1:max(AA_by_pos$pos)) {
    pos_table <- subset(AA_by_pos, pos == i)
    pre_aas <- subset(pos_table, timepoint=="W23")$AA
    rebound_aas <- subset(pos_table, timepoint=="rebound")$AA
    if (length(setdiff(rebound_aas, pre_aas)) > 0) {
      #type <- ifelse(subset(pos_table, timepoint=="rebound" & AA == setdiff(rebound_aas, pre_aas)[1])$frequency > 0.5, "major","minor")
      frequency <- subset(pos_table, timepoint=="rebound" & AA == setdiff(rebound_aas, pre_aas)[1])$frequency
      
      new_mut_pos <- rbind(new_mut_pos, data.frame(pos=i, frequency=frequency))
    }
  }
  
  return(new_mut_pos)
}



ui <- shinyUI(fluidPage(
      theme = shinytheme("cosmo"),
      title = "HIV env variation analysis",
      #shinythemes::themeSelector(),
      
      titlePanel("HIV env variation analysis"),
      hr(),
      
      # sidebar with parameters
      fluidRow(column(width = 1,
                      wellPanel(
                        selectInput(
                          "patient",
                          "Patient:",
                          as.character(c(601, 602, 603, 605, 608, 609, 610, 611, 613))
                        )
                      )),
               
               # main plot
               mainPanel(width=11,#tabsetPanel(
                 #tabPanel("main",
                          fluidRow(
                            column(width = 12, h4("variation plot"),
                                   plotlyOutput("variation_plot", height=350)
                            ), 
                          fluidRow(
                            column(width=12, h4("positions with new mutations in rebound"),
                                   plotlyOutput("new_muts_plot", height=105)),
                            code("Click on position to display AA logo plot"))
                          
                          ),
                          
                          br(),
                          fluidRow(
                            column(width=1),
                            column(width=3, h4("W23 amino acids"), tableOutput("pre_aa_table")),
                            column(width=3, h4("rebound amino acids"), tableOutput("rebound_aa_table")),
                            column(width=4, h4("AA logo plot"), plotOutput("react", height=200, width=300)),
                            column(width=1)
                          )
                          )
                 #tabPanel("new mutations")
                        #  ))
)))



# initial computation of necessary information -----------
ab_binding_pos_orig <- as.numeric(read.csv("~/data_alpha/home/jpai/julio/escape_positions_analysis/3BNC117_pos.txt", header=F))
patients <- as.character(c(601, 602, 603, 605, 608, 609, 610, 611, 613))
rebound_times <- list("W29","W28","W31","W29","W29","W29","W31","W31","W37")
names(rebound_times) = patients

variation_df_list <- list()
ab_binding_pos_list <- list()
AA_by_pos_list <- list()
pos_less_variation_list <- list()
new_mut_pos_list <- list()

# analyze for each patient
for (patient in patients) {
  fasta_file = paste("~/data_alpha/home/jpai/julio/variation/",patient,"/",patient,"aacodonalign.fasta",sep="")
  fasta <- readAAStringSet(fasta_file)
  
  rebound_wk <- rebound_times[[patient]]
  ab_binding_pos <- edit_binding_sites(ab_binding_pos_orig, fasta["HXB2_ENV_K03455.1"])
  ab_binding_pos_list[[patient]] <- ab_binding_pos
  
  # split sequences into pre and rebound
  seq_ids <- names(fasta)
  W23_seq_ids <- c()
  if (patient=="609"){
    W23_seq_ids <- seq_ids[grepl("pre", seq_ids)]
  } else {
    W23_seq_ids <- seq_ids[grepl("W23", seq_ids)]
  }
  rebound_seq_ids <- seq_ids[grepl(rebound_wk, seq_ids)]
  W23_seqs <- fasta[W23_seq_ids]
  rebound_seqs <- fasta[rebound_seq_ids]
  
  # get AA by position
  AA_by_pos <- AA_breakdown(W23_seqs, "W23")
  AA_by_pos <- rbind(AA_by_pos, AA_breakdown(rebound_seqs, "rebound"))
  AA_by_pos_list[[patient]] <- AA_by_pos
  
  new_mut_pos <- find_new_mutations(AA_by_pos)
  new_mut_pos_list[[patient]] <- new_mut_pos
  
  # generate consensus sequence from pre (W23) timepoint
  tmp_cons <- consensusString(W23_seqs)
  m <- gregexpr("\\?", tmp_cons)
  if (m[[1]][1] != -1) {
    pre_consensus <- ""
    for (i in m[[1]]) {
      aas <- table(as.character(subseq(W23_seqs, i, i)))
      mc <- names(aas[which.max(aas)])
      pre_consensus <- paste(pre_consensus, substr(tmp_cons, nchar(pre_consensus)+1,i-1), mc, sep="")
    }
  } else {
    pre_consensus <- tmp_cons
  }
  pre_consensus <- paste(pre_consensus, substr(tmp_cons, nchar(pre_consensus)+1, nchar(tmp_cons)), sep="")

  # calculate variation based on W23 consensus
  pre_df <- data.frame(do.call(rbind, strsplit(as.character(W23_seqs), ''))) 
  rownames(pre_df) <- W23_seq_ids

  rebound_df <- data.frame(do.call(rbind, strsplit(as.character(rebound_seqs), ''))) 
  rownames(rebound_df) <- rebound_seq_ids

  var_df <- data.frame(pos=seq(1, nchar(pre_consensus)))
  var_df$pre <- sapply(var_df$pos, function(x) calc_variation(pre_df[,x], substr(pre_consensus, x, x)))
  var_df$rebound <- sapply(var_df$pos, function(x) calc_variation(rebound_df[,x], substr(pre_consensus, x, x)))
  #var_df$diff <- ifelse(var_df$pre == 0, var_df$rebound, 0)
  #var_df$less_var <- with(var_df, ifelse(var_df$pre > var_df$rebound, var_df$rebound, NA))
  #var_df$more_var <- with(var_df, ifelse(var_df$pre < var_df$rebound, var_df$rebound, NA))
  
  
  pos_less_variation <- as.numeric(rownames(subset(var_df, pre > rebound)))
  pos_less_variation_list[[patient]] <- pos_less_variation
  #less_var_df <- data.frame(xmin=pos_less_variation-0.5, xmax=pos_less_variation+0.5, ymin=rep(-.05, length(pos_less_variation)), ymax= rep(0.6,length(pos_less_variation)))
  variation_df_list[[patient]] <- var_df
}


# server -----------
server <- shinyServer(function(input, output) {
  
  output$variation_plot <- renderPlotly({ 
    var_df <- variation_df_list[[input$patient]]
    ab_binding_pos <- ab_binding_pos_list[[input$patient]]
    pos_less_variation <- pos_less_variation_list[[input$patient]]
    
    print(head(var_df))
    var_melt_df <- melt(var_df, id.vars = "pos", measure.vars = c("pre", "rebound"), variable.name = "timepoint",value.name = "variation")
    
    #var_melt_df$sel_type <- with(var_melt_df, ifelse(var_melt_df$pos %in% pos_less_variation,"less","none"))
    
    p1 <- ggplot(var_melt_df, aes(x=pos, y=variation, color=timepoint)) + geom_line(size=0.4) + facet_wrap(~timepoint, nrow=2) + theme_bw() +
      #geom_point(data=subset(var_melt_df, sel_type=="less" & timepoint=="rebound"), aes(x=pos, y=variation, color=timepoint), size=1) +
      geom_point(data=data.frame(x=ab_binding_pos, y=-0.1), aes(x=x,y=y), shape=15, size=0.5, color="black") + 
      scale_color_aaas() + coord_cartesian(xlim=c(0, 900)) + ylab("variation (compared to W23 consensus)") +
      annotate("text", x=50, y=-0.1, label="3BNC117 binding sites", size=3) +
      #geom_rect(data=data.frame(xmin=ab_binding_pos-0.5, xmax=ab_binding_pos+0.5, ymin=-0.1, ymax=1.0), 
      #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="grey", alpha=0.15, inherit.aes = FALSE) +
      guides(color=FALSE, alpha=TRUE)
  
    ggplotly(p1) %>% layout(margin=list(l = 90), showlegend=FALSE) 
  
  })
  
  output$new_muts_plot <- renderPlotly({
    new_mut_pos <- new_mut_pos_list[[input$patient]]
    new_mut_pos$y <- 1
    
    
    p2 <- ggplot(new_mut_pos, aes(x=pos,y=y, color=frequency)) + geom_point() + theme_bw() + coord_cartesian(xlim=c(0, 900)) + 
      theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + guides(color=FALSE) +
      labs(x = "pos", y = "") + scale_color_material("green", limits=c(0,1)) #+ scale_color_simpsons()
    
    ggplotly(p2, tooltip = c("pos", "frequency")) %>% layout(margin=list(l = 90), showlegend=FALSE)# legend = list(orientation = "h"))
    
  })
  #margin=list(l = 120), 
  output$pre_aa_table <- renderTable({
    d <- event_data("plotly_click")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      selected_pos <- d$x
      AA_selected_pos <- subset(AA_by_pos, pos==selected_pos & timepoint=="W23")
      AA_selected_pos[with(AA_selected_pos, order(-frequency)),]
    }
  })
  
  output$rebound_aa_table <- renderTable({
    d <- event_data("plotly_click")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      selected_pos <- d$x
      AA_selected_pos <- subset(AA_by_pos, pos==selected_pos & timepoint=="rebound")
      AA_selected_pos[with(AA_selected_pos, order(-frequency)),]
    }
  })
  
  
  output$react <- renderPlot ({
    d <- event_data("plotly_click")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      print(d)
      selected_pos <- d$x
      seq_logos = list()
      
      AA_selected_pos <- subset(AA_by_pos, pos==selected_pos)
      a <- subset(AA_selected_pos, timepoint=="W23")
      b <- as.character(rep(a$AA, a$count))
      names(b) <- rep(1:length(b))
      seq_logos[['pre']] <- b
      a <- subset(AA_selected_pos, timepoint=="rebound")
      b <- as.character(rep(a$AA, a$count))
      names(b) <- rep(1:length(b))
      seq_logos[['rebound']] <- b
      
      if ("-" %in% unlist(seq_logos)) {
        a <- as.relistable(seq_logos)
        u <- unlist(a)
        u[u == "-"] <- "x"
        seq_logos <- relist(u, a)
        ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="other", namespace=c(names(AMINO_ACID_CODE), 'x')) 
      } else {
        ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="aa") 
      }
    }
  })
})



# Run the application 
shinyApp(ui = ui, server = server)
