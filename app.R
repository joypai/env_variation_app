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
library(gridExtra)
library(reshape2)

patients <- as.character(c(601:603, 605, 608:611, 613))

# helper functions ---------
col_to_hxb2_pos <- function(hxb2) {
  hxb2_pos_mapping <- data.frame(col=numeric(0), hxb2_pos=character(0))
  prev <- 1
  gaps <- 0
  
  for (x in 1:width(hxb2)) {
    if (substr(hxb2, x, x) == "-" ) { # gap
      new_row <- data.frame(col=x, hxb2_pos=paste(prev-1,letters[gaps+1], sep=""))
      gaps <- gaps + 1
    } else{
      gaps <- 0
      new_row <- data.frame(col=x, hxb2_pos=as.character(prev))
      prev <- prev + 1
    }
    hxb2_pos_mapping <- rbind(hxb2_pos_mapping, new_row)
  }
  return(hxb2_pos_mapping)
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
  check_muts <- function(i) {
    pos_table <- subset(AA_by_pos, pos == i)
    pre_aas <- subset(pos_table, timepoint=="W23")$AA
    rebound_aas <- subset(pos_table, timepoint=="rebound")$AA
    if (length(setdiff(rebound_aas, pre_aas)) > 0) {
      frequency <- subset(pos_table, timepoint=="rebound" & AA == setdiff(rebound_aas, pre_aas)[1])$frequency
      new_row <- data.frame(pos=i, frequency=frequency)
      return(new_row)
    }
  }
  
  new_mut_pos <- do.call(rbind,lapply(1:max(AA_by_pos$pos),function(x) check_muts(x)))
  return(new_mut_pos)
}

make_variation_plot <- function(var_df, ab_binding_pos, src, ab_label) {
  var_melt_df <- melt(var_df, id.vars = c("pos","hxb2_pos"), measure.vars = c("pre", "rebound"), variable.name = "timepoint",value.name = "variation")
  
  p1 <- ggplot(var_melt_df, aes(x=pos, y=variation, color=timepoint, label=hxb2_pos)) + geom_line(size=0.4) + facet_wrap(~timepoint, nrow=2) + theme_bw() +
    geom_point(data=data.frame(x=ab_binding_pos, y=-0.1, hxb2_pos=names(ab_binding_pos)), aes(x=x,y=y), shape=15, size=0.5, color="black") + 
    scale_color_aaas() + coord_cartesian(xlim=c(0, 900)) + ylab("variation (compared to W23 consensus)") +
    annotate("text", x=50, y=-0.1, label=ab_label, size=3) +
    guides(color=FALSE, alpha=TRUE)
  
  ggplotly(p1, tooltip=c("pos","hxb2_pos","variation"), source=src) %>% layout(margin=list(l = 90), showlegend=FALSE) 
}

plot_new_mutations <- function(new_mut_pos, src) {
  new_mut_pos$y <- 1
  p2 <- ggplot(new_mut_pos, aes(x=pos,y=y, color=frequency, label=hxb2_pos)) + geom_point() + theme_bw() + coord_cartesian(xlim=c(0, 900)) + 
    theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + guides(color=FALSE) +
    labs(x = "pos", y = "") + scale_color_material("green", limits=c(0,1)) #+ scale_color_simpsons()
  
  ggplotly(p2, tooltip = c("pos","hxb2_pos", "frequency"), source=src) %>% layout(margin=list(l = 90), showlegend=FALSE)
}

compute_aa_table <- function(tp, clicked_pos) {
  tables <- list()
  
  for (patient in patients) {
    selected_pos <- ab_binding_pos_list[[patient]][clicked_pos]
    AA_by_pos <- AA_by_pos_list[[patient]]
    AA_selected_pos <- subset(AA_by_pos, pos==selected_pos & timepoint==tp)
    
    tables[[patient]] <- print(xtable(AA_selected_pos[with(AA_selected_pos, order(timepoint,-frequency)),],
                                      caption=paste("patient:", patient)),
                               type='html', caption.placement="top", include.rownames = FALSE,
                               html.table.attributes='class="data table table-bordered table-condensed"')
  }
  
  return(lapply(tables, paste))
}

create_aa_table <- function(AA_by_pos, selected_pos, tp) {
  AA_selected_pos <- subset(AA_by_pos, pos==selected_pos & timepoint==tp)
  AA_selected_pos[with(AA_selected_pos, order(-frequency)),]
}

make_logo_plot <- function(AA_by_pos, selected_pos) {
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
  
  # make logo plot
  if ("-" %in% unlist(seq_logos)) {   # convert gaps to "x" so seqs are compatible with ggseqlogo
    a <- as.relistable(seq_logos)
    u <- unlist(a)
    u[u == "-"] <- "x"
    seq_logos <- relist(u, a)
    ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="other", namespace=c(names(AMINO_ACID_CODE), 'x')) 
  } else {
    ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="aa")
  }
}

analyze_variation <- function(patient, fasta_file, ab_binding_pos_orig, rebound_wk, pre_label) {
  results <- list()
  fasta <- readAAStringSet(fasta_file)
  
  #rebound_wk <- rebound_times[[patient]]
  hxb2_pos_mapping <- col_to_hxb2_pos(fasta[1])
  ab_binding_pos <- as.numeric(sapply(ab_binding_pos_orig, function(x) subset(hxb2_pos_mapping, hxb2_pos==x)$col))
  names(ab_binding_pos) <- ab_binding_pos_orig
  results[['ab_pos']] <- ab_binding_pos
  
  # split sequences into pre and rebound
  names(fasta) <- gsub("-","_",names(fasta))
  seq_ids <- names(fasta)
  W23_seq_ids <- c()
  if (patient=="609"){
    W23_seq_ids <- seq_ids[grepl("pre", seq_ids)]
  } else {
    W23_seq_ids <- seq_ids[grepl(paste("_",pre_label,"_",sep=""), seq_ids)]
  }
  W23_seqs <- fasta[W23_seq_ids]
  results[['pre']] <- W23_seqs
  
  rebound_seq_ids <- seq_ids[grepl(paste("_",rebound_wk,"_",sep=""), seq_ids)]
  rebound_seqs <- fasta[rebound_seq_ids]
  results[['rebound']] <- rebound_seqs
  
  # get AA by position
  AA_by_pos <- AA_breakdown(W23_seqs, "W23")
  AA_by_pos <- rbind(AA_by_pos, AA_breakdown(rebound_seqs, "rebound"))
  results[['AA']] <- AA_by_pos
  
  new_mut_pos <- find_new_mutations(AA_by_pos)
  new_mut_pos$hxb2_pos <- as.character(sapply(new_mut_pos$pos, function(x) subset(hxb2_pos_mapping, col==x)$hxb2_pos))
  results[['new_mut']] <- new_mut_pos
  
  # generate consensus sequence from pre (W23) timepoint
  tmp_cons <- consensusString(W23_seqs)
  m <- gregexpr("\\?", tmp_cons)    # replace ambiguities from biostring-generated consensus with majority call
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
  var_df <- cbind(var_df, hxb2_pos_mapping[,"hxb2_pos", drop=FALSE])
  var_df$pre <- sapply(var_df$pos, function(x) calc_variation(pre_df[,x], substr(pre_consensus, x, x)))
  var_df$rebound <- sapply(var_df$pos, function(x) calc_variation(rebound_df[,x], substr(pre_consensus, x, x)))
  results[['var']] <- var_df
  
  return(results)
}


# UI -----------
ui <- shinyUI(
  fluidPage(
    theme = shinytheme("cosmo"),
    title = "HIV env variation analysis",
    
    titlePanel("HIV env variation analysis"),
    hr(),
    
    # sidebar with parameters and inputs
    fluidRow(
      column(width = 2,
         wellPanel(
           selectInput("patient", "Choose patient:", patients),
           hr(),
           fileInput("new_fasta", "Upload new file:", accept=".fasta"), 
           textInput("new_pre_label", "pre"),
           textInput("new_rebound_label", "rebound"),
           checkboxInput("ab", "include 10-1074 positions", value=FALSE),
           actionButton("confirm","GO"),
           helpText("file must be in fasta format",br(),br(),
                    "pre: pattern common to all sequence names from timepoint before rebound", br(), br(),
                    "rebound: pattern common to all sequence names from rebound timepoint")
         )
      ),
      
      # main plot
      mainPanel(width = 10,
        tabsetPanel("tabs",
          tabPanel("variation",
           fluidRow(
             column(width = 12,
                    h4("variation plot"),
                    plotlyOutput("variation_plot", height = 350)
           ),
           fluidRow(
             column(width = 12,
                    h4("positions with new mutations in rebound"),
                    plotlyOutput("new_muts_plot", height = 125)
             ),
             code("Click on position to display AA logo plot")
           )),
           
           br(),
           fluidRow(
             column(width = 1),
             column(width = 3, 
                    h4("W23 amino acids"),
                    tableOutput("pre_aa_table")),
             column(width = 3,
                    h4("rebound amino acids"),
                    tableOutput("rebound_aa_table")
             ),
             column(width = 4,
                    h4("AA logo plot"),
                    plotOutput("logo", height = 200, width = 300)
             )
           )
          ),
          tabPanel("3BNC117 binding sites",
           fluidRow(
             column(width = 6, 
                    h4("variation and binding sites"), 
                    plotlyOutput("binding_pos_all_plot", height=1200)),
             column(width = 3, 
                    h4(textOutput("binding_site")),
                    uiOutput("W23_aa_all_table")),
             column(width = 3, 
                    br(), br(), 
                    uiOutput("rebound_aa_all_table"))
           )
          ),
          tabPanel("new data",
           fluidRow(
             column(width=12, 
                    h4("variation plot"),
                    plotlyOutput("variation_plot_nd", height = 350))
           ),
           fluidRow(
             column(width = 12,
                    h4("positions with new mutations in rebound"),
                    plotlyOutput("new_muts_plot_nd", height = 125)
             )
           ),
           br(),
           fluidRow(
             column(width = 1),
             column(width = 3,
                    h4("W23 amino acids"),
                    tableOutput("pre_aa_table_nd")),
             column(width = 3,
                    h4("rebound amino acids"),
                    tableOutput("rebound_aa_table_nd")
             ),
             column(width = 4,
                    h4("AA logo plot"),
                    plotOutput("logo_nd", height = 200, width = 300)
             )
           )
          )
        )
      )
    )
  )
)

# initial computation of necessary information -----------
ab_binding_pos_3BNC <- as.numeric(read.csv("~/data_alpha/home/jpai/julio/escape_positions_analysis/3BNC117_pos.txt", header=F))
ab_binding_pos_101074 <- as.numeric(read.csv("~/data_alpha/home/jpai/julio/escape_positions_analysis/3BNC117_101074_pos.txt", header=F))

rebound_times <- list("W29","W28","W31","W29","W29","W29","W31","W31","W37")
names(rebound_times) = patients

variation_df_list <- list()
ab_binding_pos_list <- list()
AA_by_pos_list <- list()
new_mut_pos_list <- list()
W23_seqs_list <- list()
rebound_seqs_list <- list()

# analyze for each patient
for (patient in patients) {
  rebound_wk <- rebound_times[[patient]]
  fasta_file = paste("~/data_alpha/home/jpai/julio/variation/variation_app/data/",patient,"aacodonalign.fasta",sep="")
  
  results <- analyze_variation(patient, fasta_file, ab_binding_pos_3BNC, rebound_wk, "W23")
  
  ab_binding_pos_list[[patient]] <- results[['ab_pos']]
  W23_seqs_list[[patient]] <- results[['pre']]
  rebound_seqs_list[[patient]] <- results[['rebound']]
  AA_by_pos_list[[patient]] <- results[['AA']]
  new_mut_pos_list[[patient]] <- results[['new_mut']]
  variation_df_list[[patient]] <- results[['var']]
}


# server -----------
server <- shinyServer(function(input, output, session) {
  
  output$variation_plot <- renderPlotly({ 
    var_df <- variation_df_list[[input$patient]]
    ab_binding_pos <- ab_binding_pos_list[[input$patient]]
    make_variation_plot(var_df, ab_binding_pos, "old_data", "3BNC117 binding sites")
    
  })
  
  analyze_new <- eventReactive(input$confirm, {
    fasta_file = input$new_fasta$datapath
    patient <- strsplit(input$new_fasta$name, "_")[[1]][1]
    print(patient)
    print(paste("pre:",input$new_pre_label))
    print(paste("rebound:",input$new_rebound_label))
    
    if(input$ab) {
      ab_binding_pos_orig <- ab_binding_pos_101074
    } else {
      ab_binding_pos_orig <- ab_binding_pos_3BNC
    }
    results <- analyze_variation(patient, fasta_file, ab_binding_pos_orig, input$new_rebound_label, input$new_pre_label)
    
    return(results)
  })
  
  output$variation_plot_nd <- renderPlotly({
    req(input$new_fasta)
    results <- analyze_new()
    make_variation_plot(results$var, results$ab, "new_data", "3BNC117 & 10-1074 binding sites")
  })
  
  output$new_muts_plot <- renderPlotly({
    new_mut_pos <- new_mut_pos_list[[input$patient]]
    plot_new_mutations(new_mut_pos, "old_data")
  })
  
  output$new_muts_plot_nd <- renderPlotly({
    results <- analyze_new()
    plot_new_mutations(results$new_mut, "new_data")
  })
  
  output$binding_pos_all_plot <- renderPlotly ({
    bs_var_all_df <- data.frame()
    
    for (patient in patients) {
      var_df <- variation_df_list[[patient]]
      ab_binding_pos <- ab_binding_pos_list[[patient]]
      bs_var_df <- subset(var_df, pos %in% ab_binding_pos)
      bs_var_df$pos <- as.character(bs_var_df$pos)
      bs_var_df$patient <- patient
      bs_var_all_df <- rbind(bs_var_all_df, bs_var_df)
    }
    
    var_melt_df <- melt(bs_var_all_df, id.vars = c("pos","hxb2_pos","patient"), measure.vars = c("pre", "rebound"), variable.name = "timepoint",value.name = "variation")
    
    p1 <- ggplot(var_melt_df, aes(x=hxb2_pos, y=variation, color=timepoint, label=pos)) + geom_point(size=1) + theme_classic() +
      scale_color_aaas() + xlab("3BNC117 binding site") + ylab("variation (compared to W23 consensus)") + facet_wrap(~patient, ncol=1)
    
    ggplotly(p1, tooltip=c("pos","hxb2_pos","variation"), source="binding_sites") %>% layout(margin=list(l = 90), legend = list(orientation = "h")) 
    
  })
  
  output$pre_aa_table <- renderTable({
    d <- event_data("plotly_click", source="old_data")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      create_aa_table(AA_by_pos, d$x, "W23")
    }
  })
  
  output$rebound_aa_table <- renderTable({
    d <- event_data("plotly_click", source="old_data")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      create_aa_table(AA_by_pos, d$x, "rebound")
    }
  })
  
  output$pre_aa_table_nd <- renderTable({
    d <- event_data("plotly_click", source="new_data")
    
    if (!is.null(d)){
      results <- analyze_new()
      create_aa_table(results$AA, d$x, "W23")
    }
  })
  
  output$rebound_aa_table_nd <- renderTable({
    d <- event_data("plotly_click", source="new_data")
    
    if (!is.null(d)){
      results <- analyze_new()
      create_aa_table(results$AA, d$x, "rebound")
    }
  })
  
  output$W23_aa_all_table <- renderUI({
    d <- event_data("plotly_click", source="binding_sites")
    
    if (!is.null(d)){
      tables_out <- paste(compute_aa_table("W23", d$x), collapse='')
      return(div(HTML(tables_out), class="shiny-html-output"))
    }
  })
  
  output$rebound_aa_all_table <- renderUI({
    d <- event_data("plotly_click", source="binding_sites")
    
    if (!is.null(d)){
      tables_out <- paste(compute_aa_table("rebound", d$x), collapse='')
      return(div(HTML(tables_out), class="shiny-html-output"))
    }
  })
  
  output$binding_site <- renderText({
    d <- event_data("plotly_click", source="binding_sites")
    if (!is.null(d)){
      return(paste("binding site:", ab_binding_pos_3BNC[d$x]))
    }
  })
  
  output$logo <- renderPlot ({
    d <- event_data("plotly_click", source="old_data")
    
    if (!is.null(d)){
      AA_by_pos <- AA_by_pos_list[[input$patient]]
      make_logo_plot(AA_by_pos, d$x)
    }
  })
  
  output$logo_nd <- renderPlot ({
    d <- event_data("plotly_click", source="new_data")
    
    if (!is.null(d)){
      results <- analyze_new()
      make_logo_plot(results$AA, d$x)
    }
  })
  
  output$logo_plots <- renderPlot ({
    logo_plots = list()
    
    for (patient in patients) {
      AA_by_pos <- AA_by_pos_list[[patient]]
      ab_binding_pos <- ab_binding_pos_list[[patient]]
      
      W23_seqs_mat <- as.matrix(W23_seqs_list[[patient]])[,ab_binding_pos]
      rebound_seqs_mat <- as.matrix(rebound_seqs_list[[patient]])[,ab_binding_pos]
      
      seq_logos = list()
      seq_logos[['pre']] <- apply(W23_seqs_mat, 1, paste, collapse="")
      seq_logos[['rebound']] <- apply(rebound_seqs_mat, 1, paste, collapse="")
      pos_labels <- as.character(ab_binding_pos_orig)
      names(pos_labels) <- as.character(c(1:31))
      
      # make logo plot
      if ("-" %in% unlist(seq_logos)) {   # convert gaps to "x" so seqs are compatible with ggseqlogo
        a <- as.relistable(seq_logos)
        u <- unlist(a)
        u[u == "-"] <- "x"
        seq_logos <- relist(u, a)
        p <- ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="other", namespace=c(names(AMINO_ACID_CODE), 'x')) + ggtitle(patient) +
          scale_x_discrete(labels=pos_labels)
      } else {
        p <- ggseqlogo(seq_logos, method='probability', font="helvetica_regular", seq_type="aa") + ggtitle(patient) +
          scale_x_discrete(labels=pos_labels)
      }
      logo_plots[[patient]] <- p
    }
    
    do.call(gridExtra::grid.arrange, c(logo_plots, ncol=1))
    
  })
  
})

# run the application -----------
shinyApp(ui = ui, server = server)