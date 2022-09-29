## script for dnds figure
library("ggplot2")
library("ggthemes")
library("gridExtra")
library("extrafont")
font_import()
# Distribution of sequence variants in alternative and principal isoforms

type_principal <- c("PI", "NPI")
for (principal in type_principal){
  path <- "/functional_isoforms/paper_set/trifid/"
  path = paste0(path,principal, '/')
  print(path)
  setwd(path)
  files <- list.files(pattern = ".txt") # get all final files from directory
  sum_variations_all <- list()
  for (f in 1:length(files)){ #read tables and add sum of variation columns in a vector
    npi_functional <- read.table(files[f], row.names = 1, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    npi_functional_sum <- npi_functional[,6:13]
    sum_variations <- apply(npi_functional_sum,2,sum)
    sum_variations_all[[f]] = sum_variations # add vector to list
  }
  rm(npi_functional, npi_functional_sum)
  # we first get the file names (NPI, PI)
  # then we get their functionality (functional, nonfunctional)
  names(sum_variations_all) = paste0(sapply(strsplit(files, "_"),'[', 1), "_",
                                     sapply(strsplit(sapply(strsplit(files, ".txt"),
                                                            '[', 1), "_"), '[', 3))
  
  ## get dnds for both rare and common variants in alternative isoforms:
  n_npi <- c("_01", "_02",  "_03", "_04", "_05","_06", "_07",  "_08", "_09", "_1")
  n_npi <- paste0(principal, n_npi)
  rare_error_min <- list()
  rare_error_max <- list()
  common_error_min <- list()
  common_error_max <- list()
  rare_dnds <- list()
  common_dnds <- list()
  total_common <- list()
  total_rare <- list()
  for (y in n_npi){
    error <- poisson.test(sum_variations_all[[y]][3], sum_variations_all[[y]][1])
    rare_error_min <- append(rare_error_min, as.vector(error$conf.int[1]))
    rare_error_max <- append(rare_error_max, as.vector(error$conf.int[2]))
    error_c <- poisson.test(sum_variations_all[[y]][4], sum_variations_all[[y]][2])
    common_error_min <- append(common_error_min, as.vector(error_c$conf.int[1]))
    common_error_max <- append(common_error_max, as.vector(error_c$conf.int[2]))
    rare_dnds <- append(rare_dnds, as.vector(sum_variations_all[[y]][3]) / as.vector(sum_variations_all[[y]][1]))
    common_dnds <- append(common_dnds, as.vector(sum_variations_all[[y]][4]) / as.vector(sum_variations_all[[y]][2]))
    total_common <- append(total_common, as.vector(sum_variations_all[[y]][8]))
    total_rare <- append(total_rare, as.vector(sum_variations_all[[y]][7]))
  }
  ## add all info from functional to df
  lista <- c('rare_dnds', 'common_dnds', 'rare_error_min', 'common_error_min',
             'rare_error_max', 'common_error_max', 'total_rare', 'total_common' )
  
  for (x in 1:length(lista)){
    l = get(lista[x])
    if (x == 1) {
      alternative <- as.data.frame(l,col.names = c('npi_01', 'npi_02', 'npi_03', 'npi_04', 'npi_05','npi_06', 'npi_07', 'npi_08', 'npi_09', 'npi_1'))
    } else {
      alt_1 <- as.data.frame(l,col.names = c('npi_01', 'npi_02', 'npi_03', 'npi_04', 'npi_05','npi_06', 'npi_07', 'npi_08', 'npi_09', 'npi_1'))
      alternative <- rbind(alternative, alt_1)
    }
  } 
  rm(alt_1, rare_dnds, common_dnds, rare_error_max, rare_error_min, common_error_max, common_error_min, total_rare, total_common)
  
  # Create a dataframe with types and frequencies as values
  type_isoform <- c(rep("0-0.1", 2) , rep("0.1-0.2", 2) , rep("0.2-0.3", 2), rep("0.3-0.4", 2), rep("0.4-0.5", 2),rep("0.5-0.6", 2) , rep("0.6-0.7", 2) , rep("0.7-0.8", 2), rep("0.8-0.9", 2), rep("0.9-1", 2))
  frequency <- rep(c("Rare" , "Common" ) , 10)
  # add values to dataframe
  value <- as.matrix(c(alternative$npi_01[1], alternative$npi_01[2], alternative$npi_02[1], alternative$npi_02[2],
                       alternative$npi_03[1], alternative$npi_03[2], alternative$npi_04[1], alternative$npi_04[2],
                       alternative$npi_05[1], alternative$npi_05[2], alternative$npi_06[1], alternative$npi_06[2],
                       alternative$npi_07[1], alternative$npi_07[2], alternative$npi_08[1], alternative$npi_08[2],
                       alternative$npi_09[1], alternative$npi_09[2], alternative$npi_1[1], alternative$npi_1[2])) 
  error_min <- as.matrix(c(alternative$npi_01[3], alternative$npi_01[4], alternative$npi_02[3], alternative$npi_02[4],
                           alternative$npi_03[3], alternative$npi_03[4], alternative$npi_04[3], alternative$npi_04[4],
                           alternative$npi_05[3], alternative$npi_05[4], alternative$npi_06[3], alternative$npi_06[4],
                           alternative$npi_07[3], alternative$npi_07[4], alternative$npi_08[3], alternative$npi_08[4],
                           alternative$npi_09[3], alternative$npi_09[4], alternative$npi_1[3], alternative$npi_1[4]))
  error_max <- as.matrix(c(alternative$npi_01[5], alternative$npi_01[6], alternative$npi_02[5], alternative$npi_02[6],
                           alternative$npi_03[5], alternative$npi_03[6], alternative$npi_04[5], alternative$npi_04[6],
                           alternative$npi_05[5], alternative$npi_05[6], alternative$npi_06[5], alternative$npi_06[6],
                           alternative$npi_07[5], alternative$npi_07[6], alternative$npi_08[5], alternative$npi_08[6],
                           alternative$npi_09[5], alternative$npi_09[6], alternative$npi_1[5], alternative$npi_1[6]))
  total <- as.matrix(c(alternative$npi_01[7], alternative$npi_01[8], alternative$npi_02[7], alternative$npi_02[8],
                       alternative$npi_03[7], alternative$npi_03[8], alternative$npi_04[7], alternative$npi_04[8],
                       alternative$npi_05[7], alternative$npi_05[8], alternative$npi_06[7], alternative$npi_06[8],
                       alternative$npi_07[7], alternative$npi_07[8], alternative$npi_08[7], alternative$npi_08[8],
                       alternative$npi_09[7], alternative$npi_09[8], alternative$npi_1[7], alternative$npi_1[8]))
  alternative_dnds <- data.frame(type_isoform, frequency, value, error_min, error_max, total)
  
  ## create a column with types as factor specifying the order to plot
  alternative_dnds$type_isoform2 <- factor(alternative_dnds$type_isoform, levels = c("0-0.1", "0.1-0.2", "0.2-0.3","0.3-0.4", "0.4-0.5","0.5-0.6", "0.6-0.7", "0.7-0.8","0.8-0.9", "0.9-1" ))
  alternative_dnds$frequency2 <- factor(alternative_dnds$frequency, levels = c("Rare", "Common"))
  # barplot. Fill by frequency, y is value and x is the column of types as facfors

   if (principal == "NPI"){
    plot_funct_dnds_bef <- ggplot(alternative_dnds, aes(fill=frequency2, y=value, x=type_isoform2)) + 
      geom_bar(alpha = 0.7, position = "dodge", stat = "identity") + 
      geom_errorbar(aes(ymin=error_min, ymax=error_max), width=.2,
                    #position=position_dodge(.9), alpha = 0.5) +
                    position=position_dodge(.9)) +
      scale_fill_manual(values = c("gold1", "darkorchid3")) + 
      geom_text(
        aes(label = total, total = total + 0.05),
        position = position_dodge(0.5),
        vjust = -0.5
      )+
      theme_minimal() + 
      theme(
        text = element_text(size = 20),
        #text = element_text(family = "Arial", colour = "grey40", size = 20),
        plot.title = element_text(family = '', face = NULL, size = 22, hjust = 0.5),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(angle = 0, size = 22),
        #axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey40"),
        #axis.title.y=element_blank(),
        legend.title=element_blank()) +
      scale_y_continuous(limits=c(0,5), expand=c(0,0)) +
      #scale_x_continuous(limits=c(0,0), expand=c(0,0)) +
      ggtitle(label = "Alternative Isoforms")+
      #subtitle = "Alternative isoforms") +
      ylab("Ratio NS/Syn")+
      xlab("Normalized TRIFID score (binned)")
    
    plot_funct_dnds <- plot_funct_dnds_bef + coord_cartesian(ylim = c(0, 8))
    plot_funct_dnds
    ggsave("dnds_norm_01_npi.pdf", plot = plot_funct_dnds, device = "pdf",
           path = "/paper_set/",
           width = 15, height = 10)
    

    
    
  } else {
    plot_funct_dnds_bef <- ggplot(alternative_dnds, aes(fill=frequency2, y=value, x=type_isoform2)) + 
      geom_bar(alpha = 0.7, position = "dodge", stat = "identity") + 
      geom_errorbar(aes(ymin=error_min, ymax=error_max), width=.2,
                    position=position_dodge(.9), alpha = .5) +
      scale_fill_manual(values = c("gold1", "darkorchid3")) + 
      geom_text(
        aes(label = total, total = total + 0.05),
        position = position_dodge(0.5),
        vjust = -0.50
      )+
      theme_minimal() + 
      theme(
        #text = element_text(family = "Arial", colour = "grey40", size = 20),
        text = element_text(size = 20),
        #plot.title = element_text(size = 20, hjust = -0.2),
        plot.title = element_text(family = '', face = NULL, size = 20, hjust = 0.5),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5),
        axis.text.y = element_text(size = 20),
        #axis.title.x=element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 0, size = 20),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        #axis.title.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line = element_line(colour = "grey40"),
        legend.title=element_blank()) +
      scale_y_continuous(limits=c(0,5), expand=c(0,0)) +
      ggtitle(label = "Principal Isoforms") + 
      ylab("Ratio NS/Syn") +
      xlab("Normalized TRIFID score (binned)")
    plot_funct_dnds_bef
    plot_funct_dnds <- plot_funct_dnds_bef + coord_cartesian(ylim = c(0, 5)) #+ 

  }
  plot_funct_dnds_bef

}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(PI_plot_funct_dnds)

all <- grid.arrange(arrangeGrob(NPI_plot_funct_dnds, PI_plot_funct_dnds + theme(legend.position = "none"), 
                                nrow = 1),
                    mylegend, nrow=2,  heights=c(10, 1), top= 'Absolute score')


