library("ggplot2")
library("ggthemes")
library("gridExtra")
library("extrafont")
font_import()

exon_type = 'longest'
path <-'/home/lmartinezg/Documents/Laura/appris_clinvar/results/'
path = paste0(path,exon_type, '/')
setwd(path)
files <- list.files(pattern = "_results.gff") # get all final files from directory
sum_variations_all <- list()
for (f in 1:length(files)){
  print(files[f])
  variations <- read.table(files[f], row.names = 3, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  variations_sum <- variations[,4:11]
  sum_variations <- apply(variations_sum,2,sum)
  sum_variations_all[[f]] = sum_variations # add vector to list
}
names(sum_variations_all) = sapply(strsplit(sapply(strsplit(files, "_wvep"), '[', 1),'v37.'), '[', 2)
rare_error_min <- list()
rare_error_max <- list()
common_error_min <- list()
common_error_max <- list()
rare_dnds <- list()
common_dnds <- list()
total_common <- list()
total_rare <- list()
for (y in names(sum_variations_all)){
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

lista <- c('rare_dnds', 'common_dnds', 'rare_error_min', 'common_error_min',
           'rare_error_max', 'common_error_max', 'total_rare', 'total_common' )

for (x in 1:length(lista)){
  l = get(lista[x])
  if (x == 1) {
    #appris
    #alternative <- as.data.frame(l,col.names = c("alt_appris_longest", "alt_appris_unique",
    #                                             "alt_no_overlap", "alt_overlap", 
    #                                             "prin_appris_longest", "prin_appris_unique", 
    #                                             "prin_no_overlap", "prin_overlap"))
    
    #mane
    # alternative <- as.data.frame(l, col.names = c("alt_mane_appris","alt_mane_no_overlap",
    #                                               "alt_mane_overlap", "alt_mane_unique",
    #                                               "mane_no_overlap", "mane_overlap",
    #                                               "prin_mane_appris", "prin_mane_unique"))
   
     # longest
    alternative <- as.data.frame(l, col.names = c("alt_longest_appris","alt_longest_no_overlap",
                                                  "alt_longest_overlap", "alt_longest_unique",
                                                  "longest_no_overlap", "longest_overlap",
                                                  "prin_longest_appris", "prin_longest_unique"))
  } else {
    #appris
    #alt_1 <- as.data.frame(l,col.names = c("alt_appris_longest", "alt_appris_unique",
    #                                       "alt_no_overlap", "alt_overlap", 
    #                                       "prin_appris_longest", "prin_appris_unique", 
    #                                       "prin_no_overlap", "prin_overlap"))
    
    #mane
    # alt_1 <- as.data.frame(l, col.names = c("alt_mane_appris","alt_mane_no_overlap",
    #                                               "alt_mane_overlap", "alt_mane_unique",
    #                                               "mane_no_overlap", "mane_overlap",
    #                                               "prin_mane_appris", "prin_mane_unique"))
    
    # longest
    alt_1 <- as.data.frame(l, col.names = c("alt_longest_appris","alt_longest_no_overlap",
                                                  "alt_longest_overlap", "alt_longest_unique",
                                                  "longest_no_overlap", "longest_overlap",
                                                  "prin_longest_appris", "prin_longest_unique"))
    
    alternative <- rbind(alternative, alt_1)
  }
} 
rm(alt_1, rare_dnds, common_dnds, rare_error_max, rare_error_min, common_error_max, common_error_min, total_rare, total_common)

# Create a dataframe with types and frequencies as values

#might have to change the next line ._.

#appris
#type_isoform <- c(rep("alt_appris_longest",2), rep("alt_appris_unique",2),
#                      rep("alt_no_overlap",2), rep("alt_overlap",2), 
#                      rep("prin_appris_longest",2), rep("prin_appris_unique", 2),
#                      rep("prin_no_overlap",2), rep("prin_overlap",2))

# mane
# type_isoform <- c(rep("alt_mane_appris", 2), rep("alt_mane_no_overlap",2),
#                       rep("alt_mane_overlap",2), rep("alt_mane_unique", 2),
#                       rep("mane_no_overlap",2), rep("mane_overlap",2),
#                       rep("prin_mane_appris",2),  rep("prin_mane_unique",2))

#longest
type_isoform <- c(rep("alt_longest_appris", 2), rep("alt_longest_no_overlap",2),
                  rep("alt_longest_overlap",2), rep("alt_longest_unique", 2),
                  rep("longest_no_overlap",2), rep("longest_overlap",2),
                  rep("prin_longest_appris",2),  rep("prin_longest_unique",2))


frequency <- rep(c("Rare" , "Common" ) , 8)
# add values to dataframe

#appris
# value <- as.matrix(c(alternative$alt_appris_longest[1], alternative$alt_appris_longest[2], alternative$alt_appris_unique[1], alternative$alt_appris_unique[2],
#                      alternative$alt_no_overlap[1], alternative$alt_no_overlap[2], alternative$alt_overlap[1], alternative$alt_overlap[2],
#                      alternative$prin_appris_longest[1], alternative$prin_appris_longest[2], alternative$prin_appris_unique[1], alternative$prin_appris_unique[2],
#                      alternative$prin_no_overlap[1], alternative$prin_no_overlap[2], alternative$prin_overlap[1], alternative$prin_overlap[2]))
# 
# error_min <- as.matrix(c(alternative$alt_appris_longest[3], alternative$alt_appris_longest[4], alternative$alt_appris_unique[3], alternative$alt_appris_unique[4],
#                      alternative$alt_no_overlap[3], alternative$alt_no_overlap[4], alternative$alt_overlap[3], alternative$alt_overlap[4],
#                      alternative$prin_appris_longest[3], alternative$prin_appris_longest[4], alternative$prin_appris_unique[3], alternative$prin_appris_unique[4],
#                      alternative$prin_no_overlap[3], alternative$prin_no_overlap[4], alternative$prin_overlap[3], alternative$prin_overlap[4]))
# 
# error_max <- as.matrix(c(alternative$alt_appris_longest[5], alternative$alt_appris_longest[6], alternative$alt_appris_unique[5], alternative$alt_appris_unique[6],
#                          alternative$alt_no_overlap[5], alternative$alt_no_overlap[6], alternative$alt_overlap[5], alternative$alt_overlap[6],
#                          alternative$prin_appris_longest[5], alternative$prin_appris_longest[6], alternative$prin_appris_unique[5], alternative$prin_appris_unique[6],
#                          alternative$prin_no_overlap[5], alternative$prin_no_overlap[6], alternative$prin_overlap[5], alternative$prin_overlap[6]))
# 
# 
# total <- as.matrix(c(alternative$alt_appris_longest[7], alternative$alt_appris_longest[8], alternative$alt_appris_unique[7], alternative$alt_appris_unique[8],
#                      alternative$alt_no_overlap[7], alternative$alt_no_overlap[8], alternative$alt_overlap[7], alternative$alt_overlap[8],
#                      alternative$prin_appris_longest[7], alternative$prin_appris_longest[8], alternative$prin_appris_unique[7], alternative$prin_appris_unique[8],
#                      alternative$prin_no_overlap[7], alternative$prin_no_overlap[8], alternative$prin_overlap[7], alternative$prin_overlap[8]))

#mane
# value <- as.matrix(c(alternative$alt_mane_appris[1], alternative$alt_mane_appris[2], alternative$alt_mane_no_overlap[1], alternative$alt_mane_no_overlap[2],
#                      alternative$alt_mane_overlap[1], alternative$alt_mane_overlap[2], alternative$alt_mane_unique[1], alternative$alt_mane_unique[2],
#                      alternative$mane_no_overlap[1], alternative$mane_no_overlap[2], alternative$mane_overlap[1], alternative$mane_overlap[2],
#                      alternative$prin_mane_appris[1], alternative$prin_mane_appris[2], alternative$prin_mane_unique[1], alternative$prin_mane_unique[2]))
# 
# error_min <- as.matrix(c(alternative$alt_mane_appris[3], alternative$alt_mane_appris[4], alternative$alt_mane_no_overlap[3], alternative$alt_mane_no_overlap[4],
#                          alternative$alt_mane_overlap[3], alternative$alt_mane_overlap[4], alternative$alt_mane_unique[3], alternative$alt_mane_unique[4],
#                          alternative$mane_no_overlap[3], alternative$mane_no_overlap[4], alternative$mane_overlap[3], alternative$mane_overlap[4],
#                          alternative$prin_mane_appris[3], alternative$prin_mane_appris[4], alternative$prin_mane_unique[3], alternative$prin_mane_unique[4]))
# 
# error_max <- as.matrix(c(alternative$alt_mane_appris[5], alternative$alt_mane_appris[6], alternative$alt_mane_no_overlap[5], alternative$alt_mane_no_overlap[6],
#                          alternative$alt_mane_overlap[5], alternative$alt_mane_overlap[6], alternative$alt_mane_unique[5], alternative$alt_mane_unique[6],
#                          alternative$mane_no_overlap[5], alternative$mane_no_overlap[6], alternative$mane_overlap[5], alternative$mane_overlap[6],
#                          alternative$prin_mane_appris[5], alternative$prin_mane_appris[6], alternative$prin_mane_unique[5], alternative$prin_mane_unique[6]))
# 
# 
# total <- as.matrix(c(alternative$alt_mane_appris[7], alternative$alt_mane_appris[8], alternative$alt_mane_no_overlap[7], alternative$alt_mane_no_overlap[8],
#                      alternative$alt_mane_overlap[7], alternative$alt_mane_overlap[8], alternative$alt_mane_unique[7], alternative$alt_mane_unique[8],
#                      alternative$mane_no_overlap[7], alternative$mane_no_overlap[8], alternative$mane_overlap[7], alternative$mane_overlap[8],
#                      alternative$prin_mane_appris[7], alternative$prin_mane_appris[8], alternative$prin_mane_unique[7], alternative$prin_mane_unique[8]))

#longest
value <- as.matrix(c(alternative$alt_longest_appris[1], alternative$alt_longest_appris[2], alternative$alt_longest_no_overlap[1], alternative$alt_longest_no_overlap[2],
                     alternative$alt_longest_overlap[1], alternative$alt_longest_overlap[2], alternative$alt_longest_unique[1], alternative$alt_longest_unique[2],
                     alternative$longest_no_overlap[1], alternative$longest_no_overlap[2], alternative$longest_overlap[1], alternative$longest_overlap[2],
                     alternative$prin_longest_appris[1], alternative$prin_longest_appris[2], alternative$prin_longest_unique[1], alternative$prin_longest_unique[2]))

error_min <- as.matrix(c(alternative$alt_longest_appris[3], alternative$alt_longest_appris[4], alternative$alt_longest_no_overlap[3], alternative$alt_longest_no_overlap[4],
                         alternative$alt_longest_overlap[3], alternative$alt_longest_overlap[4], alternative$alt_longest_unique[3], alternative$alt_longest_unique[4],
                         alternative$longest_no_overlap[3], alternative$longest_no_overlap[4], alternative$longest_overlap[3], alternative$longest_overlap[4],
                         alternative$prin_longest_appris[3], alternative$prin_longest_appris[4], alternative$prin_longest_unique[3], alternative$prin_longest_unique[4]))

error_max <- as.matrix(c(alternative$alt_longest_appris[5], alternative$alt_longest_appris[6], alternative$alt_longest_no_overlap[5], alternative$alt_longest_no_overlap[6],
                         alternative$alt_longest_overlap[5], alternative$alt_longest_overlap[6], alternative$alt_longest_unique[5], alternative$alt_longest_unique[6],
                         alternative$longest_no_overlap[5], alternative$longest_no_overlap[6], alternative$longest_overlap[5], alternative$longest_overlap[6],
                         alternative$prin_longest_appris[5], alternative$prin_longest_appris[6], alternative$prin_longest_unique[5], alternative$prin_longest_unique[6]))


total <- as.matrix(c(alternative$alt_longest_appris[7], alternative$alt_longest_appris[8], alternative$alt_longest_no_overlap[7], alternative$alt_longest_no_overlap[8],
                     alternative$alt_longest_overlap[7], alternative$alt_longest_overlap[8], alternative$alt_longest_unique[7], alternative$alt_longest_unique[8],
                     alternative$longest_no_overlap[7], alternative$longest_no_overlap[8], alternative$longest_overlap[7], alternative$longest_overlap[8],
                     alternative$prin_longest_appris[7], alternative$prin_longest_appris[8], alternative$prin_longest_unique[7], alternative$prin_longest_unique[8]))



alternative_dnds <- data.frame(type_isoform, frequency, value, error_min, error_max, total)

## create a column with types as factor specifying the order to plot

#appris
#alternative_dnds$type_isoform2 <- factor(alternative_dnds$type_isoform, levels = c("alt_appris_longest", "alt_appris_unique",
#                                                                                   "alt_no_overlap", "alt_overlap", 
#                                                                                   "prin_appris_longest", "prin_appris_unique", 
#                                                                                   "prin_no_overlap", "prin_overlap"))

#mane
# alternative_dnds$type_isoform2 <- factor(alternative_dnds$type_isoform, levels = c("alt_mane_appris","alt_mane_no_overlap",
#                                                                                     "alt_mane_overlap", "alt_mane_unique",
#                                                                                     "mane_no_overlap", "mane_overlap",
#                                                                                     "prin_mane_appris", "prin_mane_unique"))
#longest
alternative_dnds$type_isoform2 <- factor(alternative_dnds$type_isoform, levels = c("alt_longest_appris","alt_longest_no_overlap",
                                                                                   "alt_longest_overlap", "alt_longest_unique",
                                                                                   "longest_no_overlap", "longest_overlap",
                                                                                   "prin_longest_appris", "prin_longest_unique"))


alternative_dnds$frequency2 <- factor(alternative_dnds$frequency, levels = c("Rare", "Common"))

# Barplot. Fill by frequency, y is value and x is the column of types as factors

plot_funct_dnds_bef <- ggplot(alternative_dnds, aes(fill=frequency2, y=value, x=type_isoform2)) + 
  geom_bar(alpha = 0.7, position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=error_min, ymax=error_max), width=.2,
                position=position_dodge(.9)) +
  geom_text(
    aes(label = total, total = total + 0.05),
    position = position_dodge(0.5),
    vjust = -0.5, family = 'ArialMT')+ 
  theme_minimal() +
  theme_set(theme_gray(base_size = 20, base_family = 'Arial' )) +
  theme(
    text = element_text(family = 'ArialMT', size = 20),
    plot.title = element_text(family = 'ArialMT', face = NULL, size = 22, hjust = 0.5),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 25, size = 22, vjust = 0.99, hjust = 0.85),
    axis.text.y = element_text(angle = 0, size = 22),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "grey40"),
    legend.title=element_blank()) +
  ggtitle(label = "LONGEST set") +
  ylab("Ratio NS/Syn")

plot_funct_dnds_bef


plot_funct_dnds <- plot_funct_dnds_bef + coord_cartesian(ylim = c(0, 4))
plot_funct_dnds
ggsave("nsynsyn_longest.pdf", plot = plot_funct_dnds, device = "pdf",
       path = "/home/lmartinezg/Documents/Laura/appris_clinvar/results/longest",
       width = 15, height = 10)
