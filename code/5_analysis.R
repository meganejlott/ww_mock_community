#Load Libraries
library(tidyverse)
library(purrr)
library(magrittr)
library(viridis)
library(Boruta)
library(ggstatsplot)
library(ggpubr)
library(ggthemes)

#Load data 
my.data = readRDS("./data/processed_data/my.data.RDS")

################################################################################
#############################Predictive Variables###############################

#Select mock community samples 
mc.data = my.data %>% filter(community == "1"|community == "2"|community == "3"|community == "4"|community == "5")
mc.data %<>% mutate(bases_filtered = case_when(read_length == "PE150" ~ filtered_reads*300, TRUE ~ filtered_reads*500)) 
                    
#Which are the most important variables in determining R^2?
#Outcome: R_2
df = mc.data %>% 
  select(filtered_reads, 
         bases_filtered,
         reads_mapped_kallisto, 
         percent_genome, 
         depth_avg, 
         percent_genome_20X, 
         primer_set,
         dilution, 
         read_length, 
         community, 
         R_2) %>% 
  drop_na()

# Run Boruta Algorithm
set.seed(456)
boruta <- Boruta(R_2~., data = df, doTrace = 2)
print(boruta)
plot(boruta)


plot(boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta$ImpHistory),function(i)
boruta$ImpHistory[is.finite(boruta$ImpHistory[,i]),i])
names(lz) <- colnames(boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(boruta$ImpHistory), cex.axis = 0.45)


final.boruta <- TentativeRoughFix(boruta)
print(final.boruta)

getSelectedAttributes(final.boruta, withTentative = F)

boruta.df <- attStats(final.boruta)
class(boruta.df)
print(boruta.df)


#######################################################################
###########################Predictive Variables########################



#Number of raw bases 
accuracy_bases = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                                TRUE ~ filtered_reads*500)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = log10(raw_bases), y = R_2)) + 
  geom_point(aes(color = seq_approach)) + 
  #facet_wrap(~read_length, ncol = 1) + 
  xlab("Log10 (Bases)") + 
  ylab(expression(R^2)) + 
  labs(color = "Sequencing Approach") + 
  stat_cor(method = "spearman")   +
  scale_colour_brewer(palette = "Paired") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))

  
#Genome Coverage
accuracy_coverage = mc.data %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = percent_genome, y = R_2)) + 
  geom_point(aes(color = seq_approach)) + 
  #facet_wrap(~read_length, ncol = 1) + 
  xlab("Genome Coverage (%)") + 
  ylab(expression(R^2)) + 
  labs(color = "Sequencing Approach") + 
  stat_cor(method = "spearman") +
  scale_colour_brewer(palette = "Paired") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))

  

#Sequencing Depth
accuracy_depth = mc.data %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = log10(depth_avg), y = R_2)) + 
  geom_point(aes(color = seq_approach)) + 
  #facet_wrap(~read_length, ncol = 1) + 
  xlab("Log10(Sequencing Depth)") + 
  ylab(expression(R^2)) + 
  labs(color = "Sequencing Approach") + 
  stat_cor(method = "spearman") + 
  scale_colour_brewer(palette = "Paired") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))


#Genome Coverage >20X
accuracy_coverage_20x = mc.data %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = percent_genome_20X, y = R_2)) + 
  geom_point(aes(color = seq_approach)) + 
  #facet_wrap(~read_length, ncol = 1) + 
  xlab("Genome Coverage >20X (%)") + 
  ylab(expression(R^2)) + 
  labs(color = "Sequencing Approach") + 
  stat_cor(method = "spearman") + 
  scale_colour_brewer(palette = "Paired") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))




predictive_variables = ggarrange(accuracy_bases, 
                                 accuracy_coverage, 
                                 accuracy_depth, 
                                 accuracy_coverage_20x, 
                                 common.legend = TRUE, 
                                 legend = "bottom", 
                                 labels = c("A", "B", "C", "D"))

predictive_variables = ggarrange(accuracy_bases, accuracy_coverage, accuracy_depth, accuracy_coverage_20x, 
                                 align = "hv",
                                 common.legend = TRUE, 
                                 legend = "bottom", 
                                 labels = c("A", "B", "C", "D"))

#tiff('./figures/predictive_variables.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(predictive_variables)
#dev.off()



###############################################################################################
#############Rsq as a function of Filtered Reads  & Depth####################################

reads_coverage_accuracy = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  ggplot(aes(x = log10(raw_bases), y = percent_genome)) + 
  geom_point(aes(color = R_2)) + 
  facet_wrap(~primer_set) + 
  scale_color_viridis(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + 
  xlab("Log10 (Bases)") + 
  ylab("Genome Coverage (%)") + 
  labs(color = expression(R^2)) + 
  stat_cor(method = "spearman") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))

reads_depth_accuracy = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  ggplot(aes(x = log10(raw_bases), y = log10(depth_avg))) + 
  geom_point(aes(color = R_2)) + 
  facet_wrap(~primer_set) + 
  scale_color_viridis(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + 
  xlab("Log10 (Bases)") + 
  ylab("Log10 (Depth)") + 
  labs(color = expression(R^2)) + 
  stat_cor(method = "spearman") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))


reads_20x_accuracy = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  ggplot(aes(x = log10(raw_bases), y = percent_genome_20X)) + 
  geom_point(aes(color = R_2)) + 
  facet_wrap(~primer_set) + 
  scale_color_viridis(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) + 
  xlab("Log10 (Bases)") + 
  ylab("Genome Coverage >20X (%)") + 
  labs(color = expression(R^2)) + 
  stat_cor(method = "spearman") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2"))



predictive_variables2 = ggarrange(reads_coverage_accuracy + rremove("xlab"), 
                                  reads_depth_accuracy + rremove("xlab"), 
                                  reads_20x_accuracy, 
                                  ncol = 1,
                                  align = "v",
                                  common.legend = TRUE, 
                                  heights = c(1,1,1.1),
                                  legend = "right", 
                                  labels = c("A", "B", "C"))

#tiff('./figures/predictive_variables2.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(predictive_variables2)
#dev.off()




################################################################################################
##################################Community Composition#########################################

my_comparisons <- list( c("1", "2"), 
                        c("1", "3"), 
                        c("1", "4"), 
                        c("1", "5"), 
                        c("2", "3"), 
                        c("2", "4"), 
                        c("2", "5"), 
                        c("3", "4"), 
                        c("3", "5"), 
                        c("4", "5") )

compare_means(R_2 ~ community,  data = mc.data)

composition_accuracy = mc.data %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = community, y = R_2, fill = community)) + 
  geom_boxplot() + 
  facet_wrap(~primer_set, ncol = 3) + 
  xlab("Mock Community") + 
  ylab(expression(R^2)) + 
  labs(fill = "Mock Community") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  stat_compare_means() 


mc.data %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  grouped_ggbetweenstats(
    x = community, 
    y = R_2, 
    grouping.var = primer_set, 
    type = "np") 

  
mc.data %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggbetweenstats(
    x = community, 
    y = R_2, 
    type = "np")




wuhan_accuracy = mc.data %>% 
  mutate(abundance_expected_Wuhan = as.factor(abundance_expected_Wuhan)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = abundance_expected_Wuhan, y = R_2)) + 
  geom_boxplot(aes(fill = abundance_expected_Wuhan)) + 
  stat_cor(method = "spearman") +
  facet_wrap(~primer_set, ncol = 3) + 
  xlab("Expected Abundance Wuhan (%)") + 
  ylab(expression(R^2)) + 
  labs(fill = "Expected Abundance (%)") + 
  scale_fill_brewer(palette = "Reds") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) 



alpha_accuracy = mc.data %>% 
  mutate(abundance_expected_Alpha = as.factor(abundance_expected_Alpha)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = abundance_expected_Alpha, y = R_2)) + 
  geom_boxplot(aes(fill = abundance_expected_Alpha)) + 
  stat_cor(method = "spearman") +
  facet_wrap(~primer_set, ncol = 3) + 
  xlab("Expected Abundance Alpha (%)") + 
  ylab(expression(R^2)) + 
  labs(fill = "Expected Abundance (%)") + 
  scale_fill_brewer(palette = "Blues") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) 


beta_accuracy = mc.data %>% 
  mutate(abundance_expected_Beta = as.factor(abundance_expected_Beta)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = abundance_expected_Beta, y = R_2)) + 
  geom_boxplot(aes(fill = abundance_expected_Beta)) + 
  stat_cor(method = "spearman") +
  facet_wrap(~primer_set, ncol = 3) + 
  xlab("Expected Abundance Beta (%)") + 
  ylab(expression(R^2)) + 
  labs(fill = "Expected Abundance (%)") + 
  scale_fill_brewer(palette = "Greens") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  stat_compare_means()



delta_accuracy = mc.data %>% 
  mutate(abundance_expected_Delta = as.factor(abundance_expected_Delta)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = abundance_expected_Delta, y = R_2)) + 
  geom_boxplot(aes(fill = abundance_expected_Delta)) + 
  stat_cor(method = "spearman") +
  facet_wrap(~primer_set, ncol = 3) + 
  xlab("Expected Abundance Delta (%)") + 
  ylab(expression(R^2)) + 
  labs(fill = "Expected Abundance (%)") + 
  scale_fill_brewer(palette = "Purples") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) 



accuracy_variant = ggarrange(wuhan_accuracy + rremove("xlab"), 
                             alpha_accuracy + rremove("xlab"), 
                             beta_accuracy + rremove("xlab"), 
                             delta_accuracy, 
                             ncol = 1, 
                             labels = c("Wuhan", "Alpha", "Beta", "Delta"))

plot(accuracy_variant)


#######################################################################################################
########################################Rejected Variables#############################################


#Read Length

read_length_accuracy = mc.data %>% 
  ggplot(aes(x = primer_set, y = R_2, fill = read_length)) + 
  geom_boxplot() + 
  xlab("") +
  ylab(expression(R^2)) + 
  labs(fill = "Read Length") + 
  scale_fill_brewer(palette = "Set1") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  stat_compare_means()

#Concentration

concentration_accuracy = mc.data %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = primer_set, y = R_2, fill = dilution)) + 
  geom_boxplot() + 
  #facet_wrap(~read_length, ncol = 1) + 
  xlab("Primer Set") + 
  ylab(expression(R^2)) + 
  labs(fill = "Concentration \n (cp/uL)") + 
  scale_fill_brewer(palette = "Set1") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  stat_compare_means()



