#Load Libraries
library(tidyverse)
library(magrittr)
library(ggstatsplot)
library(viridis)
library(ggpubr)
library(ggpmisc)
library(ggstatsplot)
library(ggthemes)

#Load Data
my.data = readRDS("./data/processed_data/my.data.RDS")
abundance.estimates = readRDS("./data/processed_data/abundance.estimates.RDS")
primer.set.error = read_csv("./data/processed_data/primer_set_error.csv")

#Seperate into samples, positive controls, and negative controls
mc.data = my.data %>% filter(community == "1"|community == "2"|community == "3"|community == "4"|community == "5")
pos.control.data = my.data %>% filter(community == "twist"|community == "heat_sars_pbs"|community == "heat_sars_wastewater")
neg.control.data = my.data %>% filter(community == "pcr_water"|community == "wastewater")

mc.abundance.estimates = abundance.estimates %>% filter(community == "1"|community == "2"|community == "3"|community == "4"|community == "5")
#############################################################################
######################Read Volume############################################

#Number of reads before/after filtering, overall. 
my.data %>% 
  group_by(read_length) %>%
  summarize(total_raw_reads = sum(raw_reads), 
            total_filtered_reads = sum(filtered_reads)) %>%
  mutate(total_filtered_reads/total_raw_reads*100)

#Summarize mc.data
summary(mc.data)

#Create a summary table
mc.stat.table = mc.data %>% 
  group_by(primer_set, read_length) %>% 
  summarize(
    raw_reads_total = sum(raw_reads), 
    filtered_reads_total = sum(filtered_reads), 
    reads_mapped_kallisto_total = sum(reads_mapped_kallisto), 
    raw_reads.avg = mean(raw_reads), 
    raw_reads.sd = sd(raw_reads), 
    filtered_reads.avg = mean(filtered_reads), 
    filtered_reads.sd = sd(filtered_reads), 
    reads_mapped_kallisto.avg = mean(reads_mapped_kallisto), 
    reads_mapped_kallisto.avg = sd(reads_mapped_kallisto), 
    filtered_raw_reads.avg = mean(filtered_raw_reads), 
    filtered_raw_reads.sd = sd(filtered_raw_reads), 
    assigned_filtered_reads.avg = mean(assigned_filtered_reads),
    assigned_filtered_reads.sd = sd(assigned_filtered_reads), 
    assigned_raw_reads.avg = mean(assigned_raw_reads), 
    assigned_raw_reads.sd = sd(assigned_raw_reads),
    genome_coverage.avg = mean(percent_genome), 
    genome_coverage.sd= sd(percent_genome),
    depth.avg = mean(depth_avg), 
    depth.sd = sd(depth_avg), 
    percent_genome_20X.avg = mean(percent_genome_20X), 
    percent_genome_20X.sd = sd(percent_genome_20X))
    
write.csv(mc.stat.table, "./data/processed_data/mc.stat.table.csv") 

############################################################################
#######################Read Retention#######################################

mc.data %>% 
  grouped_ggbetweenstats(
    x = read_length, 
    y = filtered_raw_reads, 
    grouping.var = primer_set, 
    type = "np"
  )


mc.data %>% 
  ggbetweenstats(
    x = primer_set, 
    y = filtered_raw_reads, 
    type = "np"
  )

#############################################################################
###########################Genome Coverage###################################


#How many samples have >20X coverage on average? 
mc.data %>% 
  filter(depth_avg > 19) %>% 
  tally()

#How many samples have an average sequencing depth of >100X?
mc.data %>% 
  filter(depth_avg > 99) %>% 
  tally()



#How does coverage observed compare to coverage expected?

#expected coverage
expected_coverage = seq(10000, 100000000, by= 500)
expected_coverage = as.data.frame(expected_coverage)
expected_coverage$raw_bases = expected_coverage$expected_coverage

expected_coverage = expected_coverage %>%
  mutate(expected_coverage = raw_bases/29901*100) %>%
  mutate(expected_coverage = case_when(expected_coverage > 100 ~ 100, 
                                       TRUE ~ expected_coverage)) %>% 
  mutate(expected_coverage_20X = raw_bases/598020*100) %>%
  mutate(expected_coverage_20X = case_when(expected_coverage_20X > 100 ~ 100, 
                                           TRUE ~ expected_coverage_20X)) 






#How does sequencing effort affect genomic coverage? 
bases_coverage = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot() + 
  geom_point(aes(x = log10(raw_bases), y = percent_genome, color = seq_approach)) + 
  geom_line(data = expected_coverage, aes(x = log10(raw_bases), y = expected_coverage), linetype = "dotted") +
  xlab("Log10 (Bases)") + 
  ylab("Genomic Coverage (%)") + 
  labs(color = "Sequencing Approach") + 
  stat_cor(aes(x = log10(raw_bases), y = percent_genome), method = "spearman")   +
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_colour_brewer(palette = "Paired") 



#How does sequencing effort affect depth? 
bases_depth = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  mutate(expected_depth = raw_bases/29903) %>%
  mutate(observed_expected_depth = depth_avg/expected_depth) %>%
  ggplot() + 
  geom_line(aes(x = log10(raw_bases), y = log10(expected_depth)), linetype = "dotted") +
  geom_point(aes(x = log10(raw_bases), y = log10(depth_avg), color = seq_approach)) + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_colour_brewer(palette = "Paired") + 
  xlab("Log10 (Bases)") + 
  ylab("Log10 (Depth)") + 
  labs(color = "Sequencing Approach") + 
  stat_cor(aes(x = log10(raw_bases), y = log10(depth_avg)), method = "spearman")   


#How does sequencing effort affect depth? 
bases_20x = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot() + 
  geom_point(aes(x = log10(raw_bases), y = percent_genome_20X, color = seq_approach)) + 
  geom_line(data = expected_coverage, aes(x = log10(raw_bases), y = expected_coverage_20X), linetype = "dotted") +
  xlab("Log10 (Bases)") + 
  ylab("Genomic Coverage \n >20X (%)") + 
  labs(color = "Sequencing Approach") + 
  stat_cor(aes(x = log10(raw_bases), y = percent_genome), method = "spearman")   +
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_colour_brewer(palette = "Paired") +
  ylim(0,100) 


seq_stats = ggarrange(bases_coverage, 
                      bases_depth, 
                      bases_20x, 
                      ncol = 1,
                      common.legend = TRUE, 
                      legend = "bottom", 
                      labels =  c("A", "B", "C"), 
                      align = "hv")

#tiff('./figures/seq_stats_alt.tiff', units="in", width = 7, height = 9, res=600, compression = 'lzw', pointsize = 12)
plot(seq_stats)
#dev.off()

#########################################################################################################################
##########################################################Depth Expected#################################################

mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  mutate(expected_depth = raw_bases/30000) %>%
  mutate(observed_expected_depth = depth_avg/expected_depth) %>%
  ggplot(aes(x = primer_set, y = observed_expected_depth)) + 
  geom_boxplot() + 
  stat_compare_means()


  
#######################################################

mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  #filter(raw_bases > 10000000) %>% 
  ggplot(aes(x = log10(raw_bases), y = percent_genome, color = primer_set)) + 
  geom_point() + 
  geom_smooth()


mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  filter(raw_bases > 10000000) %>% 
  ggplot(aes(x = log10(raw_bases), y = percent_genome_20X, color = primer_set)) + 
  geom_point() + 
  geom_smooth()

mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  filter(raw_bases > 10000000) %>% 
  group_by(primer_set) %>% 
  summarize(mean = mean(percent_genome), 
            sd = sd(percent_genome))

mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  filter(raw_bases > 10000000) %>% 
  group_by(primer_set) %>% 
  summarize(mean = mean(percent_genome_20X), 
            sd = sd(percent_genome_20X))

########################################################################################
#########################RNA Concentration vs. Coverage#################################

dilution_sequences = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  ggplot(aes(x = dilution, y = log10(raw_bases))) + 
  geom_boxplot(aes(fill = dilution)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Concentration (cp/uL)") + 
  ylab("Log10 (Bases)") + 
  xlab("") 

dilution_depth = mc.data %>% 
  ggplot(aes(x = dilution, y = log10(depth_avg))) + 
  geom_boxplot(aes(fill = dilution)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Concentration (cp/uL)") + 
  ylab("Log10 (Depth)") + 
  xlab("")

dilution_coverage = mc.data %>% 
  ggplot(aes(x = dilution, y = percent_genome)) + 
  geom_boxplot(aes(fill = dilution)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Concentration (cp/uL)") + 
  ylab("Genome Coverage (%)") + 
  xlab("") + 
  ylim(0,100)

dilution_20X = mc.data %>% 
  ggplot(aes(x = dilution, y = percent_genome_20X)) + 
  geom_boxplot(aes(fill = dilution)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Concentration (cp/uL)") + 
  ylab("Genome Coverage \n >20X (%)") + 
  xlab("") + 
  ylim(0,100)

dilution_stats = ggarrange(dilution_sequences + rremove("x.text"), 
                      dilution_depth + rremove("x.text"), 
                      dilution_coverage + rremove("x.text"), 
                      dilution_20X + rremove("x.text"), 
                      common.legend = TRUE, 
                      legend = "bottom", 
                      ncol = 1, 
                      labels =  c("A", "B", "C", "D"), 
                      align = "hv")

#tiff('./figures/dilution_stats.tiff', units="in", width = 7, height = 9, res=600, compression = 'lzw', pointsize = 12)
plot(dilution_stats)
#dev.off()


########################################################################################
########################Mock Community vs. Coverage#####################################

community_sequences = mc.data %>% 
  mutate(raw_bases = case_when(read_length == "PE150" ~ filtered_reads*300, 
                               TRUE ~ filtered_reads*500)) %>%
  ggplot(aes(x = community, y = log10(raw_bases))) + 
  geom_boxplot(aes(fill = community)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Mock Community No.") + 
  ylab("Log10 (Bases)") + 
  xlab("") 

community_depth = mc.data %>% 
  ggplot(aes(x = community, y = log10(depth_avg))) + 
  geom_boxplot(aes(fill = community)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Mock Community No.") + 
  ylab("Log10 (Depth)") + 
  xlab("")

community_coverage = mc.data %>% 
  ggplot(aes(x = community, y = percent_genome)) + 
  geom_boxplot(aes(fill = community)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Mock Community No.") + 
  ylab("Genome Coverage (%)") + 
  xlab("") + 
  ylim(0,100)

community_20X = mc.data %>% 
  ggplot(aes(x = community, y = percent_genome_20X)) + 
  geom_boxplot(aes(fill = community)) + 
  facet_wrap(~primer_set, ncol = 3) +
  stat_compare_means(label.y.npc = "bottom") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Mock Community No.") + 
  ylab("Genome Coverage \n >20X (%)") + 
  xlab("") + 
  ylim(0,100)

community_stats = ggarrange(community_sequences + rremove("x.text"), 
                           community_depth + rremove("x.text"), 
                           community_coverage + rremove("x.text"), 
                           community_20X + rremove("x.text"), 
                           common.legend = TRUE, 
                           legend = "bottom", 
                           ncol = 1, 
                           labels =  c("A", "B", "C", "D"), 
                           align = "hv")

#tiff('./figures/community_stats.tiff', units="in", width = 10, height = 9, res=600, compression = 'lzw', pointsize = 12)
plot(community_stats)
#dev.off()

########################################################################################
############################Assignments in Kallisto#####################################

shapiro.test(mc.data$assigned_filtered_reads)

mc.data %>% 
 ggstatsplot::ggbetweenstats(
   x = primer_set, 
   y = assigned_filtered_reads, 
   type = "np") 


#########################################################################################
#################################Abundance Estimates#####################################

#Abundance estimates overall 
abundance_primerset = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes( fill = lineage_name), size = 2, shape = 21) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = primer.set.error, aes(x = 25, y = 95, label = paste0("Rsq = ", round(R_2, 2)))) +
  facet_wrap(~primer_set) + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(fill = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_few() 



#Abundance estimates overall 
error_primerset = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  ggplot(aes(x = lineage_name, y = abundance_observed-abundance_expected)) +
  geom_boxplot(aes( fill = lineage_name)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~primer_set) + 
  ylab("Error (%)") + 
  xlab("") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(fill = "Variant") +
  theme_few()

abundance_estimates = ggarrange(abundance_primerset, 
          error_primerset,
          ncol = 1, 
          common.legend = TRUE, 
          align = "hv", 
          labels =  c("A", "B") 
          
)

#tiff('./figures/abundance_estimates.tiff', units="in", width = 8, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(abundance_estimates)
#dev.off()


####################################################################################
##########################Wuhan vs. Variants########################################


mc.abundance.estimates %>% 
  mutate(detection = case_when(abundance_observed > 0 ~ "yes", 
                               TRUE ~ "no")) %>% 
  filter(primer_set == "Midnight") %>%
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  ggstatsplot::grouped_ggbarstats(
    y = abundance_expected, 
    x = detection, 
    grouping.var = lineage_name) 
  

#################################################################################################
##########################################Abundances Stacked#####################################

observed = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  ggplot(aes(x = seq_approach, y = abundance_observed)) + 
  geom_col(aes(fill = lineage_name)) + 
  facet_wrap(~community + dilution, ncol = 3) + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Sequencing Approach") + 
  ylab("Abundance (%)") + 
  labs(fill = "Variant")




expected = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  mutate(seq_approach = paste0(primer_set, "-", read_length)) %>%
  filter(primer_set == "V4" & read_length == "PE250" & dilution == "10^3") %>%
  ggplot(aes(x = "Expected", y = abundance_expected)) + 
  geom_col(aes(fill = lineage_name)) + 
  facet_wrap(~community, ncol = 1) + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  xlab("") + 
  ylab("Abundance (%)") + 
  labs(fill = "Variant")


abundance_stacked = ggarrange(expected, 
                              observed + rremove("y.text") + rremove("ylab") + rremove("xlab"),
                              ncol = 2, 
                              common.legend = TRUE, 
                              legend = "bottom", 
                              widths = c(0.4, 3),
                              align = "h")


#tiff('./figures/abundance_stacked.tiff', units="in", width = 9, height = 9, res=600, compression = 'lzw', pointsize = 12)
plot(abundance_stacked)
#dev.off()

#################################################################################################
##########################################Variants: Other########################################


#Which variants are "other"? Which is the most common?

mc.abundance.estimates %>% 
  filter(lineage_name != "Wuhan" & lineage_name != "Alpha" & lineage_name != "Beta" & lineage_name != "Delta") %>%
  ggplot(aes(x = lineage_name, y = abundance_observed)) +
  geom_boxplot() + 
  xlab("Variant") + 
  ylab("Abundance Observed (%)") + 
  theme_ggstatsplot() + stat_compare_means()


mc.abundance.estimates %>% 
  ungroup() %>%
  filter(lineage_name != "Wuhan" & lineage_name != "Alpha" & lineage_name != "Beta" & lineage_name != "Delta") %>%
  filter(abundance_observed > 0) %>% 
  group_by(primer_set) %>% 
  tally()


mc.abundance.estimates %>% 
  filter(lineage_name != "Wuhan" & lineage_name != "Alpha" & lineage_name != "Beta" & lineage_name != "Delta") %>%
  ggplot(aes(x = community, y = abundance_observed)) +
  geom_boxplot() + 
  facet_wrap(~lineage_name) +
  xlab("Variant") + 
  ylab("Abundance Observed (%)") + 
  theme_ggstatsplot() + 
  stat_compare_means()


mc.abundance.estimates %>% 
  filter(lineage_name != "Wuhan" & lineage_name != "Alpha" & lineage_name != "Beta" & lineage_name != "Delta") %>%
  ggplot(aes(x = community, y = abundance_observed)) +
  geom_boxplot() + 
  facet_wrap(~lineage_name) +
  xlab("Variant") + 
  ylab("Abundance Observed (%)") + 
  theme_ggstatsplot() + 
  stat_compare_means()







  


