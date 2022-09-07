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

#Separate into samples, positive controls, and negative controls
mc.data = my.data %>% filter(community == "1"|community == "2"|community == "3"|community == "4"|community == "5")
pos.control.data = my.data %>% filter(community == "twist"|community == "heat_sars_pbs"|community == "heat_sars_wastewater")
neg.control.data = my.data %>% filter(community == "pcr_water"|community == "wastewater")

mc.abundance.estimates = abundance.estimates %>% filter(community == "1"|community == "2"|community == "3"|community == "4"|community == "5")
pos.control.abundance.estimates = abundance.estimates %>% filter(community == "twist"|community == "heat_sars_pbs"|community == "heat_sars_wastewater")
neg.control.abundance.estimates = abundance.estimates %>% filter(community == "pcr_water"|community == "wastewater")


#############################################################################
##########################Positive Controls##################################

twist_abundance = pos.control.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(community == "twist") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = pos.control.data %>% filter(community == "twist"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_wrap(primer_set~read_length, ncol = 2) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()

#tiff('./figures/twist_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(twist_abundance)
#dev.off()


sars_pbs_abundance = pos.control.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(community == "heat_sars_pbs") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = pos.control.data %>% filter(community == "heat_sars_pbs"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_wrap(primer_set~read_length, ncol = 2) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()

#tiff('./figures/sars_pbs_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(sars_pbs_abundance)
#dev.off()

sars_ww_abundance = pos.control.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(community == "heat_sars_wastewater") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = pos.control.data %>% filter(community == "heat_sars_wastewater"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_wrap(primer_set~read_length, ncol = 2) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()

#tiff('./figures/sars_ww_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(sars_ww_abundance)
#dev.off()

#####################################################################################################
#############################Negative Controls#######################################################

neg.control.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(community == "pcr_water") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = neg.control.data %>% filter(community == "pcr_water"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_wrap(primer_set~read_length, ncol = 2) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


neg.control.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(community == "wastewater") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = neg.control.data %>% filter(community == "wastewater"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_wrap(primer_set~read_length, ncol = 2) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()



#############################################################################
######################Mock Communities#######################################


#Abundance estimates by sample

#Midnight annotations
Midnight_abundance_estimates = 
  mc.data %>% select(sample_id, primer_set, community, dilution, read_length, R_2) %>% 
  filter(primer_set == "Midnight")


midnight_estimates_PE150 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "Midnight" & read_length == "PE150") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = Midnight_abundance_estimates %>% filter(read_length == "PE150"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()



midnight_estimates_PE250 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "Midnight" & read_length == "PE250") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = Midnight_abundance_estimates %>% filter(read_length == "PE250"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


Midnight_abundance = ggarrange(midnight_estimates_PE150, 
                               midnight_estimates_PE250, 
                               common.legend = TRUE, 
                               legend = "right",
                               ncol = 1,
                               labels = c("PE150", "PE250"), 
                               align = "hv") 

#tiff('./figures/Midnight_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(Midnight_abundance)
#dev.off()


#V4 annotations
V4_abundance_estimates = 
  mc.data %>% select(sample_id, primer_set, community, dilution, read_length, R_2) %>% 
  filter(primer_set == "V4")


V4_estimates_PE150 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "V4" & read_length == "PE150") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = V4_abundance_estimates %>% filter(read_length == "PE150"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


V4_estimates_PE250 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "V4" & read_length == "PE250") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = V4_abundance_estimates %>% filter(read_length == "PE250"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


V4_abundance = ggarrange(V4_estimates_PE150, 
                               V4_estimates_PE250, 
                               common.legend = TRUE, 
                               legend = "right",
                               ncol = 1,
                               labels = c("PE150", "PE250"), 
                               align = "hv") 

#tiff('./figures/V4_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(V4_abundance)
#dev.off()




#VarSkip annotations
VarSkip_abundance_estimates = 
  mc.data %>% select(sample_id, primer_set, community, dilution, read_length, R_2) %>% 
  filter(primer_set == "VarSkip")


VarSkip_estimates_PE150 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "VarSkip" & read_length == "PE150") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = VarSkip_abundance_estimates %>% filter(read_length == "PE150"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


VarSkip_estimates_PE250 = mc.abundance.estimates %>% 
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(primer_set == "VarSkip" & read_length == "PE250") %>%
  ggplot(aes(x = abundance_expected, y = abundance_observed)) +
  geom_point(aes(color = lineage_name)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_text(data = VarSkip_abundance_estimates %>% filter(read_length == "PE250"), 
            aes(x = 25, y = 95, label = paste0("Rsq =", round(R_2, 2)))) +
  facet_grid(dilution~community) + 
  scale_colour_brewer(palette = "Dark2") + 
  xlab("Abundance Expected (%)") + 
  ylab("Abundance Observed (%)") + 
  labs(color = "Variant", shape = "Read Length") + 
  xlim(0,100) + 
  ylim(0,100) + 
  theme_few()


VarSkip_abundance = ggarrange(VarSkip_estimates_PE150, 
                         VarSkip_estimates_PE250, 
                         common.legend = TRUE, 
                         legend = "right",
                         ncol = 1,
                         labels = c("PE150", "PE250"), 
                         align = "hv") 

#tiff('./figures/VarSkip_abundance.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(VarSkip_abundance)
#dev.off()



###########################################################################################################################
#####################################Abundance vs. Composition############################################################


accuracy_abundance = mc.abundance.estimates %>% 
  left_join(mc.data %>% select(sample_id, read_length, R_2)) %>%
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(lineage_name != "Other") %>%
  ggplot(aes(x = abundance_expected, y = R_2, color = lineage_name)) + 
  geom_point() +
  facet_grid(lineage_name~primer_set) + 
  ylab(expression(R^2)) + 
  scale_color_brewer(palette = "Dark2") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  xlab("Abundance Expected (%)") +
  stat_cor(method = "spearman", label.y.npc = "bottom") +
  theme(legend.position="none")

tiff('./figures/accuracy_abundance.tiff', units="in", width = 7, height = 7, res=600, compression = 'lzw', pointsize = 12)
plot(accuracy_abundance)
dev.off()



accuracy_abundance2 = mc.abundance.estimates %>% 
  left_join(mc.data %>% select(sample_id, read_length, R_2)) %>%
  mutate(lineage_name = case_when(lineage_name == "Wuhan" ~ "Wuhan", 
                                  lineage_name == "Alpha" ~ "Alpha", 
                                  lineage_name == "Beta" ~ "Beta", 
                                  lineage_name == "Delta" ~ "Delta", 
                                  TRUE ~ "Other")) %>%
  mutate(lineage_name = factor(lineage_name, ordered = TRUE, 
                               levels = c("Wuhan", "Alpha", "Beta", "Delta", "Other"))) %>% 
  filter(lineage_name != "Other") %>%
  mutate(abundance_expected = as.factor(abundance_expected)) %>%
  ggplot() + 
  geom_boxplot(aes(x = abundance_expected, y = R_2, color = lineage_name)) +
  geom_point(aes(x = abundance_expected, y = R_2, color = lineage_name)) +
  facet_grid(lineage_name~primer_set) + 
  ylab(expression(R^2)) + 
  scale_color_brewer(palette = "Dark2") + 
  theme_few() + 
  theme(panel.grid = element_line(colour = "snow2")) + 
  xlab("Abundance Expected (%)") + 
  theme(legend.position="none") 

tiff('./figures/accuracy_abundance2.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(accuracy_abundance2)
dev.off()
