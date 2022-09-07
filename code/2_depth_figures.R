#Load Libraries
library(tidyverse)
library(purrr)
library(reshape2)
library(viridis)
library(ggpubr)
library(ggthemes)


#####################################################################################
###############################Sequencing Depth######################################


#Load Data

#add in sample details
sample_details = read_csv("./data/raw_data/sample_details.csv")

#add in sequencing stats
PE150 = read_csv("./data/raw_data/sequencing_stats_PE150.csv")
PE250 = read_csv("./data/raw_data/sequencing_stats_PE250.csv")

##PE150 Depth Data##

#read in file names for HiSeq (PE150)
PE150_files = fs::dir_ls("./data/raw_data/HiSeq/HiSeq_MC_Amp_03.15.2022", glob="*.depth")

#create a data frame with all files, named by their "source" file
depth_data_PE150 = PE150_files %>% 
  purrr::map_dfr(read_table, .id = "source", col_names = FALSE) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "./data/raw_data/HiSeq/HiSeq_MC_Amp_03.15.2022/", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, ".depth", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "-RBD-", "_Amp_")) 

#add in column names
names(depth_data_PE150) = c("sample_id", "ref_strain", "variant", "nucleotide", "depth")


##PE250 Depth Data##

#read in file names
PE250_files = fs::dir_ls("./data/raw_data/MiSeq/MiSeq_MC_Amp_02.14.2022", glob="*.depth")

#create a data frame with all files, named by their "source" file
depth_data_PE250 = PE250_files %>% 
  purrr::map_dfr(read_table, .id = "source", col_names = FALSE) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "./data/raw_data/MiSeq/MiSeq_MC_Amp_02.14.2022/", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, ".depth", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "-RBD-", "_Amp_")) 


#add column names
names(depth_data_PE250) = c("sample_id", "ref_strain", "variant", "nucleotide", "depth")

#add in sample data 
depth_data_PE150 = depth_data_PE150 %>% 
  left_join(sample_details, by = "sample_id") %>% 
  mutate(read_length = "PE150")

depth_data_PE250 = depth_data_PE250 %>% 
  left_join(sample_details, by = "sample_id") %>% 
  mutate(read_length = "PE250")

depth_data = bind_rows(depth_data_PE250, depth_data_PE150)

PE150$read_length = "PE150"
PE250$read_length = "PE250"

stats = bind_rows(PE150, PE250)


##################################################################
################DATA VISUALIZATION################################


#Let's start with the controls. 

twist = depth_data %>% 
  filter(community == "twist") %>% 
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  facet_wrap(primer_set~read_length, ncol = 2) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylab("Log10(Depth)") + 
  xlab("Position") + 
  theme_few() 

#tiff('./figures/twist_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(twist)
#dev.off()

#SARS-CoV-2 in PBS

heat_sars_pbs = depth_data %>% 
  filter(community == "heat_sars_pbs") %>% 
  ggplot(aes(x = nucleotide, y = log10(depth))) + 
  geom_line() + 
  facet_wrap(primer_set~read_length, ncol = 2) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylab("Log10(Depth)") + 
  xlab("Position") + 
  theme_few() 

#tiff('./figures/sars_pbs_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(heat_sars_pbs)
#dev.off()

#SARS-CoV-2 in Wastewater

heat_sars_wastewater = depth_data %>% 
  filter(community == "heat_sars_wastewater") %>% 
  ggplot(aes(x = nucleotide, y = log10(depth))) + 
  geom_line() + 
  facet_wrap(primer_set~read_length, ncol = 2) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylab("Log10(Depth)") + 
  xlab("Position") + 
  theme_few() 

#tiff('./figures/sars_ww_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(heat_sars_wastewater)
#dev.off()



##############Samples############


#How many nucleotides > 20X? 
depth_20X = depth_data %>% 
  filter(depth > 19) %>% 
  group_by(sample_id, read_length) %>% 
  tally() %>% 
  mutate(percent_genome_20X = n/29903*100)  %>% 
  select(-n) %>% 
  right_join(stats, by = c("sample_id", "read_length")) %>% 
  mutate(percent_genome_20X = replace_na(percent_genome_20X, 0)) %>% 
  left_join(sample_details %>% select(-details), by = c("sample_id"))


###Annotations###

##Midnight

depth_20X_Midnight = 
  depth_20X %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>%
  filter(primer_set == "Midnight")



##VarSkip

depth_20X_VarSkip = 
  depth_20X %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>%
  filter(primer_set == "VarSkip")


##V4

depth_20X_V4 = 
  depth_20X %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>%
  filter(primer_set == "V4")



#############Figures###################

midnight_depth_PE150 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "Midnight") %>%
  filter(read_length == "PE150") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_Midnight %>% filter(read_length == "PE150"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  ylim(0,5) + 
  theme_few()


midnight_depth_PE250 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "Midnight") %>%
  filter(read_length == "PE250") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_Midnight %>% filter(read_length == "PE250"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylim(0,5) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  theme_few()

midnight_depth = ggarrange(midnight_depth_PE150, 
                           midnight_depth_PE250, 
                           ncol = 1, 
                           align = "hv", 
                           labels = c("PE150", "PE250"))

#tiff('./figures/midnight_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(midnight_depth)
#dev.off()

##V4##

V4_depth_PE150 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "V4") %>%
  filter(read_length == "PE150") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_V4 %>% filter(read_length == "PE150"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylim(0,5) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  theme_few()


V4_depth_PE250 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "V4") %>%
  filter(read_length == "PE250") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_V4 %>% filter(read_length == "PE250"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylim(0,5) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  theme_few()


V4_depth = ggarrange(V4_depth_PE150, 
                     V4_depth_PE250, 
                           ncol = 1, 
                           align = "hv", 
                           labels = c("PE150", "PE250"))

#tiff('./figures/V4_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(V4_depth)
#dev.off()




#VarSkip

VarSkip_depth_PE150 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "VarSkip") %>%
  filter(read_length == "PE150") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_VarSkip %>% filter(read_length == "PE150"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylim(0,5) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  theme_few()


VarSkip_depth_PE250 = 
  depth_data %>% 
  filter(community == "1" | 
           community == "2" | 
           community == "3" | 
           community == "4" |
           community == "5") %>% 
  filter(primer_set == "VarSkip") %>%
  filter(read_length == "PE250") %>%
  mutate(depth = ifelse(depth == 0, 1, depth)) %>%
  ggplot() + 
  geom_line(aes(x = nucleotide, y = log10(depth))) + 
  geom_text(data = depth_20X_VarSkip %>% filter(read_length == "PE250"), 
            aes(x = 15000, y = 4.5, label = paste0(round(percent_genome_20X), "%")), 
            color = "#2297E6") + 
  facet_grid(dilution ~ community) + 
  geom_hline(yintercept = 1.30, linetype='dashed', color = "#2297E6", size = 1.1) + 
  ylim(0,5) + 
  xlab("Position") + 
  ylab ("Log10 (Depth)") + 
  theme_few()

VarSkip_depth = ggarrange(VarSkip_depth_PE150, 
                     VarSkip_depth_PE250, 
                     ncol = 1, 
                     align = "hv", 
                     labels = c("PE150", "PE250"))

#tiff('./figures/VarSkip_depth.tiff', units="in", width = 9, height = 8, res=600, compression = 'lzw', pointsize = 12)
plot(VarSkip_depth)
#dev.off()
