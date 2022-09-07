#Load Libraries
library(tidyverse)
library(magrittr)
library(ggstatsplot)
library(purrr)
library(reshape2)
library(viridis)
library(ggpubr)


##########################################################################
#####################Sequencing Statistics################################

#Load Data
PE150 = read_csv("./data/raw_data/sequencing_stats_PE150.csv")
PE250 = read_csv("./data/raw_data/sequencing_stats_PE250.csv")
sample_details = read_csv("./data/raw_data/sample_details.csv")

#Summary Statistics
PE150 %<>% 
  left_join(sample_details, by = "sample_id") %>% 
  select(-details.x) %>% 
  mutate(read_length = "PE150")

PE250 %<>% 
  left_join(sample_details, by = "sample_id") %>% 
  select(-details.x) %>% 
  mutate(read_length = "PE250")


#Combine all data 
my.data = rbind(PE150, PE250)

my.data %<>%
  mutate(filtered_raw_reads = filtered_reads/raw_reads*100) %>% 
  mutate(assigned_filtered_reads = reads_mapped_kallisto/filtered_reads*100) %>% 
  mutate(assigned_raw_reads = reads_mapped_kallisto/raw_reads*100) %>%
  mutate(details = details.y) %>% 
  select(-details.y) %>%
  mutate(
    conc_cp_ul = ifelse(dilution == "10^3", 1000,
                        ifelse(dilution == "10^2", 100, 
                               ifelse(dilution == "10^1", 10, "NA" 
                               ))))


#####################################################################################
###############################Sequencing Depth######################################


#Load Data

##PE150 Data##

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


##PE250 Data##

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


###SUMMARY STATISTICS###

#How much of the genome is sequenced?
genome_coverage_PE150 = depth_data_PE150 %>% 
  filter(depth > 0) %>% 
  group_by(details) %>% 
  tally() %>% 
  mutate(percent_genome = n/29903*100) %>% 
  select(-n) %>% 
  right_join(sample_details, by = "details") %>% 
  mutate(read_length = "PE150") %>% 
  mutate(percent_genome = replace_na(percent_genome, 0))

genome_coverage_PE250 = depth_data_PE250 %>% 
  filter(depth > 0) %>% 
  group_by(details) %>% 
  tally() %>% 
  mutate(percent_genome = n/29903*100) %>% 
  select(-n) %>% 
  right_join(sample_details, by = "details") %>% 
  mutate(read_length = "PE250") %>% 
  mutate(percent_genome = replace_na(percent_genome, 0))

genome_coverage = rbind(genome_coverage_PE150, genome_coverage_PE250)


#What is the average depth of coverage for each sample? 
depth_avg_PE150 = depth_data_PE150 %>% 
  group_by(sample_id, read_length) %>% 
  summarise(depth_avg = mean(depth), 
            depth_sd = sd(depth))

depth_avg_PE250 = depth_data_PE250 %>% 
  group_by(sample_id, read_length) %>% 
  summarise(depth_avg = mean(depth), 
            depth_sd = sd(depth))


depth_avg = rbind(depth_avg_PE150, depth_avg_PE250)


#How much of the genome is covered with >20X depth? 
depth_20X_PE150 = depth_data_PE150 %>% 
  filter(depth > 19) %>% 
  group_by(details) %>% 
  tally() %>% 
  mutate(percent_genome_20X = n/29903*100)  %>% 
  select(-n) %>% 
  right_join(sample_details, by = "details") %>% 
  mutate(read_length = "PE150") %>% 
  mutate(percent_genome_20X = replace_na(percent_genome_20X, 0))


depth_20X_PE250 = depth_data_PE250 %>% 
  filter(depth > 19) %>% 
  group_by(details) %>% 
  tally() %>% 
  mutate(percent_genome_20X = n/29903*100)  %>% 
  select(-n) %>% 
  right_join(sample_details, by = "details") %>% 
  mutate(read_length = "PE250") %>% 
  mutate(percent_genome_20X = replace_na(percent_genome_20X, 0))


depth_20X = rbind(depth_20X_PE150, depth_20X_PE250)


#Join to my.data 
my.data %<>% 
left_join(genome_coverage, by = c("sample_id", "details", "read_length", "primer_set","community","dilution")) %>%
left_join(depth_avg, by = c("sample_id", "read_length")) %>%
left_join(depth_20X, by = c("sample_id", "details", "read_length", "primer_set","community","dilution")) 


###########################################################################################
########################Abundance Estimates################################################



#Load Data


##Reference Genomes##
reference_genomes = read_tsv("./data/raw_data/reference_genomes.tsv")

#Add WHO Label (Lineage Name) 
reference_genomes %<>% 
  mutate(lineage_name = case_when(
    .$pangolin_lineage == "B.1.617.2" ~ "Delta",
    .$pangolin_lineage == "B" ~ "Wuhan",
    .$pangolin_lineage == "B.1.1.7" ~ "Alpha",
    .$pangolin_lineage == "B.1.351" ~ "Beta",
    .$pangolin_lineage == "P.1" ~ "Gamma",
    .$pangolin_lineage == "B.1.427" | .$pangolin_lineage == "B.1.429" ~ "Epsilon",
    .$pangolin_lineage == "B.1.1.529" ~ "Omicron",
    .$pangolin_lineage == "B.1.526" ~ "Iota",
    .$pangolin_lineage == "B.1.525" ~ "Eta",
    .$pangolin_lineage == "B.1.617.1" ~ "Kappa",
    .$pangolin_lineage == "B.1.621" | .$pangolin_lineage == "B.1.621.1" ~ "Mu",
    .$pangolin_lineage == "P.2" ~ "Zeta",
    .$pangolin_lineage == "B.1.617.3" ~ "Unamed", 
    TRUE ~ "Other"
  ))

#Update to match kallisto output 
reference_genomes %<>% 
  mutate(date = as.Date(date, "%m/%d/%y")) %>% 
  mutate(date = as.character(date)) %>%
  mutate(date = case_when(
    .$strain == "hCoV-19/USA/FL-STP-Saliva-412-VAC/2021" ~ "2021-02", 
    .$strain == "hCoV-19/USA/FL-Shands-VTM-53/2021" ~ "2021-01", 
    .$strain == "hCoV-19/USA/FL-Path-Saliva-872/2021" ~ "2021-01",
    TRUE ~ date)) %>%
  unite(col = variant_id, strain, gisaid_epi_isl, date, sep = "|") %>%
  mutate(variant_id = case_when(variant_id == "Junk|OL672836.1|NA" ~ "Junk/betacoronavirus/OL672836.1", 
                                variant_id == "Junk|OL717063.1|NA" ~ "Junk/betacoronavirus/OL717063.1", 
                                TRUE ~ variant_id))



#Load expected proportions
variant_props_mc = read_csv("./data/raw_data/variant_props_mc.csv")

##Convert "extected prop" to expected tpm 
variant_props_mc %<>%
  mutate(tpm_expected = abundance_expected/100*1e6)



##Load in kallisto output - PE150 Data##


#read in file names for HiSeq PE150 Run
PE150_files = fs::dir_ls("./data/raw_data/HiSeq/abundance_PE150", glob="*.tsv")

#create a data frame with all files, named by their "source" file
abundance_PE150 = 
  PE150_files %>% 
  purrr::map_dfr(read_table, .id = "source", col_names = FALSE) %>%
  dplyr::mutate(source = stringr::str_replace(source, "./data/raw_data/HiSeq/abundance_PE150/", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "_abundance.tsv", "")) %>%
  dplyr::mutate(source = stringr::str_replace(source, "-Amp-", "_Amp_")) 


names(abundance_PE150) = c("sample_id", "variant_id", "length", "eff_length", "est_counts", "tpm")

abundance_PE150 %<>% 
  mutate(est_counts = as.numeric(est_counts)) %>%
  drop_na() %>% 
  left_join(
    reference_genomes %>%
      select(variant_id, lineage_name)) %>%
  left_join(
    sample_details %>% 
      select(sample_id, details)) %>%
  mutate(tpm = as.numeric(tpm)) 






##Load Kallisto output - PE250 Data##

#read in file names for MiSeq PE250 Run
PE250_files = fs::dir_ls("./data/raw_data/MiSeq/abundance_PE250", glob="*.tsv")

#create a data frame with all files, named by their "source" file
abundance_PE250 = 
  PE250_files %>% 
  purrr::map_dfr(read_table, .id = "source", col_names = FALSE) %>%
  dplyr::mutate(source = stringr::str_replace(source, "./data/raw_data/MiSeq/abundance_PE250/", "")) %>% 
  dplyr::mutate(source = stringr::str_replace(source, "_abundance.tsv", "")) %>%
  dplyr::mutate(source = stringr::str_replace(source, "-Amp-", "_Amp_")) 


names(abundance_PE250) = c("sample_id", "variant_id", "length", "eff_length", "est_counts", "tpm")

abundance_PE250 %<>% 
  mutate(est_counts = as.numeric(est_counts)) %>%
  drop_na() %>% 
  left_join(
    reference_genomes %>%
      select(variant_id, lineage_name)) %>%
  left_join(
    sample_details %>% 
      select(sample_id, details)) %>%
  mutate(tpm = as.numeric(tpm))



###########################################################################

#Filtering Abundance estimates based on 0.1% threshold as described in Bajeens et al. 
#0.1% == 1,000 tpm
#Sum tpm per lineage

abundance_PE150 %<>% 
  mutate(tpm = case_when(tpm < 1000 ~ 0, TRUE ~ tpm)) %>%
  group_by(sample_id, lineage_name) %>% 
  summarise(tpm_observed = sum(tpm)) %>% 
  left_join(sample_details) %>% 
  left_join(variant_props_mc) %>% 
  mutate(abundance_expected = replace_na(abundance_expected, 0)) %>%
  mutate(tpm_expected = replace_na(tpm_expected, 0)) %>%
  mutate(abundance_observed = tpm_observed/1e6*100) %>% 
  mutate(read_length = "PE150")


abundance_PE250 %<>% 
  mutate(tpm = case_when(tpm < 1000 ~ 0, TRUE ~ tpm)) %>%
  group_by(sample_id, lineage_name) %>% 
  summarise(tpm_observed = sum(tpm)) %>% 
  left_join(sample_details) %>% 
  left_join(variant_props_mc) %>% 
  mutate(abundance_expected = replace_na(abundance_expected, 0)) %>%
  mutate(tpm_expected = replace_na(tpm_expected, 0)) %>%
  mutate(abundance_observed = tpm_observed/1e6*100) %>% 
  mutate(read_length = "PE250")


abundance_estimates = rbind(abundance_PE150, abundance_PE250)
write_rds(abundance_estimates, "./data/processed_data/abundance.estimates.RDS")

abundance_estimates_wide = 
  abundance_estimates %>%
  select(sample_id, primer_set, community, dilution, read_length, 
         lineage_name, abundance_expected, abundance_observed) %>% 
  mutate(error = abundance_observed - abundance_expected) %>%
  pivot_wider(names_from = lineage_name, values_from = c(error,abundance_expected, abundance_observed))


##############################################################################

RMSE_Sample = 
  abundance_estimates %>% 
  group_by(sample_id, read_length) %>% 
  summarize(RMSE = sqrt(mean((abundance_observed - abundance_expected)^2)))


SSR_Sample = 
  abundance_estimates %>% 
  mutate(squares = (abundance_expected - abundance_observed)^2) %>% 
  group_by(sample_id, read_length) %>% 
  summarize(SSR = sum(squares))

SST_Sample = 
  abundance_estimates %>% 
  group_by(sample_id, read_length) %>% 
  summarize(abundance_observed_mean = mean(abundance_observed)) %>% 
  right_join(abundance_estimates) %>% 
  mutate(squares = (abundance_observed - abundance_observed_mean)^2) %>% 
  group_by(sample_id, read_length) %>%
  summarise(SST = sum(squares))

R_2_Sample = 
  left_join(SSR_Sample, SST_Sample) %>% 
  mutate(R_2 = (1 - (SSR/SST))) %>% 
  select(sample_id, read_length, R_2)

#########################################################################################
#Estimate error for each sequencing approach, overall.

abundance_estimates_mc = abundance_estimates %>% filter(dilution != "NA")

RMSE_seq_approach =
  abundance_estimates_mc %>% 
  unite(col = seq_approach, primer_set, read_length, sep = "_") %>%
  group_by(seq_approach) %>% 
  summarize(RMSE = sqrt(mean((abundance_observed - abundance_expected)^2)))


SSR_seq_approach = 
  abundance_estimates_mc %>% 
  mutate(squares = (abundance_expected - abundance_observed)^2) %>% 
  unite(col = seq_approach, primer_set, read_length, sep = "_") %>%
  group_by(seq_approach) %>% 
  summarize(SSR = sum(squares))

SST_seq_approach = 
  abundance_estimates_mc %>% 
  group_by(primer_set, read_length) %>% 
  summarize(abundance_observed_mean = mean(abundance_observed)) %>% 
  right_join(abundance_estimates_mc) %>% 
  mutate(squares = (abundance_observed - abundance_observed_mean)^2) %>% 
  group_by(primer_set, read_length) %>% 
  summarise(SST = sum(squares)) %>% 
  unite(col = seq_approach, primer_set, read_length, sep = "_") 

R_2_seq_approach = 
  left_join(SSR_seq_approach, SST_seq_approach) %>% 
  mutate(R_2 = (1 - (SSR/SST))) %>%
  select(seq_approach, R_2)

seq_approach_error = left_join(RMSE_seq_approach, R_2_seq_approach)
write.csv(seq_approach_error, "./data/processed_data/seq_approach_error.csv")

################################################################################################

#Estimate error for each primer set.

RMSE_primer_set =
  abundance_estimates_mc %>% 
  group_by(primer_set) %>% 
  summarize(RMSE = sqrt(mean((abundance_observed - abundance_expected)^2)))


SSR_primer_set = 
  abundance_estimates_mc %>% 
  mutate(squares = (abundance_expected - abundance_observed)^2) %>% 
  group_by(primer_set) %>% 
  summarize(SSR = sum(squares))

SST_primer_set = 
  abundance_estimates_mc %>% 
  group_by(primer_set) %>% 
  summarize(abundance_observed_mean = mean(abundance_observed)) %>% 
  right_join(abundance_estimates_mc) %>% 
  mutate(squares = (abundance_observed - abundance_observed_mean)^2) %>% 
  group_by(primer_set) %>% 
  summarise(SST = sum(squares)) 

R_2_primer_set = 
  left_join(SSR_primer_set, SST_primer_set) %>% 
  mutate(R_2 = (1 - (SSR/SST))) %>%
  select(primer_set, R_2)

primer_set_error = left_join(RMSE_primer_set, R_2_primer_set)
write.csv(primer_set_error, "./data/processed_data/primer_set_error.csv")




#Add in kallisto data to my.data df

my.data %<>% left_join(RMSE_Sample) %>% left_join(R_2_Sample) %>% left_join(abundance_estimates_wide)
write_rds(my.data, "./data/processed_data/my.data.RDS")