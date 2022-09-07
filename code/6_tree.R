#Load libraries
library(tidyverse)

#Load data 

my.data = readRDS("./data/processed_data/my.data.RDS")

mc.data = my.data %>% filter(community == "1"|
                               community == "2"|
                               community == "3"|
                               community == "4"|
                               community == "5")


#Clean and tidy data 

mc.data = mc.data %>% select(R_2, 
                             filtered_reads, 
                             reads_mapped_kallisto, 
                             primer_set, 
                             community, 
                             read_length, 
                             conc_cp_ul, 
                             percent_genome, 
                             depth_avg, 
                             percent_genome_15X) %>% 
  mutate(community = as.factor(community)) %>%
  mutate_if(is.character, factor) 
                             


## Build a model 

library(tidymodels)

set.seed(123)
mc_split <- initial_split(mc.data, strata = R_2)
mc_train <- training(mc_split)
mc_test <- testing(mc_split)

set.seed(234)
mc_folds <- vfold_cv(mc_train, strata = R_2)
mc_folds


#Create tunable decision tree
tree_spec <- decision_tree(
  cost_complexity = tune(),
  tree_depth = tune(),
  min_n = tune()
) %>%
  set_engine("rpart") %>%
  set_mode("regression")

tree_spec

#Create a grid (set) of possible parameters to try
tree_grid <- grid_regular(cost_complexity(), tree_depth(), min_n(), levels = 4)

tree_grid

#Tune the model directly
doParallel::registerDoParallel()


set.seed(345)
tree_rs <- tune_grid(
  tree_spec,
  R_2 ~ .,
  resamples = mc_folds,
  grid = tree_grid,
  metrics = metric_set(rmse, rsq, mae, mape)
)

tree_rs

## Explore Results 

collect_metrics(tree_rs)

autoplot(tree_rs) + theme_light(base_family = "IBMPlexSans")

final_tree <- finalize_model(tree_spec, select_best(tree_rs, "rsq"))

#Finalize model
final_fit <- fit(final_tree, R_2 ~ ., mc_train)
final_rs <- last_fit(final_tree, R_2 ~ ., mc_split)


#What are the most important variables? 

library(vip)

final_fit %>%
  vip(geom = "col", aesthetics = list(fill = "midnightblue", alpha = 0.8)) +
  scale_y_continuous(expand = c(0, 0))

library(parttree)

ex_fit <- fit(
  final_tree,
  R_2 ~ log10(filtered_reads) + primer_set,
  mc_train
)

mc_train %>%
  ggplot(aes(log10(filtered_reads), primer_set)) +
  geom_parttree(data = ex_fit, aes(fill = R_2), alpha = 0.3) +
  geom_point(alpha = 0.7, width = 1, height = 0.5, aes(color = R_2)) +
  scale_colour_viridis_c(aesthetics = c("color", "fill"))

library(rpart.plot)  # for visualizing a decision tree

final_fit %>% 
  extract_fit_engine() %>%
  rpart.plot(roundint=FALSE)







multi.class.model <- rpart(R_2 ~ ., data = mc.data)
rpart.plot(multi.class.model)

