################################ Libraries
library(tidyverse)
library(tidymodels)
library(xgboost)

# Read data
tr <- read_csv("mouse_tr.csv") %>%
  mutate(celltype = factor(celltype))

ts <- read_csv("mouse_ts_mask.csv")


################################ remove outliers through mahalanobis distance
tr_mat <- as.matrix(tr[, 2:(ncol(tr) - 1)])

# covariance matrix
cov_mat <- cov(tr[, 2:(ncol(tr) - 1)])

# get diagonal indices
cov_mat_diag <- row(cov_mat) - col(cov_mat) == 0

# add a very small amount to the diagonal to prevent it from becoming a singular matrix
cov_mat[cov_mat_diag] <- cov_mat[cov_mat_diag] + 1e-10

# add mahalanobis distance
tr <- tr %>% mutate(mah = mahalanobis(tr[, 2:(ncol(tr) - 1)],
  center = colMeans(tr_mat),
  cov = cov_mat
))

mah_hist <- ggplot(tr, aes(x = mah)) +
  geom_histogram()

# remove observations over 2,000 mahalanobis distance
tr_rem_out <- tr %>%
  filter(mah < 2000) %>%
  select(-mah)

# split tr into train/test sets
set.seed(1)
split <- initial_split(tr_rem_out, 0.8, strata = celltype)
tr_train <- training(split)
tr_test <- testing(split)


################################ boosting
# create grid of tree and learning rate values
boost_grid <- expand_grid(
  trees = seq(50, 200, 5),
  learn = c(0.2, 0.3, 0.4, 0.5, 0.6)
) %>%
  mutate(accuracy = 0)

# try different combinations of the values
for (i in 1:nrow(boost_grid)) {
  fit <- boost_tree(
    tree_depth = 1,
    trees = boost_grid[i, ]$trees,
    learn_rate = boost_grid[i, ]$learn
  ) %>%
    set_engine("xgboost") %>%
    set_mode("classification") %>%
    fit(celltype ~ ., tr_train[, -1])

  pred <- augment(fit, tr_test)

  boost_grid[i, "accuracy"] <- accuracy(pred, truth = celltype, estimate = .pred_class)$.estimate
}

# get minimum number of trees per accuracy (less variance)
max_acc_min_trees <- boost_grid %>%
  group_by(accuracy) %>%
  filter(trees == min(trees)) %>%
  arrange(-accuracy) %>%
  ungroup()

# best validation accuracy model
best_fit <- boost_tree(
  tree_depth = 1,
  trees = 110,
  learn_rate = 0.5
) %>%
  set_engine("xgboost") %>%
  set_mode("classification") %>%
  fit(celltype ~ ., tr_train[, -1])

best_pred <- augment(best_fit, tr_test)

# accuracy per class
best_pred_acc <- best_pred %>%
  count(celltype, .pred_class) %>%
  group_by(celltype) %>%
  mutate("accuracy(%)" = n / sum(n) * 100) %>%
  filter(celltype == .pred_class) %>%
  select(celltype, "accuracy(%)") %>%
  arrange(-`accuracy(%)`)

################################ fit on entire training set and predict test set
test_fit <- boost_tree(
  tree_depth = 1,
  trees = 135,
  learn_rate = 0.3
) %>%
  set_engine("xgboost") %>%
  set_mode("classification") %>%
  fit(celltype ~ ., tr_rem_out[, -1])

ts_boost <- augment(test_fit, ts) %>%
  select(location, .pred_class) %>%
  rename("celltype" = .pred_class)

write_csv(ts_boost, file = "mouse_mypred.csv")