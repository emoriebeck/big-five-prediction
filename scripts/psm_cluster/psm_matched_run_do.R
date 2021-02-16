dir.create("~/psm/models", recursive = T)
dir.create("~/psm/psm_combined")

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC",
    "CXX14=g++", # or clang++ but you may need a version postfix
    file = M, sep = "\n", append = TRUE)

install.packages("brms")

# setwd("/scratch/edbeck")
library(rstan)
library(brms)
# library(tidybayes)
library(plyr)
library(tidyverse)
sessionInfo()

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC",
    "CXX14=g++", # or clang++ but you may need a version postfix
    file = M, sep = "\n", append = TRUE)

# jobid = as.integer(Sys.getenv("PBS_ARRAYID"))
# print(jobid)
args <- read.table("~/psm_matched_args_old.txt", header = F, stringsAsFactors = F)

# print(args)

trait <- args[1,1]; outcome <- args[1,2]; mod <- args[1,3]; chain <- args[1,4]
m <- if(mod == "SES") c("parEdu", "grsWages", "parOccPrstg")  else mod
print(paste("m =", m))
d <- if(mod %in% c("reliability", "predInt")){"none"} else mod
print(paste("d =", d))

# load data & sample model 
load(sprintf("~/psm/psm_combined/%s_%s_%s.RData", outcome, trait, d))

# clean data & keep only needed columns and a subset of the used variables
d1 <- df_l[[1]] %>%
  group_by(study, o_value) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = map(data, ~(.) %>% filter(row_number() %in% sample(1:nrow(.), 50, replace = F)))) %>%
  unnest(data) %>%
  select(study, SID, p_value, o_value, one_of(d)) %>%
  filter(complete.cases(.))

# set priors & model specifications 
Prior <-  c(set_prior("cauchy(0,1)", class = "sd"),
            set_prior("student_t(3, 0, 2)", class = "b"),
            set_prior("student_t(3, 0, 5)", class = "Intercept"))
Iter <- 30; Warmup <- 21; treedepth <- 20
f <- formula(paste("o_value ~ p_value + ", paste("p_value*", m, collapse = " + "), "+ (", paste("p_value*", m, collapse = " + "), " | study)", sep = ""))
fit2 <- brm(formula = f
            # , data = df_l
            , data = d1
            , prior = Prior
            , iter = Iter
            , warmup = Warmup
            , family = bernoulli(link = "logit")
            # , control = list(adapt_delta = 0.99, max_treedepth = treedepth)
            , cores = 4)

save(fit2, file = "~/matched_compiled_small.RData")
# brms_matched_fun <- function(trait, outcome, mod, chain, ncores){
rm(list = ls())

brms_matched_fun <- function(i){
  trait <- args[i,1]; outcome <- args[i,2]; mod <- args[i,3]; chain <- args[i,4]
  print(sprintf("outcome = %s, trait = %s, mod = %s, chain = %s", outcome, trait, mod, chain))
  
  # setup
  m <- if(mod == "SES") c("parEdu", "grsWages", "parOccPrstg")  else mod
  print(paste("m =", m))
  d <- if(mod %in% c("reliability", "predInt")){"none"} else mod
  print(paste("d =", d))
  
  # load data & sample model 
  load(sprintf("~/psm/psm_combined/%s_%s_%s.RData", outcome, trait, d))
  load(sprintf("~/matched_compiled_small.RData"))
  
  # clean data to keep only needed columns
  df_l <- map(df_l, ~(.) %>%
                # comment out here
                # group_by(study, o_value) %>%
                # nest() %>%
                # ungroup() %>%
                # mutate(data = map(data, ~(.) %>% filter(row_number() %in% sample(1:nrow(.), 50, replace = F)))) %>%
                # unnest(data) %>%
                # to here
                select(study, SID, p_value, o_value, one_of(m)) %>%
                filter(complete.cases(.)))
  d1 <- df_l[chain]
  
  # formula 
  if(mod == "none"){f <- formula(o_value ~ p_value + (p_value | study))}
  else {f <- formula(paste("o_value ~ p_value + ", paste("p_value*", m, collapse = " + "), "+ (", paste("p_value*", m, collapse = " + "), " | study)", sep = ""))}
  
  Prior <-  c(set_prior("cauchy(0,1)", class = "sd"),
              set_prior("student_t(3, 0, 2)", class = "b"),
              set_prior("student_t(3, 0, 5)", class = "Intercept"))
  Iter <- 2000; Warmup <- 1000; treedepth <- 20
  
  # run the models using update and previously compiled C++ stan code
  start.tmp <- Sys.time()
  # plan(multiprocess)
  # fit <- future_map(df_l, function(x){
  fit <- #map(df_l, function(x){
    # tmp <- 
    update(fit2
           , formula = f
           , newdata = d1#x
           , iter = Iter
           , warmup = Warmup
           , cores = 4
           )
    # class(fit) <- c("brmsfit_multiple", class(fit))
    # return(tmp)
  # })
  # }, .progress = T)
  print(end.tmp <- Sys.time() - start.tmp)
  
  # combine models 
  # rhats <- map_df(fit, function(x)data.frame(as.list(rhat(x))))
  # fit <- combine_models(mlist = fit, check_data = FALSE)
  # fit$data.name <- "df_l"
  # fit$rhats <- rhats
  # class(fit) <- c("brmsfit_multiple", class(fit))
  
  # fit2 <- brm_multiple(formula = f
  #                      # , data = df_l
  #                      , data = d1
  #                      , prior = Prior
  #                      , iter = Iter
  #                      , warmup = Warmup
  #                      , family = bernoulli(link = "logit")
  #                      , control = list(adapt_delta = 0.99, max_treedepth = treedepth))
  #                      # , cores = 4)
  
  # extract key parameters
  # fixed effects
  # fx <- fixef(fit, probs = c(0.055, 0.945)) %>% data.frame %>% 
  #   rownames_to_column("names") %>%
  #   mutate_at(vars(-names), lst(OR = inv_logit_scaled)) %>%
  #   tbl_df 
  # # random effects
  # rx <- ranef(fit, probs = c(0.055, 0.945))[[1]] %>% array_tree(3) %>% 
  #   tibble(names = names(.), data = .) %>% 
  #   mutate(data = map(data, ~(.) %>% data.frame %>% 
  #                       rownames_to_column("study"))) %>% 
  #   unnest(data) %>%
  #   mutate_at(vars(-names, -study), lst(OR = inv_logit_scaled)) %>%
  #   tbl_df
  # 
  # # samples
  # fx.draws <- fit %>% tidy_draws() %>% 
  #   select(.chain:.draw, matches("^b_"), matches("p_value]$")) %>%
  #   mutate_at(vars(matches("p_value]$")), ~(.) + b_p_value) %>%
  #   gather(key = item, value = s_value, -(.chain:.draw))
  # 
  # tau.draws <- fit %>% tidy_draws() %>% 
  #   select(.chain:.draw, matches("^sd"), matches("^cor"))
  # 
  # if(mod != "none"){
  #   pred.fx <- fx_pred_fun(fit, m)
  #   pred.rx <- rx_pred_fun(fit, m)
  #   save(pred.fx, pred.rx, file = sprintf("/scratch/edbeck/psm/matched/predicted/matched_pred_%s_%s_%s", trait, outcome, mod))
  #   rm(c("pred.fx", "pred.rx"))
  # }
  # 
  save(fit, file = sprintf("~/psm/models/matched_%s_%s_%s_%s.RData", trait, outcome, mod, chain))
  # save(fx, rx, file = sprintf("/scratch/edbeck/psm/matched/summary/matched_%s_%s_%s", trait, outcome, mod))
  # save(fx.draws, tau.draws, file = sprintf("/scratch/edbeck/psm/matched/draws/matched_%s_%s_%s", trait, outcome, mod))
  # rm(c("fit", "fx", "rx", "fx.draws", "rx.draws", "df"))
  rm(list = c("fit", "fit2", "df_l", "d1"))
  gc()
}

map(1:nrow(args), brms_matched_fun)
# brms_matched_fun(args[,1], args[,2], args[,3], args[,4], args[,5])