setwd("/scratch/edbeck")
library(rstan)
library(brms)
library(tidybayes)
library(plyr)
library(tidyverse)
sessionInfo()

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))
print(jobid)
args <- read.table("/scratch/edbeck/psm/psm_matched_args_pe1_24.txt", header = F, stringsAsFactors = F)[jobid,]
print(args)

# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, "Makevars")
# if (!file.exists(M)) file.create(M)
# cat("\nCXX14FLAGS=-O3 -std=c++14 -march=native -mtune=native -fPIC",
#     "CXX14=g++",
#     "CXX14= /export/gcc-5.5.0/bin/g++", # or clang++ but you may need a version postfix
#     "CXX14FLAGS=-std=c++14 -fPIC",
#     "CXX14FLAGS += -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION",
#     file = M, sep = "\n", append = TRUE)

def_hyp_fun <- function(study, mean, sd){
  h <- c(
    paste("p_value:", mod, "-", sd, " = 0", sep = ""),
    paste("p_value:", mod, " = 0", sep = ""),
    paste("p_value:", mod, "+", sd, " = 0", sep = "")
  )
  hypothesis(fit, h, class = "study")
}

# function to test hypotheses for moderators
hyp_fun <- function(fit, mod){
  fit$data %>% select(study, one_of(mod)) %>% 
    tbl_df() %>%
    group_by(study) %>%
    summarize_all(lst(mean, sd))
  print("hyp_fun")
  me <- 
    fixef(fit) %>% data.frame %>% mutate(term = rownames(.)) %>%
    filter(term %in% c("Intercept", "new.wave")) %>%
    mutate(term.type = mapvalues(term, unique(term), c("Intercept", "Slope")),
           term = "No Event") %>%
    setNames(c("b", "SE", "lower", "upper", "term", "term.type"))
  if(event != "none"){
    
    h <- hypothesis(fit, h, alpha = .05)
    me <- h$hypothesis %>% data.frame %>% 
      mutate(term = "Event", term.type = c("Slope", "Intercept")) %>% 
      dplyr::select(-Evid.Ratio, -Star) %>%
      setNames(c("b", "SE", "lower", "upper", "term", "term.type")) %>%
      full_join(me)
  }
  return(me)
}

fx_pred_fun <- function(fit, mod){
  cl <- class(fit$data[, mod])
  
  if(cl %in% c("numeric", "integer")){
    msd <- fit$data %>% select(one_of(mod)) %>% tbl_df() %>%
      summarize_all(lst(mean, sd))
    
    pred.fx <- crossing(
      p_value = seq(0,10,.25),
      mod_value = with(msd, c(mean-sd, mean,mean+sd))
    ) %>%
      setNames(c("p_value", mod))
  } else {
    pred.fx <- crossing(
      p_value = seq(0,10,.25),
      mod_value = factor(levels(fit$data[, mod]))
    ) %>%
      setNames(c("p_value", mod))
  }
  
  pred.fx <- bind_cols(
    pred.fx, 
    fitted(fit, newdata = pred.fx, probs = c(0.055, 0.945), re_formula = NA) %>% data.frame
  ) %>%
    mutate_at(vars(Estimate, Q5.5, Q94.5), lst(OR = inv_logit_scaled))
}

crossing_fun <- function(df, mod){
  pred.rx <- crossing(
    p_value = seq(0,10,.25),
    mod_value = with(df, c(mean-sd, mean,mean+sd))
  ) %>%
    setNames(c("p_value", mod))
  return(pred.rx)
}

rx_pred_fun <- function(fit, mod){
  studies <- unique(fit$data$study)
  cl <- class(fit$data[, mod])
  
  if(cl %in% c("numeric", "integer")){
    pred.rx <- fit$data %>% select(study, one_of(mod)) %>% 
      tbl_df() %>%
      group_by(study) %>%
      summarize_all(lst(mean, sd)) %>%
      group_by(study) %>%
      nest() %>%
      ungroup() %>%
      mutate(pred = map2(data, mod, crossing_fun)) %>%
      unnest(pred) %>%
      select(-data)
  } else {
    pred.rx <- crossing(
      study = studies,
      p_value = seq(0,10,.25),
      mod_value = factor(levels(fit$data[, mod]))
    ) %>%
      setNames(c("study", "p_value", mod))
  }
  
  pred.rx <- bind_cols(
    pred.rx, 
    fitted(fit, newdata = pred.rx, probs = c(0.055, 0.945)) %>% data.frame
  ) %>%
    mutate_at(vars(Estimate, Q5.5, Q94.5), lst(OR = inv_logit_scaled))
  return(pred.rx)
}

brms_matched_fun <- function(trait, outcome, mod, chain, ncores){
  print(sprintf("outcome = %s, trait = %s, mod = %s, chain = %s", outcome, trait, mod, chain))
  
  # setup
  m <- if(mod == "SES") c("parEdu", "grsWages", "parOccPrstg")  else mod
  print(paste("m =", m))
  d <- if(mod %in% c("reliability", "predInt")){"none"} else mod
  print(paste("d =", d))
  
  # load data & sample model 
  load(sprintf("/scratch/edbeck/data/psm_combined/%s_%s_%s.RData", outcome, trait, d))
  load(sprintf("/scratch/edbeck/data/matched_compiled_small.RData"))
  
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
  else if(mod %in% c("reliability", "predInt")){f <- formula(paste("o_value ~ p_value + ", paste("p_value*", m, collapse = " + "), " + (p_value | study)", sep = ""))}
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
           , cores = ncores
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
  save(fit, file = sprintf("/scratch/edbeck/psm/matched/models/matched_%s_%s_%s_%s.RData", trait, outcome, mod, chain))
  # save(fx, rx, file = sprintf("/scratch/edbeck/psm/matched/summary/matched_%s_%s_%s", trait, outcome, mod))
  # save(fx.draws, tau.draws, file = sprintf("/scratch/edbeck/psm/matched/draws/matched_%s_%s_%s", trait, outcome, mod))
  # rm(c("fit", "fx", "rx", "fx.draws", "rx.draws", "df"))
  rm(list = c("fit", "fit2", "df_l", "d1"))
  gc()
}

brms_matched_fun(args[,1], args[,2], args[,3], args[,4], args[,5])