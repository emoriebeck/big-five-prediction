sessionInfo()
wd <- "~/selection"
library(rstan)
library(brms)
library(tidybayes)
library(plyr)
library(tidyverse)

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))
print(jobid)
args <- read.table("~/psm_unmatched_args.txt", header = F)[jobid,]
print(args)

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

brms_unmatched_fun <- function(trait, outcome, mod){
  print(paste("outcome =", outcome))
  print(paste("trait =", trait))
  print(paste("mod =", mod))
  # load data 
  m <- if(mod == "SES") c("parEdu", "grsWages", "parOccPrstg")  else mod
  print(paste("m =", m))
  d <- if(mod %in% c("reliability", "predInt")){"none"} else mod
  print(paste("d =", d))
  load(sprintf("/scratch/edbeck/data/psm_unmatched/%s_%s.RData", outcome, trait))
  
  df <- df %>% mutate(o_value = as.numeric(as.character(o_value)))
  
  # set priors
  Prior <-  c(set_prior("cauchy(0,1)", class = "sd"))
  # formula 
  if(mod == "none"){f <- formula(o_value | trials(1) ~ p_value + (p_value | study))}
  else {f <- formula(paste("o_value | trials(1) ~ p_value + ", paste("p_value*", m, collapse = " + "), "+ (", paste("p_value*", m, collapse = " + "), " | study)", sep = ""))}
  
  # run the model 
  Iter <- 2000; Warmup <- 1000; treedepth <- 20
  start.tmp <- Sys.time()
  fit <- brm(formula = f
             , data = df
             , prior = Prior
             , iter = Iter
             , warmup = Warmup
             , family = binomial(link = "logit")
             , control = list(adapt_delta = 0.99, max_treedepth = treedepth))
             # , cores = 4)
  print(Sys.time() - start.tmp)
  
  # extract key parameters
  # fixed effects
  fx <- fixef(fit) %>% data.frame %>% 
    rownames_to_column("names") %>%
    mutate_at(vars(-names), lst(OR = inv_logit_scaled)) %>%
    tbl_df 
  # random effects
  rx <- ranef(fit)[[1]] %>% array_tree(3) %>% 
    tibble(names = names(.), data = .) %>% 
    mutate(data = map(data, ~(.) %>% data.frame %>% 
                        rownames_to_column("study"))) %>% 
    unnest(data) %>%
    mutate_at(vars(-names, -study), lst(OR = inv_logit_scaled)) %>%
    tbl_df
  
  # samples
  fx.draws <- fit %>% tidy_draws() %>% 
    select(.chain:.draw, matches("^b_"), matches("p_value]$")) %>%
    mutate_at(vars(matches("p_value]$")), ~(.) + b_p_value) %>%
    gather(key = item, value = s_value, -(.chain:.draw))
  
  tau.draws <- fit %>% tidy_draws() %>% 
    select(.chain:.draw, matches("^sd"), matches("^cor"))
  
  if(mod != "none"){
    pred.fx <- fx_pred_fun(fit, m)
    pred.rx <- rx_pred_fun(fit, m)
    save(pred.fx, pred.rx, file = sprintf("/scratch/edbeck/psm/predicted/unnmatched_pred_%s_%s_%s", trait, outcome, mod))
  }
  
  save(fit, file = sprintf("/scratch/edbeck/psm/models/unmatched_%s_%s_%s", trait, outcome, mod))
  save(fx, rx, file = sprintf("/scratch/edbeck/psm/summary/unmatched_%s_%s_%s", trait, outcome, mod))
  save(fx.draws, tau.draws, file = sprintf("/scratch/edbeck/psm/draws/unmatched_%s_%s_%s", trait, outcome, mod))
  rm(list("fit", "fx", "rx", "fx.draws", "rx.draws", "df"))
  gc()
}

brms_unmatched_fun(args[,1], args[,2], args[,3])