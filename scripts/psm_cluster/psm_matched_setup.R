psm_setup_fun <- function(d){
  sessionInfo()
  wd <- "~/selection"
  library(rstan)
  library(brms)
  library(tidybayes)
  library(plyr)
  library(tidyverse)
  trait = "C"; outcome = "retired"; mod = "age"
  m <- if(mod == "SES") c("parEdu", "grsWages", "parOccPrstg")  else mod
  print(paste("m =", m))
  d <- if(mod %in% c("reliability", "predInt")){"none"} else mod
  print(paste("d =", d))
  
  # load data & sample model 
  load(sprintf("/scratch/edbeck/data/psm_combined/%s_%s_%s.RData", outcome, trait, d))
  
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
  
  save(fit2, file = "scratch/edbeck/data/matched_compiled_small.RData")
  # save(fit, file = sprintf("/scratch/edbeck/psm/models/matched_%s_%s_%s", trait, outcome, mod))
}