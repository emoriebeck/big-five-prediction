wd <- "~/selection"
library(lme4)
library(fixest)
library(haven)
library(readxl)
library(parallel)
library(future)
library(future.apply)
library(listenv)
library(plyr)
library(tidyverse)

source(sprintf("%s/mumin_fun.R", wd))

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))
print(jobid)
args <- read.table("~/sca_args.txt", header = F)[jobid,]
print(args[,4])

sample_fun <- function(i, dfx){
  (dfx %>% select(study, SID, p_value) %>%
     group_by(study) %>%
     mutate(p_value = sample(p_value, n(), replace = F)) %>%
     ungroup())$p_value 
}

perm_fun <- function(i, df3, f1, po){
  n <- df3$n
  po <- po[n]
  x <- df3 %>% select(-p_value) %>% mutate(p_value = po)
  start <- Sys.time()
  fit2  <- femlm(formula(f1), data = x, family = "logit")
  print(Sys.time() - start)
  par2 <- coef(fit2)["p_value"]
  ci2 <- as.matrix(confint(fit2, parm="p_value")[1,])
  # tidy2 <- tidy(fit2, effects = "fixed", conf.int = T) %>% 
  #   filter(term == "p_value") %>%
  #   select(estimate, conf.low, conf.high) %>%
  #   as.matrix()
  tidy2 <- c(par2, ci2)
  names(tidy2) <- c("p_value", "ci_lower", "ci_upper")
  rm(list = c("x", "df3", "fit2"))
  gc()
  return(tidy2)
}

sca_term_fun <- function(trait, outcome){
  terms <- c("p_value", (specifications %>% select(name, Effect, one_of(outcome)) %>% filter(complete.cases(.)) %>%
                           mutate(name = ifelse(Effect == "moderator", paste(name, ": p_value", sep = " "), name)))$name)
}

sca_formula_fun <- function(terms){
  f <- paste("o_value ~ ", # outcome 
             paste(terms, collapse = " + "), # fixed effects 
             "+ (p_value | study)", collapse = "") # random effects
}

sca_run_fun <- function(trait, outcome, min, max){
  load(sprintf("/scratch/edbeck/sca_workspace/sca_workspace_%s_%s.RData", outcome, trait))
  k <- 0L
  seq_comb <- seq(min, max, 1)
  # seq_comb <- seq_comb[1:2] # remove after testing
  
  res <- lapply(seq_comb, function(iComb){
    varComb <- iComb%%nVariants
    jComb <- (iComb - varComb)%/%nVariants
    if (varComb == 0L) {
      isok <- TRUE
      comb <- c(as.logical(intToBits(jComb)[comb.seq]), 
                comb.sfx)
      nvar <- sum(comb) - 1
      if (!formula_margin_check(comb, deps)) {
        isok <- FALSE
        res <- NA
      }
      new_terms <- allTerms[comb]
      if(sum(grepl(":", new_terms)) > 1){
        isok <- FALSE
        res <- NA
      }
    }
    if (!isok) {
      res <- NA
    } else{
      f1 <- paste("o_value ~ ", # outcome
                  paste(new_terms[new_terms != "(Intercept)"], collapse = " + "), # fixed effects
                  " | study", collapse = "")
      std <- df1 %>%
        group_by(study, o_value) %>%
        tally() %>%
        full_join(crossing(study = unique(.$study), o_value = c(0,1)))
      std <- unique(std$study[std$n < 50 | is.na(std$n)])
      df3 <- df1 %>%
        filter(!study %in% std) %>% 
        select(study:o_value, n, one_of(new_terms)) %>%
        mutate(o_value = as.numeric(as.character(o_value))) %>%
        filter(complete.cases(.)) 
      #   filter(study %in% std) 
      fit <- femlm(formula(f1), data = df3, family = "logit")
      mci1 <- confint(fit, parm = "p_value")
      mcoef1 <- matchCoef(fit, all.terms = gmCoefNames)
      ll1    <- logLik(fit)
      bic1 <- BIC(fit)
      nobs1  <- nobs(fit)
      psr1 <- fit$pseudo_r2
      row1   <- c(mcoef1[gmCoefNames], ci.lower = mci1[1,1], ci.upper = mci1[1,2],
                  bic = bic1, psr2 = psr1, ll = ll1)
      perm1 <- t(sapply(1:500, function(ii) perm_fun(ii, df3, f1, as.numeric(pval_out[,ii]))))
      print(Sys.time() - start)
      
      k <- k + 1L
      print(sprintf("%s %s %s", iComb, f1, k))
      res <- list(rval = row1, perm_res = perm1, ord = iComb)
      rm(list = c("row1", "perm1", "fit", "mci1", "mcoef1", "ll1", "nobs1"))
      gc()
    }
    return(res)
  })
  
  res <- res[!is.na(res)]
  
  rval <- ldply(res, `[[`, 1)
  perm_res <- llply(res, `[[`, 2)
  ord <- laply(res, `[[`, 3)
  perm_res <- abind::abind(perm_res, along = 3)
  
  # colnames(rval) <- c(gmCoefNames, "ci.lower", "ci.upper", "df", "ll")
  row.names(rval) <- ord
  # rval[, seq_along(gmCoefNames)] <- rval[, v <- order(termsOrder)]
  save(rval, file = sprintf("/scratch/edbeck/raw/raw_%s_%s_%s_%s.RData", trait, outcome, min, max))
  save(perm_res, file = sprintf("/scratch/edbeck/perm/perm_%s_%s_%s_%s.RData", trait, outcome, min, max))
}


# sca_nested <- crossing(
#   Trait = unique(p_waves$p_item),
#   Outcome = (codebook$codebook[[2]] %>% filter(category == "out"))$name
# ) %>%
#   filter(Outcome == "frstjob" & Trait == "C") %>%
#   mutate(terms = map2(Trait, Outcome, sca_term_fun),
#          formula = map(terms, sca_formula_fun),
#          ncomb = pmap(list(Trait, Outcome, formula, "p_value"), sca_setup_fun))

# args <- commandArgs(trailingOnly = TRUE)
sca_run_fun(args[,1], args[,2], args[,3], args[,4])