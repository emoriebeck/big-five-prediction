wd <- "~/selection"
lib_loc <- "~/My_R_Libraries"
library(lme4)
library(haven)
library(readxl)
# library(parallel)
library(future, lib.loc = lib_loc)
library(future.apply)
library(future.batchtools, lib.loc = lib_loc)
library(batchtools)
library(listenv, lib.loc = lib_loc)
library(plyr)
library(tidyverse)

source(sprintf("%s/mumin_fun.R", wd))
source(sprintf("%s/sca_run.R", wd))
spl_fun <- function(trait, outcome){
  load(sprintf("%s/results/sca_workspace/sca_workspace_%s_%s.RData", wd, outcome, trait))
  x <- seq(0,ncomb,1)
  seq_comb = split(x, sort(x%%250))
}

jobs_nested <- crossing(
  Trait = unique(p_waves$p_item), 
  Outcome = (codebook$codebook[[2]] %>% filter(category == "out"))$name
) %>%
  filter(Outcome == "married") %>%
  mutate(data = map2(Trait, Outcome, spl_fun)) %>%
  unnest(data) %>%
  mutate(min = map_dbl(data, min),
         max = map_dbl(data, max)) %>%
  select(trait = Trait, outcome = Outcome, min, max) %>%
  write.table(., file = sprintf("%s/scripts/sca_args_married.txt", wd), row.names = F, col.names = F)

sca_comb_fun <- function(trait, outcome){
  load(sprintf("%s/results/sca_workspace/sca_workspace_%s_%s.RData", wd, outcome, trait))
  seq_comb <- seq(0,ncomb,1)
  
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
      if(any(c("parEdu", "parOccPrstg", "grsWages") %in% new_terms)){
        if(!all(c("parEdu", "parOccPrstg", "grsWages") %in% new_terms)){
          isok <- FALSE
          res <- NA
        }
      }
    }
    if (!isok) {
      res <- NA
    } else{
      res <- 1
    }
    return(res)
  })
  res <- unlist(res)
  res <- which(!is.na(res))-1
}

sca_arrays <- tibble::tibble(
  Trait = unique(p_waves$p_item), 
  Outcome = "retired",
  data = map2(Trait, Outcome, ~sca_comb_fun(.x, .y))
  ) %>%
  unnest(data) 

save_fun <- function(trait, outcome, df){
  tibble(Trait = trait, Outcome = outcome, comb = df) %>%
  write.table(., file = sprintf("%s/scripts/%s_%s_sca_args.txt", 
        wd, trait, outcome), row.names = F, col.names = F)
}

sca_arrays %>%
  mutate(data = pmap(list(Trait, Outcome, data), save_fun))
  select(trait = Trait, outcome = Outcome, min, max) %>%
  write.table(., file = sprintf("%s/scripts/sca_args.txt", wd), row.names = F, col.names = F)

resrc <- list(nodes = 1, walltime = 24*60*60, memory = 6144)

reg = makeRegistry(NA)
reg$cluster.functions = makeClusterFunctionsTORQUE(template = "torque-lido", scheduler.latency = 1)
reg$packages = c("lme4", "haven", "readxl", "listenv", "plyr", "tidyverse")
# jobs_nested <- sca_nested %>%
jobs_nested <- crossing(
  Trait = "C", 
  Outcome = "frstjob"
  ) %>%
  
  mutate(data = map(ncomb, spl_fun)) %>%
  unnest(data) %>%
  mutate(min = map_dbl(data, min),
         max = map_dbl(data, max)) %>%
  select(trait = Trait, outcome = Outcome, min, max) %>%
  write.table(., file = sprintf("%s/scripts/sca_args.txt", wd), row.names = F, col.names = F)
  # unite(job, trait, outcome, sep = " ", remove = F) %>%
  # mutate(map2(row_number(), seq_comb, save_fun)
  # filter(row_number() %in% 1) %>%
  # group_by(job) %>%
  # nest() %>%
  # ungroup() %>%
  # mutate(ids = map(data, ~batchMap(sca_run_fun, args = ., reg = reg)),
  #        ids = map(ids, ~getJobPars(reg = reg)),
  #        sub = map(ids, ~submitJobs(ids = ., resources = resrc, reg = reg)))
  mutate(map)

## if it all goes badly wrong run this to delete and start over
removeRegistry(tmp)


# trait = "C"
# outcome = "frstjob"
# fixed = "p_value"
# f = sca_formula_fun(sca_term_fun("C", "frstjob"))