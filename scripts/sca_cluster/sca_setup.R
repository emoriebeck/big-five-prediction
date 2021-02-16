wd <- "~/selection"
lib_loc <- "~/My_R_Libraries"
library(lme4)
library(haven)
library(readxl)
# library(parallel)
library(future, lib.loc = lib_loc)
library(future.apply)
library(future.batchtools, lib.loc = lib_loc)
library(listenv, lib.loc = lib_loc)
library(plyr)
library(tidyverse)

sheets <- sprintf("%s/codebooks/master_codebook_01.24.20.xlsx", wd) %>% excel_sheets()

# function for reading in sheets
read_fun <- function(x){
  sprintf("%s/codebooks/master_codebook_01.24.20.xlsx", wd) %>% read_xlsx(., sheet = x)
}

# read in sheets and index source
codebook <- tibble(
  study = sheets,
  codebook = map(study, read_fun)
)

p_waves <- sprintf("%s/codebooks/personality_waves.xlsx", wd) %>% read_xlsx()
# used covariates for specifications 
specifications <- sprintf("%s/codebooks/specifications.xlsx", wd) %>% read_xlsx()


get_call <- function (x) 
{
  rval <- if (isS4(x)) {
    if (any(i <- (sln <- c("call", "CALL", "Call")) %in% 
            methods::slotNames(x))) 
      slot(x, sln[i][1L])
    else if (!is.null(attr(x, "call"))) 
      attr(x, "call")
    else NULL
  }
  else {
    if (!is.atomic(x) && (i <- match("call", names(x), nomatch = 0L)) != 
        0L) {
      x[[i]]
    }
    else if (!is.null(attr(x, "call"))) {
      attr(x, "call")
    }
    else NULL
  }
  if (is.null(rval)) 
    stats::getCall(x)
  else rval
}


fixCoefNames <- function (x, peel = TRUE) 
{
  if (!length(x)) 
    return(x)
  ox <- x
  ia <- grep(":", x, fixed = TRUE)
  if (!length(ia)) 
    return(structure(x, order = rep.int(1L, length(x))))
  x <- ret <- x[ia]
  if (peel) {
    if (all(substr(x, 1L, pos <- regexpr("_", x, fixed = TRUE)) %in% 
            c("count_", "zero_"))) {
      ret <- substr(ret, pos + 1L, 256L)
      k <- TRUE
      suffix <- ""
    }
    else {
      k <- grepl("^\\w+\\(.+\\)$", x, perl = TRUE)
      fname <- substring(x[k], 1L, attr(regexpr("^\\w+(?=\\()", 
                                                x[k], perl = TRUE), "match.length"))
      k[k] <- !vapply(fname, exists, FALSE, mode = "function", 
                      envir = .GlobalEnv)
      if (any(k)) {
        pos <- vapply(x[k], function(z) {
          parens <- lapply(lapply(c("(", ")"), function(s) gregexpr(s, 
                                                                    z, fixed = TRUE)[[1L]]), function(y) y[y > 
                                                                                                             0L])
          parseq <- unlist(parens, use.names = FALSE)
          p <- cumsum(rep(c(1L, -1L), sapply(parens, 
                                             length))[order(parseq)])
          if (any(p[-length(p)] == 0L)) 
            -1L
          else parseq[1L]
        }, 1L, USE.NAMES = FALSE)
        k[k] <- pos != -1L
        pos <- pos[pos != -1]
        if (any(k)) 
          ret[k] <- substring(x[k], pos + 1L, nchar(x[k]) - 
                                1L)
      }
      suffix <- ")"
    }
  }
  else k <- FALSE
  spl <- expr.split(ret, ":", prepare = function(x) gsub("((?<=:):|:(?=:))", 
                                                         "_", x, perl = TRUE))
  ret <- vapply(lapply(spl, base::sort), paste0, "", collapse = ":")
  if (peel && any(k)) 
    ret[k] <- paste0(substring(x[k], 1L, pos), ret[k], suffix)
  ox[ia] <- ret
  ord <- rep.int(1, length(ox))
  ord[ia] <- sapply(spl, length)
  structure(ox, order = ord)
}

std.coef <- function (x, partial.sd, ...) 
{
  coefmat <- coefTable(x, ...)
  mm <- model.matrix(x)
  if (partial.sd) {
    bx <- .partialsd(coefmat[, 1L], apply(mm, 2L, sd), .vif(x), 
                     nobs(x), sum(attr(mm, "assign") != 0))
  }
  else {
    response.sd <- sd(model.response(model.frame(x)))
    bx <- apply(mm, 2L, sd)/response.sd
  }
  coefmat[, 1L:2L] <- coefmat[, 1L:2L] * bx
  colnames(coefmat)[1L:2L] <- c("Estimate*", "Std. Error*")
  return(coefmat)
}

matchCoef <- function (m1, m2, all.terms = getAllTerms.fixest(m2, terms), 
                       beta = 0L, terms1 = names(coef(m1)), coef1 = NULL, 
                       allCoef = FALSE, ...) 
{
  if (is.null(coef1)) {
    ct <- if (beta != 0L) 
      std.coef(m1, beta == 2L, ...)
    else coefTable.fixest(m1)
    coef1 <- ct[, 1L]
    names(coef1) <- rownames(ct)
  }
  else if (allCoef) 
    stop("'coef1' is given and 'allCoef' is not FALSE")
  if (any((terms1 %in% all.terms) == FALSE)) 
    stop("'m1' is not nested within 'm2'")
  row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)
  fxdCoefNames <- fixCoefNames(names(coef1))
  row[terms1] <- NaN
  pos <- match(terms1, fxdCoefNames, nomatch = 0L)
  row[fxdCoefNames[pos]] <- coef1[pos]
  if (allCoef) {
    i <- match(names(coef1), rownames(ct))
    j <- !is.na(i)
    rownames(ct)[i[j]] <- fxdCoefNames[j]
    attr(row, "coefTable") <- ct
  }
  row
}


formula_margin_check <- function (j, m) {
  stopifnot(is.logical(j))
  !any(m[!j, j], na.rm = TRUE)
}

`getAllTerms.default` <-
  #function(x, ...) getAllTerms.formula(as.formula(formula(x)), ...)
  function(x, ...) getAllTerms.terms(terms(as.formula(formula(x))), ...)

`getAllTerms.gam` <-
  function(x, intercept = FALSE, offset = TRUE, ...)
    getAllTerms.terms(terms(formula(x), ...), intercept = intercept, offset = offset)

`getAllTerms.lm` <-
  function(x, intercept = FALSE, offset = TRUE, ...)
    getAllTerms.terms(terms(x, ...), intercept = intercept, offset = offset)

`getAllTerms.fixest` <-
  function(m, terms, x = terms(as.formula(formula(m))), intercept = FALSE, offset = TRUE){
    interceptLabel <- "(Intercept)"
    variables <- terms
    # variables <- names(coef(m))
    ans <- terms
    # ans <- names(coef(m))
    nvars <- length(variables)
    
    factors <- matrix(rep(0,nvars^2), nrow = nvars)
    colnames(factors) <- ans; rownames(factors) <- ans
    diag(factors) <- 1
    factors <- factors[order(rownames(factors)), , drop = FALSE]
    v <- rownames(factors)
    ans <- apply(factors != 0L, 2L, function(x) paste0(v[x], collapse = ":"))
    
    deps <- if (length(ans) > 0L) termdepmat(reformulate(ans))
    
    dimnames(deps) <- list(ans, ans)
    diag(deps) <- NA
    
    if(intercept && attr(x, "intercept")) {
      ans <- c(interceptLabel, ans)
      ord <- c(1L, ord + 1L)
    }
    
    attr(ans, "intercept") <- attr(x, "intercept")
    attr(ans, "interceptLabel") <- interceptLabel
    
    attr(ans, "deps") <- deps
    ans
  }

`getAllTerms.terms` <-
  function(x, intercept = FALSE, offset = TRUE, ...) {
    #function(x, offset = TRUE, intercept = FALSE, ...) { # XXX!
    
    interceptLabel <- "(Intercept)"
    variables <- attr(x, "variables")[-1L]
    # variables <- names(coef(x2))
    
    if (!is.null(attr(x, "offset"))){
      offs <- sapply(variables[attr(x, "offset")], deparse)
    } else offs <- NULL
    
    ans <- attr(x, "term.labels")
    # ans <- names(coef(x2))
    
    # Get term names, with higher order term components arranged alphabetically
    if (length(ans) > 0L) {
      factors <- attr(x, "factors")
      factors <- factors[order(rownames(factors)), , drop = FALSE]
      v <- rownames(factors)
      ans <- apply(factors != 0L, 2L, function(x) paste0(v[x], collapse = ":"))
    }
    
    # Leave out random terms (lmer type)
    #ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
    
    .is.re <- function(x) {
      n <- length(x)
      if(n == 3L && x[[1L]] == "|") return(1L)
      if(n == 2 && is.call(x[[2L]]) && x[[2L]][[1L]] == "|") return(2L)
      return(0L)
    }
    
    reType <- vapply(variables, .is.re, 0L)
    # 1 -> (terms|group), 2 -> struc(terms|group)
    
    ran <- as.character(variables[reType != 0L])
    ifx <- !(ans %in% ran)
    
    ans <- ans[ifx] # ifx - indexes of fixed terms
    #retUnsorted <- ans
    
    # finally, sort by term order and then alphabetically
    #ans <- unname(ans[order(attr(x, "order")[ifx], ans)])
    ord <- order(attr(x, "order")[ifx], gsub("I\\((.*)\\)", "\\1", ans))
    ans <- unname(ans[ord])
    
    deps <- if (length(ans) > 0L) termdepmat(reformulate(ans)) else
      matrix(FALSE, 0L, 0L)
    
    dimnames(deps) <- list(ans, ans)
    diag(deps) <- NA
    
    if(intercept && attr(x, "intercept")) {
      ans <- c(interceptLabel, ans)
      ord <- c(1L, ord + 1L)
    }
    
    if (!is.null(offs[1L])) {
      if (offset) {
        ans <- c(ans, offs)
        ord <- c(ord, length(ord) + 1L)
      }
      attr(ans, "offset") <- offs
    }
    attr(ans, "intercept") <- attr(x, "intercept")
    attr(ans, "interceptLabel") <- interceptLabel
    
    if (length(ran) > 0L) {
      attr(ans, "random.terms") <- ran
      i <- reType[reType != 0L] == 1L
      ran1 <- ran
      ran1[i] <- paste0("(", ran1[i], ")")
      f.random <- reformulate(c(".", ran1), response = ".")
      environment(f.random) <- environment(x)
      attr(ans, "random") <- f.random
    }
    
    response <- attr(x, "response")
    response <- if(response == 0L) NULL else variables[[response]]
    attr(ans, "response") <- response
    attr(ans, "order") <- order(ord)
    attr(ans, "deps") <- deps
    ans
  }

`getAllTerms.formula` <-
  function(x, ...) getAllTerms.terms(terms.formula(x), ...)

`getAllTerms.lme` <-
  function(x, ...) {
    termsobj <- if(inherits(x, "glmmPQL"))
      terms(formula(x), data = x$data) else
        terms(x)
    
    ret <- getAllTerms.terms(termsobj, ...)
    attr(ret, "random") <- . ~ .
    
    # Code from nlme:::print.reStruct, modified slightly
    reStruct <- x$modelStruct$reStruct
    nobj <- length(reStruct)
    if (is.null(namx <- names(reStruct)))
      names(reStruct) <- nobj:1L
    aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
    aux[lower.tri(aux)] <- ""
    reStruct[] <- rev(reStruct)
    aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
    aux[lower.tri(aux)] <- ""
    attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
                                              "[[", 2L), "|",
                                       rev(apply(aux, 1L, function(z) paste(z[z != ""], collapse = " %in% "))))
    
    return(ret)
  }

# Apparently there is no (explicit) intercept in coxph, but 'terms' gives
# attr(,"intercept") == 1.
`getAllTerms.coxph` <- function (x, ...) {
  ret <- getAllTerms.default(x, ...)
  attr(ret, "intercept") <- 0L
  attr(ret, "interceptLabel") <- NULL
  return(ret)
}

`getAllTerms.glmmML` <- function (x, ...) {
  ret <- getAllTerms.terms(terms(x), ...)
  attr(ret, "random.terms") <-  paste("1 |",  x$call$cluster)
  return(ret)
}

#`getAllTerms.hurdle` <- function(x, intercept = FALSE, ...) {
#	f <- as.formula(formula(x))
#	# to deal with a dot in formula (other classes seem to expand it)
#	if("." %in% all.vars(f))
#		getAllTerms.terms(terms.formula(f, data = eval(x$call$data, envir = environment(f)))
#			
#			, intercept = intercept)
#	else getAllTerms.formula(f, intercept = intercept)
#}

split_formula_by_bar <- function(f) {
  n <- length(f)
  ans <- if(length(f[[n]]) != 1L && f[[n]][[1L]] == "|") {
    f1 <- vector("list", 2L)
    for(i in 1L:2L) {
      f1[[i]] <- f
      f1[[i]][[n]] <- f[[n]][[i+ 1]]
    }
    f1
  } else list(f)
  ans
}

`getAllTerms.hurdle` <- 
  `getAllTerms.zeroinfl` <-
  function(x, intercept = FALSE, ...) {
    
    formList <- split_formula_by_bar(formula(x))
    formList <- lapply(lapply(formList, terms.formula, data = eval(x$call$data)),
                       formula)
    z <- lapply(formList, getAllTerms, intercept = TRUE)
    
    if(oneform <- length(formList) == 1L) z <- c(z, z)
    
    deps <- termdepmat_combine(lapply(z, attr, "deps"))
    
    ord <- unlist(lapply(z, attr, "order"))
    n <- sapply(z, length)
    if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
    
    zz <- unlist(z)
    interceptIdx <- zz == "(Intercept)"
    offsetIdx <- match(zz, unique(unlist(lapply(z, attr,"offset"))), nomatch = 0) != 0
    termIdx <- !(offsetIdx | interceptIdx)
    
    zz <- paste0(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
                 "_", zz)
    
    dimnames(deps) <- list(zz[termIdx], zz[termIdx])
    
    if(oneform) { # dependency of count_X and zero_X
      k <- length(zz[termIdx]) / 2
      deps[c(seq(k + 1L, by = 2L * k + 1L, length.out = k),
             seq((2L * k * k) + 1L, by = 2L * k + 1L, length.out = k))] <- TRUE
    }
    
    ret <- if(!intercept) zz[!interceptIdx] else zz
    if(any(offsetIdx)) attr(ret, "offset") <- zz[offsetIdx]
    attr(ret, "intercept") <- pmin(which(interceptIdx), 1)
    attr(ret, "interceptLabel") <- zz[interceptIdx]
    attr(ret, "response") <- attr(z[[1L]], "response")
    attr(ret, "order") <- if(!intercept) order(ord[!interceptIdx]) else ord
    attr(ret, "deps") <- deps
    ret
  }

## TODO: test with offsets
`getAllTerms.betareg` <-
  function(x, intercept = FALSE, ...) {
    formList <- split_formula_by_bar(formula(x))
    formList <- lapply(lapply(formList, terms.formula, data = model.frame(x)),
                       formula)
    oneform <- length(formList) == 1L
    z <- lapply(formList, getAllTerms, intercept = TRUE)
    
    deps <- termdepmat_combine(lapply(z, attr, "deps"))
    
    ord <- unlist(lapply(z, attr, "order"))
    n <- sapply(z, length)
    if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
    zz <- unlist(z)
    interceptIdx <- zz == "(Intercept)"
    offsetIdx <- match(zz, unique(unlist(lapply(z, attr,"offset"))), nomatch = 0) != 0
    termIdx <- !(offsetIdx | interceptIdx)
    
    if(!oneform && n[2L] != 0L) {
      i.phi <- -seq.int(n[1L])
      zz[i.phi] <- paste("(phi)", zz[i.phi], sep = "_")
    }
    dimnames(deps) <- list(zz[termIdx], zz[termIdx])
    
    ret <- if(!intercept) zz[!interceptIdx] else zz
    if(any(offsetIdx)) attr(ret, "offset") <- zz[offsetIdx]
    attr(ret, "intercept") <- pmin(which(interceptIdx), 1)
    attr(ret, "interceptLabel") <- zz[interceptIdx]
    attr(ret, "response") <- attr(z[[1L]], "response")
    attr(ret, "order") <- if(!intercept) order(ord[!interceptIdx]) else ord
    attr(ret, "deps") <- deps
    ret
  }




`getAllTerms.glimML` <- function(x, intercept = FALSE, ...) {
  ret <- getAllTerms.default(x, intercept = intercept, ...)
  ttran <- terms.formula(x@random)
  ran <- attr(ttran, "term.labels")
  if(length(ran)) attr(ret, "random.terms") <- paste("1 |", ran)
  ret
}

`getAllTerms.coxme` <-
  function(x, ...)  {
    ret <- getAllTerms.terms(terms(x))
    random <- x$formulaList$random
    attr(ret, "random.terms") <- as.character(random)
    f <- as.name(".")
    for(f1 in random) f <- call("+", f, f1)
    attr(ret, "random") <- call("~", as.name("."), f)
    attr(ret, "intercept") <- 0L
    attr(ret, "interceptLabel") <- NULL
    ret
  }

`getAllTerms.MCMCglmm` <- 
  function (x, ...) {
    res <- getAllTerms.default(x, ...) 
    attr(res, "random") <- .formulaEnv(.~., environment(formula(x)))
    attr(res, "random.terms") <- asChar(x$Random$formula)[1L]
    res
  }

`getAllTerms.gamm` <-
  function (x, ...) getAllTerms(x$gam, ...)


`getAllTerms.mark` <- 
  function (x, intercept = FALSE, ...) {
    
    f <- formula(x, expand = FALSE)[[2L]]
    formlist <- list()
    while(length(f) == 3L && f[[1L]] == "+") {
      formlist <- c(f[[3L]], formlist)
      f <- f[[2L]]
    }
    formlist <- append(f, formlist)
    
    wrapfunc <- function(x, func) if(length(x) == 0L) x else paste0(func, "(", x, ")")
    
    alltermlist <- lapply(formlist, function(x, intercept) {
      func <- asChar(x[[1L]])
      at <- getAllTerms(terms(eval(call("~", x[[2L]]))), intercept = intercept)
      at[] <- wrapfunc(at, func)
      dn <- wrapfunc(rownames(attr(at, "deps")), func)
      attr(at, "interceptLabel") <- wrapfunc(attr(at, "interceptLabel"), func)
      dimnames(attr(at, "deps")) <- list(dn, dn)
      at
    }, intercept)
    
    retval <- unlist(alltermlist, recursive = TRUE)
    for(a in c("intercept", "interceptLabel")) {
      attr(retval, a) <-	unlist(sapply(alltermlist, attr, a))
    }
    attr(retval, "order") <- order(rep(seq_along(alltermlist), vapply(alltermlist, length, 1L)),
                                   unlist(lapply(alltermlist, attr, "order")))
    attr(retval, "deps") <- termdepmat_combine(lapply(alltermlist, attr, "deps"))
    retval
  }

`getAllTerms.asreml`  <-
  function(x, intercept = FALSE, ...)
    getAllTerms.terms(terms(formula(x)), intercept = intercept, ...)

`getAllTerms.cpglmm` <-
  function (x, intercept = FALSE, ...) 
    getAllTerms(x@formula, intercept = intercept, ...)

`getAllTerms` <-
  function(x, ...)
    UseMethod("getAllTerms")

# TODO: return object of class 'allTerms'
print.allTerms <-
  function(x, ...) {
    cat("Model terms: \n")
    if(!length(x)) {
      cat("<None> \n")
    } else {
      print.default(as.vector(x), quote = TRUE)
    }
    ints <- attr(x, "interceptLabel")
    if(!is.null(ints)) {
      cat(ngettext(n = length(ints), "Intercept:", "Intercepts:"), "\n")
      print.default(ints,quote = TRUE)
    }
  }

termlist <- function(x) {
  is.plus <- function(x) is.call(x) && x[[1L]] == "+"
  ## parses interaction expression into list: a:b:c --> list(a,b,c)
  intr <- function(x) {
    # is it an expression for interaction? (e.g. a:b:c)
    is.intr <- function(x) is.call(x) && x[[1L]] == ":"
    if(is.intr(x)) {
      res <- list()
      repeat {
        res <- c(x[[3L]], res)
        x <- x[[2L]]
        if(!is.intr(x)) break
      }
      list(c(x, res))
    } else x
  }
  if(x[[1L]] == "~") x <- x[[length(x)]]
  res <- list()
  while(is.plus(x)) {
    res <- c(intr(x[[3L]]), res)
    x <- x[[2L]]
  }
  res <- c(intr(x), res)
  res	
}

# calculates all lower order term names:
# expandintr(1:3) --> c("1", "2", "1:2", "3", "1:3", "2:3", "1:2:3")
expandintr <- function(x) {
  asstr <- function(x) asChar(x, backtick = TRUE)
  if(!is.language(x)) {
    a <- sapply(x, asstr)
    k <- length(a)
    vapply(seq.int(2L^k - 1L), function(y) paste0(a[as.logical(intToBits(y)[1L:k])],
                                                  collapse = ":"), "")
  } else asstr(x)
}

# evaluate 'expr' in 'env' after adding variables passed as '...'
evalExprInEnv <- function(expr, env, enclos, ...) {
  list2env(list(...), env)
  eval(expr, envir = env, enclos = enclos)
}

# change `names[]` for varName[1], varName[2], ... in expression
# Not using `substitute` anymore to omit function calls.
# Ignore also expressions within I(), elements extracted with $ or @
`.subst.names.for.items` <-
  function(expr, names, varName, n = length(names), fun = "[") {
    exprApply(expr, names, symbols = TRUE,
              function(x, v, fun, varName, parent) {
                if(is.call(parent) && any(parent[[1L]] == c("I", "$", "@")))
                  return(x)
                if(length(x) == 1L)
                  return(call(fun, varName, match(asChar(x), v)))
                x
              }, v = names, fun = fun, varName = as.name(varName))
  }

# like substitute, but does evaluate 'expr'.
subst <-
  function(expr, envir = NULL, ...) {
    eval.parent(call("substitute", expr, c(envir, list(...))))
  }

asChar <- function(x, control = NULL, nlines = 1L, ...)
  if(is.character(x)) x[1L:nlines] else
    deparse(x, control = control, nlines = nlines, ...)

## .sub_* functions used with '.exprapply' as 'func'
.subst.term <- function(x) {
  if(length(x) < 2L) cry(x, "'Term' needs one argument")
  as.name(asChar(x[[2L]]))
}

.subst.with <- function (x, fac, allTerms, vName, envir = parent.frame()) {
  if (length(x) > 4L) cry(x, "too many arguments [%d]", length(x) - 1L)
  if (length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
    fun <- "all"
    sx <- asChar(x[[2L]][[2L]], backtick = FALSE)
  } else {
    fun <- "any"
    sx <- asChar(x[[2L]], backtick = FALSE)
  }
  dn <- dimnames(fac)
  if (!(sx %in% dn[[2L]])) cry(x, "unknown variable name '%s'", sx)
  xorder <- if(length(x) >= 3L) as.integer(eval(x[[3L]], envir))
  else unique(rowSums(fac))
  i <- which(fac[, sx])
  j <- which(is.element(rowSums(fac[i, , drop = FALSE]), xorder))
  if(length(j) == 0L) cry(x, "no terms match the criteria")    
  as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), match(dn[[1L]][i[j]], allTerms))))))
}

.subst.vars.for.args <- function(e) {
  for(i in 2L:length(e))
    if(!is.name(e[[i]]))
      e[[i]] <- as.name(asChar(e[[i]]))
    e
}

.subst.has <- function(e) {
  n <- length(e)
  for(i in seq.int(2L, n)) {
    ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
      call("is.na", e[[i]][[2L]]) else
        call("!", call("is.na", if(is.name(e[[i]])) e[[i]] else
          as.name(asChar(e[[i]]))))
    res <- if(i == 2L) ex else call("&", res, ex)
  }
  call("(", res)
}

.subst.has.dc <- function(e) {
  for(i in 2L:length(e)) e[[i]] <- call("has", e[[i]])
  e
}

.subst.v <- function(x, cVar, fn) {
  if(length(x) > 2L) cry(x, "discarding extra arguments", warn = TRUE)
  i <- which(fn == x[[2L]])[1L]
  if(is.na(i)) cry(x, "'%s' is not a valid name of 'varying' element",
                   as.character(x[[2L]]), warn = TRUE)
  call("[[", cVar, i)
}

# substitute function calls in 'e'. 'func' must take care of the substitution job.
`exprapply0` <- function(e, name, func, ...)
  exprApply(e, name, func, ..., symbols = FALSE)

`exprApply` <-
  function (expr, what, FUN, ..., symbols = FALSE) {
    FUN <- match.fun(FUN)
    funcl <- as.call(c(as.name("FUN"), as.name("expr"), list(...)))
    if(all(names(formals(FUN)) != "parent"))
      formals(FUN)[["parent"]] <- NA
    .exprapply(expr, what, FUN, ..., symbols = symbols)
  }

`.exprapply` <-
  function (expr, what, FUN, ..., symbols = FALSE, parent = NULL) {
    self <- sys.function()
    if((ispairlist <- is.pairlist(expr)) || is.expression(expr)) {
      for (i in seq_along(expr))	expr[i] <-
          list(self(expr[[i]], what, FUN, ..., symbols = symbols, parent = expr))
      return(if(ispairlist) as.pairlist(expr) else expr)
    }
    n <- length(expr)
    if (n == 0L)
      return(expr) else
        if (n == 1L) {
          if (!is.call(expr)) {
            if (symbols && (anyNA(what) || any(expr == what)))
              expr <- FUN(expr, ..., parent = parent)
            return(expr)
          }
        } else {
          if(expr[[1L]] == "function") {
            if(n == 4L) {
              n <- 3L
              expr[[4L]] <- NULL ## remove srcref
            }
          }
          for (i in seq.int(2L, n)) {
            y <- self(expr[[i]], what, FUN, ..., symbols = symbols, parent = expr)
            if(!missing(y)) expr[i] <- list(y)
          }
        }
    if (anyNA(what) || (length(expr[[1L]]) == 1L && any(expr[[1L]] == what)))
      expr <- FUN(expr, ..., parent = parent)
    return(expr)
  }

`tTable` <-
  function (model, ...) 	{
    .Deprecated("coefTable")
    coefTable(model, ...)
  }

`coefTable` <-
  function (model, ...) UseMethod("coefTable")

.makeCoefTable <- 
  function(x, se, df = NA_real_, coefNames = names(x)) {
    if(n <- length(x)) {
      xdefined <- !is.na(x)
      ndef <- sum(xdefined)
      if(ndef < n) {
        if(length(se) == ndef) {
          y <- rep(NA_real_, n); y[xdefined] <- se; se <- y
        }
        if(length(df) == ndef) {
          y <- rep(NA_real_, n); y[xdefined] <- df; df <- y
        }
      }
    }
    if(n && n != length(se)) stop("length(x) is not equal to length(se)")
    ret <- matrix(NA_real_, ncol = 3L, nrow = length(x),
                  dimnames = list(coefNames, c("Estimate", "Std. Error", "df")))
    if(n) ret[, ] <- cbind(x, se, rep(if(is.null(df)) NA_real_ else df,
                                      length.out = n), deparse.level = 0L)
    class(ret) <- c("coefTable", "matrix")
    ret
  }

`coefTable.default` <-
  function(model, ...) {
    dfs <- tryCatch(df.residual(model), error = function(e) NA_real_)
    cf <- summary(model, ...)$coefficients
    .makeCoefTable(cf[, 1L], cf[, 2L], dfs, coefNames = rownames(cf))
  }

`coefTable.fixest` <-
  function(model, ...) {
    dfs <- tryCatch(df.residual(model), error = function(e) NA_real_)
    cf <- model$coeftable
    .makeCoefTable(cf[, 1L], cf[, 2L], dfs, coefNames = rownames(cf))
  }

`coefTable.lm` <-
  function(model, ...)
    .makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))), model$df.residual)


`coefTable.survreg` <- 
  function(model, ...) {
    .makeCoefTable(
      coeffs(model), 
      sqrt(diag(vcov(model, ...))), 
      NA,
    )
  }

`coefTable.coxph` <-
  function(model, ...) {
    .makeCoefTable(coef(model), if(all(is.na(model$var))) 
      rep(NA_real_, length(coef(model))) else sqrt(diag(model$var)), 
      model$df.residual)
  }

`coefTable.glmmML` <- function(model, ...)
  .makeCoefTable(model$coefficients, model$coef.sd)

`coefTable.gls` <-
  function (model, ...)
    .makeCoefTable(coef(model), sqrt(diag(as.matrix(model$varBeta))),
                   model$dims$N - model$dims$p)

`coefTable.lme` <-
  function(model, adjustSigma = TRUE, ...) {
    se <- sqrt(diag(as.matrix(model$varFix)))
    if (adjustSigma && model$method == "ML")
      se <- se * sqrt(model$dims$N / (model$dims$N - length(se)))
    .makeCoefTable(nlme::fixef(model), se, model$fixDF[["X"]])
  }

`coefTable.multinom` <- 
  function (model, ...) {
    .makeCoefTable(coeffs(model), sqrt(diag(vcov(model, ...))))
  }

`coefTable.sarlm` <-
  `coefTable.spautolm` <-
  function(model, ...) {
    x <- coef(model)
    .makeCoefTable(x, sqrt(diag(summary(model, ...)$resvar))[names(x)])
  }

`coefTable.coxme` <-
  `coefTable.lmekin` <-
  function(model, ...)  {
    # code from coxme:::print.coxme
    beta <- model$coefficients # for class coxme:
    if(is.list(beta) && !is.null(beta$fixed))
      beta <- beta$fixed # for class lmekin and older coxme
    nvar <- length(beta)
    if(nvar) {
      nfrail <- nrow(model$var) - nvar
      se <- sqrt(get("diag", getNamespace("Matrix"))(model$var)[nfrail + 1L:nvar])
    } else se <- NULL
    .makeCoefTable(beta, se)
  }

`coefTable.rq` <- 
  function(model, ...)
    .makeCoefTable(model$coefficients, rep(NA_real_, length(model$coefficients)))

`coefTable.zeroinfl` <-
  function(model, ...)
    .makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.hurdle` <- 
  function(model, ...) {
    cts <- summary(model)$coefficients
    ct <- do.call("rbind", unname(cts))
    cfnames <- paste0(rep(names(cts), vapply(cts, nrow, 1L)), "_", rownames(ct))
    .makeCoefTable(ct[, 1L], ct[, 2L], coefNames = cfnames)
    #.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))
  }

`coefTable.aodql` <-
  `coefTable.betareg` <- 
  `coefTable.glimML` <-
  `coefTable.unmarkedFit` <- 
  function(model, ...)
    .makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.gee` <-
  `coefTable.geeglm` <-
  function(model, ..., type = c("naive", "robust")) {
    cf <- summary(model, ...)$coefficients
    j <- if(match.arg(type) == "naive") 2L else 4L
    .makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
  }

`coefTable.geem` <-
  function(model, ..., type = c("naive", "robust")) {
    smr <- summary(model)
    .makeCoefTable(smr$beta, smr[[if(match.arg(type) == "naive")
      "se.model" else "se.robust"]],
      coefNames = smr$coefnames)
  }

`coefTable.geese` <-
  function(model, ..., type = c("naive", "robust")) {
    cf <- summary(model, ...)$mean
    type <- match.arg(type)
    j <- if(type == "naive") 2L else 4L
    .makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
  }

`coefTable.yagsResult` <-
  function(model, ..., type = c("naive", "robust")) {
    type <- match.arg(type)
    vcv <- slot(model, if(type == "naive") "naive.parmvar" else "robust.parmvar")
    .makeCoefTable(model@coefficients, sqrt(diag(vcv)), coefNames = model@varnames)
  }

`coefTable.splm` <- 
  function (model, ...) {
    cf <- sapply(c("coefficients", "arcoef", "errcomp"), function(i)
      if(is.matrix(model[[i]])) model[[i]][, 1L] else model[[i]],
      simplify = FALSE)
    
    ncf <- sapply(cf, length)
    vcovlab <- c(coefficients = "vcov", arcoef = "vcov.arcoef", errcomp = "vcov.errcomp")
    se <- sqrt(unlist(lapply(names(vcovlab), function(i) {
      vcv2 <- diag(model[[vcovlab[i]]])
      c(vcv2, rep(NA_real_, ncf[[i]] - length(vcv2)))
    })))
    
    .makeCoefTable(unlist(cf, use.names = FALSE), se,
                   coefNames = unlist(lapply(cf, names), use.names = FALSE))
  }

`coefTable.MCMCglmm` <-
  function (model, ...) {
    cf <- coeffs(model)
    .makeCoefTable(cf, se = rep(NA_real_, length.out = length(cf)))
  }

`coefTable.gamm` <-
  function (model, ...) coefTable.lm(model$gam, ...)

`coefTable.mark` <- 
  function (model, orig.names = FALSE, ...) {
    dfs <- model$results[['n']] - model$results[['npar']]
    beta <- model$results[['beta']]
    .makeCoefTable(beta[, 1L], beta[, 2L], dfs,
                   coefNames = if(orig.names) rownames(beta) else
                     gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)", rownames(beta), perl = TRUE))
  }

`coefTable.logistf` <-
  function (model, ...)
    .makeCoefTable(model$coefficients, sqrt(diag(model$var)))

`coefTable.aodml` <-
  function (model, ...) {
    .makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))
    #.makeCoefTable(coeffs(model), sqrt(diag(model$varparam)))
  }

## XXX: fixed effects coefficients only
`coefTable.asreml` <- 
  function (model, ...)  {
    .makeCoefTable(
      x = model$coefficients$fixed, 
      se = sqrt(model$vcoeff$fixed * model$sigma2) ## ?
    ) 
  }

`coefTable.cplm` <-
  function (model, ...) 
    .makeCoefTable(coef(model), sqrt(diag(vcov(model))),
                   model@df.residual)

`coefTable.cpglmm` <-
  function (model, ...) 
    .makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))


`coefTable.maxlikeFit` <-
  function (model, ...)
    .makeCoefTable(model$Est[, 1L], model$Est[, 2L])


coefTable.bic.glm <-
  function (model, ...) {
    .makeCoefTable(model$condpostmean, model$condpostsd, NA_integer_,
                   dimnames(model$mle)[[2L]])
  }

# coefTable methods:

`print.coefTable` <-
  function (x, ...)
    stats::printCoefmat(x[, if(all(is.na(x[, 3L]))) -3L else TRUE, drop = FALSE],
                        has.Pvalue = FALSE)

summary.coefTable <-
  function (object, ...) {
    tvalue <- object[, 1L] / object[, 2L]
    if (all(is.na(object[, 3L]))) {
      pvalue <- 2 * pnorm(-abs(tvalue))
      rval <- cbind(object, tvalue, pvalue)
      cn <- c("z value", "Pr(>|z|)")
    } else if (any(is.finite(tvalue))) {
      pvalue <- 2 * pt(-abs(tvalue), object[, 3L])
      cn <- c("t value", "Pr(>|t|)")
    } else {
      pvalue <- tvalue <- NaN
      cn <- c("t value", "Pr(>|t|)")
    }
    rval <- cbind(object, tvalue, pvalue)
    colnames(rval)[4L:5L] <- cn
    class(rval) <- c("summary.coefTable", class(object))
    rval
  }

`print.summary.coefTable` <-
  function (x, signif.stars = getOption("show.signif.stars"), ...) {
    j <- if(all(is.na(x[, 3L]))) -3L else TRUE
    stats::printCoefmat(x[, j],
                        has.Pvalue = any(is.finite(x[, 5L])),
                        signif.stars = signif.stars,
                        ...)
  }



plot.mcoefTable <-
  function (x, y, labAsExpr = FALSE, n = 101, w = 5, ...) {
    lab_as_expr <- function(x) {
      x <- gsub(":", "%*%", x, perl = TRUE)
      x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
      parse(text = x)
    }
    
    xd <- function(z, n, w, x = NULL) {
      rval <- matrix(NA_real_, ncol = 3L, nrow = n)
      rval[, 1L] <- if(is.null(x)) seq(z[1L] - (w * z[2L]), z[1L] + (w * z[2L]),
                                       length.out = n) else x
      if(!is.na(z[3L]))
        rval[, 2L] <- dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
      rval[, 3L] <- dnorm(rval[, 1L], z[1L], z[2L])
      rval
    }
    
    m <- nrow(x)
    lab <- if(labAsExpr) lab_as_expr(rownames(x)) else rownames(x)
    nmodels <- dim(x)[3L]
    
    col <- 1L:nmodels
    
    par(mfrow = n2mfrow(m), ...)
    for(i in 1L:m) {
      #cat("--", dimnames(x)[[1]][i], "--\n")
      
      
      mat <- matrix(0, n, 2L * nmodels + 1L)
      for(k in 1L:nmodels) {
        j <- seq.int(length.out = 2 + (k == 1), from = 1 + (k - 1)* 2 + (k != 1))
        xlim <- c(min(x[i, 1L, ]) - w * min(x[i, 2L, ]), max(x[i, 1L,]) + w * max(x[i,2L,]))
        vx <- seq(xlim[1L], xlim[2L], length.out = n)
        v <- xd(x[i, , k], n = n, w = w, x = vx)
        mat[, j] <- v[, (2 - (k == 1)):3]
      }
      
      plot.new()
      plot.window(xlim = range(mat[, 1L]), ylim = range(mat[,-1L], na.rm = TRUE))
      j <- seq.int(2, length.out = nmodels, by = 2)
      if(any(!is.na(mat[, j])))
        matplot(mat[, 1], mat[, j], type = "l", lty = 2, add = TRUE, col = col)
      matplot(mat[, 1L], mat[, j + 1L], type = "l", lty = 1, add = TRUE, col = col)
      abline(v = x[i, 1L, ], lty = 1L, col = col)
      abline(v = 0, lty =3L, col = 8)
      axis(1L)
      axis(2L)
      box()
      
      title(lab[[i]])
      
    }
    invisible()
  }

plot.coefTable <-
  function (x, y, labAsExpr = FALSE, n = 101, w = 5,...) {
    lab_as_expr <- function(x) {
      x <- gsub(":", "%*%", x, perl = TRUE)
      x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
      parse(text = x)
    }
    xd <- function(z, n, w) {
      rval <- matrix(NA_real_, ncol = 3L, nrow = n)
      rval[, 1L] <- seq(z[1L] - w * z[2L], z[1L] + w * z[2L],
                        length.out = n)
      if(!anyNA(z[3L]))
        rval[, 2L] <- dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
      rval[, 3L] <- dnorm(rval[, 1L], z[1L], z[2L])
      rval
    }
    
    m <- nrow(x)
    lab <- if(labAsExpr) lab_as_expr(rownames(x)) else rownames(x)
    i <- 1
    nmodels <- dim(x)[3L]
    
    par(mfrow = n2mfrow(m))
    for(i in 1L:m) {	
      v <- xd(x[i, ], n = n, w = w)
      plot.new()
      plot.window(xlim = range(v[, 1L]), ylim = range(v[,-1L], na.rm = TRUE))
      lines(v[, 1L], v[, 2L])
      lines(v[, 1L], v[, 3L], col = "red")
      abline(v = c(x[i, 1L], 0), lty = c(1L, 3L))
      axis(1L)
      axis(2L)
      box()
      title(subst(expression(A == B %+-% C), A = lab[[i]], B = round(x[i, 1L], 1L), C = round(x[i, 2L], 2L)))
    }
    invisible()
  }


# given a formula, 'term dependency matrix', i.e. dependency of higer
# order terms on other lower order terms
termdepmat <- function(f) {
  trm <- terms(f, simplify = TRUE)
  tl <- termlist(trm)
  v <- attr(trm, "term.labels")
  n <- length(v)
  mat <- matrix(FALSE, n, n, dimnames = list(v, v))
  for(i in seq.int(n)) mat[match(expandintr(tl[[i]]), v), i] <- TRUE
  mat
}

# alternative to 'termdepmat', gives matrix dimension names as numbers
# so a,b,a:b  become 1,2,1:2 
termdepmat2 <- function(f) {
  filist <- formula2idx(f, asCall = FALSE)
  n <- length(filist)
  v <- vapply(filist, paste0, "", collapse = ":")
  mat <- matrix(FALSE, n, n, dimnames = list(v, v))
  for(i in seq.int(n)) mat[match(expandintr(filist[[i]]), v), i] <- TRUE
  mat
}

## combines term-dependency-matrices
#termdepmat_list <- function(fl) 
#	termdepmat_combine(lapply(fl, termdepmat))

termdepmat_combine <- function(x) {
  dm <- sum(vapply(x, nrow, 1L))
  mat <- matrix(FALSE, dm, dm)
  j <- 1L
  for(i in seq_along(x)) {
    n <- nrow(x[[i]])
    k <- seq.int(j, length.out = n)
    mat[k, k] <- x[[i]]
    j <- j + n
  }
  dn <- unlist(lapply(x, rownames))
  dimnames(mat) <- list(dn, dn)
  mat
}

# converts formula to a(n unevaluated) list of numeric indices
# e.g. a*b --> list(1,2,1:2)
formula2idx <- function(x, asCall = TRUE) {
  if(!is.call(x) || !inherits(x, "formula")) stop("'x' is not a formula")
  fac <- attr(delete.response(terms(x)), "factors")
  dimnames(fac) <- NULL
  ret <- apply(fac > 0L, 2L, which)
  if(asCall) as.call(c(as.name("list"), ret)) else ret 
}


sample_fun <- function(i, dfx){
  (dfx %>% select(study, SID, p_value) %>%
     group_by(study) %>%
     mutate(p_value = sample(p_value, n(), replace = F)) %>%
     ungroup())$p_value 
}

perm_fun <- function(i, df2, f1, po){
  print(i)
  x <- df2 %>% select(-o_value) %>% mutate(o_value = factor(po))
  start <- Sys.time()
  fit2  <- glmer(formula(f1), data = x, na.action = "na.omit", family = "binomial")
  print(Sys.time() - start)
  par2 <- fixef(fit2)["p_value"]
  ci2 <- confint.merMod(fit2, method = "Wald")["p_value",]
  tidy2 <- c(par2, ci2)
  print(paste(i, "again!"))
  rm(list = c("x", "df2", "fit2"))
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
             "| study", collapse = "") # random effects
}

sca_setup_fun <- function(trait, outcome, f, fixed){
  terms <- sca_term_fun(trait, outcome)
  f <- sca_formula_fun(terms)
  # wd <- "~/selection"
  load(sprintf("%s/data/sca/sca_%s_%s.RData", wd, trait, outcome))
  df1 <- df1 %>% mutate(n = 1:n())
  std <- df1 %>%
    group_by(study, o_value) %>%
    tally() %>%
    full_join(crossing(study = unique(.$study), o_value = factor(c(0,1))))
  std <- unique(std$study[std$n < 50 | is.na(std$n)])
  df1 <- df1 %>% 
    filter(!study %in% std) %>%
    mutate(o_value = as.numeric(as.character(o_value))) 
  df2 <- df1 %>%
    select(study, SID, p_value, o_value, n, one_of(terms)) %>%
    filter(complete.cases(.)) %>%
    mutate_if(is.factor, factor)
  
  
  start <- Sys.time()
  cntrl <- glmerControl(optimizer = "bobyqa", boundary.tol = 0)#,
  # optCtrl=list(maxfun=2e5))
  full_mod <- femlm(formula(f), data = df2, family = "logit")
  # full_mod <- glmer(formula(f), data = df1, family = "binomial", na.action = "na.omit", control = cntrl)
  print(Sys.time() - start)
  gmCall <- get_call(full_mod)
  gmEnv <- parent.frame()
  gmFormulaEnv <- environment(as.formula(formula(full_mod), env = gmEnv))
  
  # set up terms
  
  allTerms <- allTerms0 <- getAllTerms.fixest(full_mod, terms)
  deps <- attr(allTerms0, "deps")
  fixed <- union(fixed, rownames(deps)[rowSums(deps, na.rm = TRUE) == 
                                         ncol(deps)])
  fixed <- c(fixed, allTerms[allTerms %in% "(Intercept)"])
  nFixed <- length(fixed)
  gmCoefNames <- gmCoefNames0 <- names(coef(full_mod))
  
  # this chunk makes sure that interactions include main effects
  varsOrder <- order(allTerms %in% fixed)
  termsOrder <- order(gmCoefNames %in% fixed)
  allTerms <- allTerms[varsOrder]
  gmCoefNames <- gmCoefNames[termsOrder]
  di <- match(allTerms, rownames(deps))
  deps <- deps[di, di, drop = FALSE]
  
  # setup output
  nVars <- length(allTerms)
  nTerms <- length(gmCoefNames)
  # rvNcol <- nVars + 3L + 2
  # rtNcol <- nTerms + 3L + 2
  # resultChunkSize <- 1000000L
  # rval <- matrix(NA_real_, ncol = rtNcol, nrow = resultChunkSize)
  
  nFixed <- 2  # p_value & intercept
  nVariants <- 1L
  nov <- as.integer(nVars - nFixed)
  ncomb <- (2L^nov) 
  comb.sfx <- rep(TRUE, nFixed)
  comb.seq <- if (nov != 0L) {seq_len(nov)}
  # ord <- integer(resultChunkSize)
  
  pval_out <- sapply(1:500, function(x){sample_fun(x, df1)}) # change back to 500
  # perm_res <- array(NA_real_, c(10, 3, resultChunkSize)) # change back to 500
  # no_cores <- detectCores() - 1
  
  k <- 0L
  save(list = ls(all.names = TRUE), 
       file = sprintf("%s/results/sca_workspace/sca_workspace_%s_%s.RData", wd, outcome, trait))
  return(ncomb)
}

sca_nested <- crossing(
  Trait = unique(p_waves$p_item),
  Outcome = (codebook$codebook[[2]] %>% filter(category == "out"))$name
) %>%
  filter(Outcome == "frstjob") %>%
  mutate(terms = map2(Trait, Outcome, sca_term_fun),
         formula = map(terms, sca_formula_fun),
         sca = pmap(list(Trait, Outcome, formula, "p_value"), sca_setup_fun))
