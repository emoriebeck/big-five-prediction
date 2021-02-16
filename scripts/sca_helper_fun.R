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

matchCoef <- function (m1, m2, all.terms = getAllTerms(m2, intercept = TRUE), 
                       beta = 0L, terms1 = getAllTerms(m1, intercept = TRUE), coef1 = NULL, 
                       allCoef = FALSE, ...) 
{
  if (is.null(coef1)) {
    ct <- if (beta != 0L) 
      std.coef(m1, beta == 2L, ...)
    else coefTable(m1, ...)
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
