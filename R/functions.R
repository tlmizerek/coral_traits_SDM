change_CTDB <- function(x) {
  if (length(x) > 1) {      # Does a species by trait combination have more than 1 value?
    x <- type.convert(x, as.is=TRUE)
    if (is.character(x)) {
      return(x[1])          # If values are strings (characters), then return the first value
    } else {
      return(as.character(max(x)))  # If values are numbers, then return the max (converted back to character). for depth range, we're interested in max.
    }
  } else {
    return(x) # If a species by trait combination has 1 value, then just return that value 
  }
}


#build.database <- function(veron, filename.CTDB) {
#CTDB <- read.csv(filename.CTDB, as.is=TRUE)
#CTDB <- read.csv("data/ctdb_20140414.csv", as.is=TRUE)
#"data/ctdb_20140414.csv"}
#CTDB[is.na(CTDB) == T] = 0
#CTDB[CTDB == "unknown"] <- NA

my_aggregate_rules <- function(x) {
    if (length(x) > 1) {               # Does a species by trait combination have more than 1 value?
        x <- type.convert(x, as.is=TRUE)
        if (is.character(x)) {
            return(x[1])                   # If values are strings (characters), then return the first value
        } else {
            return(as.character(max(x)))  # If values are numbers, then return the max (converted back to character). for depth range, we're interested in max.
        }
    } else {
        return(x)                        # If a species by trait combination has 1 value, then just return that value
    }
}

classify <- function(x, tr) {
  unname(tr[match(x, names(tr))])
}


make.transparent  <-  function(col, opacity=0.5) {
  if(length(opacity) > 1 && any(is.na(opacity))) {
    n       <-  max(length(col), length(opacity))
    opacity <-  rep(opacity, length.out=n)
    col     <-  rep(col, length.out=n)
    ok      <-  !is.na(opacity)
    ret     <-  rep(NA, length(col))
    ret[ok] <-  Recall(col[ok], opacity[ok])
    ret
  } else {
    tmp  <-  col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
  }
}

to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

to.pdf <- function(expr, filename, ...) {
  to.dev(expr, pdf, filename, ...)
}

label <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
  usr  <-  par("usr")
  x.p  <-  usr[1] + px*(usr[2] - usr[1])
  y.p  <-  usr[3] + py*(usr[4] - usr[3])
  if(log=="x"){x.p<-10^(x.p)}
  if(log=="y"){y.p<-10^(y.p)}
  if(log=="xy"){x.p<-10^(x.p);y.p<-10^(y.p)}
  if(text){
    text(x.p, y.p, lab, adj=adj, ...)
  } else {
    points(x.p, y.p, ...)
  }
}

#create new region column
assignNewRegion  <-  function(data, chosenColumns, newRegion, tag, finalValue, startValue=NA) {
  newColumn          <-  paste0(tag, newRegion)
  data[[newColumn]]  <-  startValue
  i  <-  rowSums(data[chosenColumns]) >= 1
  if(any(i))
    data[[newColumn]][i]  <-  finalValue
  data
}


#variance inflation factor mixed models
vif.mer <- function (fit) {
    ## adapted from rms::vif
    v <- vcov(fit)
    nam <- names(fixef(fit))
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
}

ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
  require(ggplot2)
  f <- function(x) {
    pv   <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
    pDf  <- data.frame(y=unlist(x)[ord],
                       ci=1.96*se[ord],
                       nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                       ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                       ind=gl(ncol(x), nrow(x), labels=names(x)))
    
    if(QQ) {  ## normal QQ-plot
      p <- ggplot(pDf, aes(nQQ, y))
      p <- p + facet_wrap(~ ind, scales="free")
      p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
    } else {  ## caterpillar dotplot
      p <- ggplot(pDf, aes(ID, y)) + coord_flip()
      if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
        p <- p + facet_wrap(~ ind)
      } else {           ## different scales for random effects
        p <- p + facet_grid(ind ~ ., scales="free_y")
      }
      p <- p + xlab("Levels") + ylab("Random effects")
    }
    
    p <- p + theme(legend.position="none")
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
    p <- p + geom_point(aes(size=1.2), colour="blue") 
    return(p)
  }
  
  lapply(re, f)
}