# routines for quantization and one-dimensional clustering
# Tom Minka 11/30/01

eps <- 1e-15

# returns the elements of x not indexed by i
not <- function(x,i) {
  if(is.numeric(i)) return(x[-i])
  if(is.character(i)) return(x[setdiff(names(x),i)])
  if(is.logical(i)) return(x[!i])
}

##############################################################################

ncut <- function(x,b) {
  cut(x,b,labels=FALSE,include.lowest=TRUE)
}

# reduce data to n points, by local averaging
# use quantiles so size of each bin is equal
squash <- function(x,n=1000) {
  b <- floor(length(x)/n)
  rb <- rep(b,n)
  # spread the overflow broadly over the range
  overflow <- length(x) - n*b
  rb <- rb + c(0,diff(trunc(overflow/n*(1:n))))
  # f is 1 1 1 2 2 2 3 3 3 ...
  f <- rep(1:n, rb)
  f <- f[1:length(x)]
  x <- sort(x)
  tapply(x,f,mean)
}

plot.breaks <- function(x,b,scale=1) {
  hb <- break.equal(x,40)
  h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
  top <- max(h$density)
  par(lwd=2)
  segments(b,rep(0,length(b)), b,rep(top*scale,length(b)), col="blue")
  par(lwd=1)
}

plot.ts.breaks <- function(x,b,scale=1) {
  plot(x,type="l")
  r <- range(x)*scale
  par(lwd=2)
  segments(b,rep(r[1],length(b)), b,rep(r[2],length(b)), col="blue")
  par(lwd=1)
}

plot.hclust.breaks <- function(hc,x,n=2:5) {
  hb <- break.equal(x,40)
  h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
  top <- max(h$density)
  par(lwd=3)
  for(i in 1:length(n)) {
    q <- cutree(hc,n[i])
    b <- break.cut(x,q)
    scale <- (length(n)-i+1)/length(n)
    segments(b,rep(0,length(b)), b,rep(top*scale,length(b)), col="blue")
  }
  par(lwd=1)
}

# evenly-spaced breaks
break.equal <- function(x,n=2) {
  r <- range(x)
  seq(0,1,by=1/n)*diff(r) + r[1]
}

# this is "equal-frequency" binning
break.quantile <- function(x,n=2,plot=TRUE) {
  b <- as.numeric(quantile(x, seq(0,1,by=1/n)))
  if(plot) plot.breaks(x,b)
  b
}

# convert a vector of cluster clusters into a break vector
break.centers <- function(x,m) {
  m <- sort(m)
  n <- length(m)
  b <- (m[1:n-1]+m[2:n])/2
  r <- range(x)
  c(r[1]-eps,b,r[2])
}

library(mva)
break.kmeans <- function(x,n=2,plot=T) {
  km <- kmeans(x,n)
  m <- as.numeric(km$centers)
  ss <- sum(km$withinss)
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.centers(x,m)
  if(plot) plot.breaks(x,b)
  b
}

# divide the range of x into n pieces, based on the largest spaces
# returns a vector of breaks, suitable for input to "hist" or "cut"
# tends to overdivide low-density regions
break.diff <- function(x,n=2) {
  n <- n-1
  len <- length(x)
  x <- sort(x)
  diff <- x[2:len] - x[1:len-1]
  i <- order(-diff)
  mid <- (x[2:len] + x[1:len-1])/2
  b <- sort(mid[i[1:n]])

  split.screen(c(2,1))
  split.screen(c(1,2),screen=1)
  screen(4)
  par(mar=c(2,4,1,1))
  last <- n+10
  plot(diff[i[1:last]],ylab="diff")
  m <- (diff[i[n]] + diff[i[n+1]])/2
  #arrows(1,m, last,m, col="black",code=0,lty=2)
  segments(n+0.5, diff[i[last]], n+0.5, diff[i[1]], col="black",lty=2)

  screen(2)
  par(mar=c(2,4,1,1))
  plot(x)
  segments(0,b, len+1,b, col="blue",lty=2)

  screen(3)
  par(mar=c(2,4,1,1))
  hist(x,30,col="bisque",main="")
  segments(b,rep(0,n), b,rep(100,n), col="blue")
  close.screen(all = TRUE)

  c(x[1]-eps,b,x[len])
}

# convert a cut of x into breaks
# f is a factor which divides x into disjoint groups
break.cut <- function(x,f) {
  low <- sort(as.numeric(tapply(x,f,min)))
  high <- sort(as.numeric(tapply(x,f,max)))
  n <- length(high)
  b <- (low[2:n]+high[1:n-1])/2
  c(low[1]-eps,b,high[n])
}

scatter <- function(x) {
  sum((x - mean(x))^2)
}

# "ward" is best
# "ave","centroid" are good
# "single", "complete","median","mcq" are bad
break.hclust2 <- function(x,n=2,method="ward") {
  h <- hclust(dist(x)^2,method)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  split.screen(c(2,1))
  plot.hclust.trace(h)
  screen(2)
  #plot.breaks(x,b)
  plot.hclust.breaks(h,x,n=2:n)
  close.screen(all=TRUE)
  b
}

# Create a hierarchy by Ward's method
# Merge clusters to minimize the sum of squares criterion
# x is a vector (one-dimensional clustering only)
# if x is a vector of batch means, n is the size of each batch
# and s is the sum of squares in each batch (only needed if same.var=F)
# if same.var=F, the adjacency constraint may give suboptimal results
sstree <- function(x,n=rep(1,length(x)),s=rep(1,length(x)),
                   sortx=TRUE,same.var=T) {
  # don't sort and it becomes time-series clustering
  if(sortx) {
    ord <- order(x)
    x <- sort(x)
  } else {
    ord <- 1:length(x)
  }
  height <- c()
  label <- -(1:length(x))
  label <- label[ord]
  merge <- matrix(0,length(x)-1,2)
  if(length(x) == 1) ks <- c()
  else ks <- 1:(length(x)-1)
  for(k in ks) {
    len <- length(x)
    # here we exploit the fact that the data is one-dimensional.
    # only clusters which are adjacent in sorted order need be considered
    # for merging.
    x1 <- x[1:len-1]
    x2 <- x[2:len]
    n1 <- n[1:len-1]
    n2 <- n[2:len]
    # merging cost
    cost <- n1*n2/(n1+n2)*(x1-x2)^2
    if(!same.var) {
      s1 <- s[1:len-1]
      s2 <- s[2:len]
      cost <- (n1+n2)*log((s1+s2+cost)/(n1+n2)) - n1*log(s1/n1) - n2*log(s2/n2)
    }
    j <- which.min(cost)

    if(!same.var) {
      s[j] <- s[j] + s[j+1] +
              n[j]*n[j+1]/(n[j]+n[j+1])*(x[j]-x[j+1])^2
      s <- not(s,j+1)
    }

    x[j] <- (x[j]*n[j] + x[j+1]*n[j+1])/(n[j] + n[j+1])
    n[j] <- n[j] + n[j+1]
    x <- not(x,j+1)
    n <- not(n,j+1)

    height[k] <- min(cost)
    merge[k,1] <- label[j]
    merge[k,2] <- label[j+1]
    label[j] <- k
    label <- not(label,j+1)
  }
  hc <- list(merge=merge,height=cumsum(height),order=ord,labels=NULL,
             method="ward")
  class(hc) <- "hclust"
  hc
}

plot.hclust.trace <- function(h,k=1:10) {
  g <- rev(h$height)
  k <- k[k < length(g)]
  #g <- sqrt(g)
  if(h$method != "single") {
    g <- -diff(g)
  }
  plot(k,g[k],ylab="merging cost",xlab="clusters")
}

break.hclust <- function(x,n=2,plot=T) {
  h <- sstree(x)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  if(plot) {
    split.screen(c(1,2))
    plot.hclust.trace(h)
    screen(2)
    #plot.breaks(x,b)
    plot.hclust.breaks(h,x,n=2:n)
    close.screen(all=TRUE)
  }
  b
}

# convert a cut of x into a changepoint plot
# f is a factor which divides x into disjoint groups
plot.chpt <- function(x,f) {
  plot(x,type="l",xlab="time",ylab="level")
  i <- 1:length(x)
  m <- tapply(x,f,mean)
  low <- tapply(i,f,min)
  high <- tapply(i,f,max)
  segments(low, m, high, m, col="blue")
}

# Apply Ward's method to find changepoints in a time-series
# x is a vector
break.chpt <- function(x,n=2,trace=TRUE,same.var=T) {
  h <- sstree(x,sortx=FALSE,same.var=same.var)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  cat("sum of squares =", format(ss), "\n")

  if(length(h$height) < 2) trace<-F
  if(trace) {
    split.screen(c(2,1))
    plot.hclust.trace(h)
    screen(2)
  }
  plot.chpt(x,q)
  points(x)
  if(trace) close.screen(all=TRUE)
  c(0,which(diff(q)>0),length(x)+1)
}
