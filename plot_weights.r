information = function(x) {
  return((-1)*x*log2(x))
}

weighted = function(x) {
  return(x)
}

unweighted = function(x) {
  if (x == 0) {
    return(0)
  }
  else {
    return(1)
  }
}

# geometric mean is assumed to be 0.01

centered_ratio = function(x) {
  return(log2(x) - log2(0.01))
}

xvals <- c(0:1000)
xvals <- xvals / 1000

plot.new()

lines(xvals, sapply(xvals,unweighted,simplify=TRUE), col = "black",xlab = "proportional abundance")
lines(xvals, sapply(xvals,weighted,simplify=TRUE), col = "red",xlab = "proportional abundance")
lines(xvals, sapply(xvals,information,simplify=TRUE), col = "blue",xlab = "proportional abundance")
lines(xvals, sapply(xvals,centered_ratio,simplify=TRUE), col = "orange",xlab = "proportional abundance")
axis(1)
axis(2)