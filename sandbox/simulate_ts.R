###############################
#Simulate binary time series.
###############################

library(simcausal)
options(simcausal.verbose=FALSE)

t.end <- 1000

#No missing values, keep it very simple
D <- DAG.empty()
D <- D +

  #Initialize first points with probability 0.5
  node("W", t = 0, distr = "rbern",prob = 0.5) +
  node("A", t = 0, distr = "rbern",prob = 0.5) +
  node("Y", t = 0, distr = "rbern",prob = 0.5)

# Define very simple dependence-
# W(t) should depend on Y(t-1) and A(t-1)
# A(t) should depend on W(t) and Y(t-1)
# Y(t) should depend on A(t) and W(t)

D <- D +
  node("W", t = 1:t.end, distr = "rbern",
       prob = plogis(-0.3 + 1.2*A[t-1] - 0.55*Y[t-1])) +
  node("A", t = 1:t.end, distr = "rbern",
       prob = plogis(0.4 - 0.95*W[t] + 1.3*Y[t-1])) +
  node("Y", t = 1:t.end, distr = "rbern",
       prob = plogis(1.2 - 1.8*W[t]+ 1.4*A[t]))

lDAG <- set.DAG(D)

# take a draw of 1 from this DGP
data <- sim(DAG = lDAG, n = 1, rndseed = 123)
sim_ts_1000 <- data.frame(data=t(data[-1]))



