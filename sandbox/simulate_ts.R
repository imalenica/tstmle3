###############################
#Simulate time-series data
###############################

#Scenario 1a:
#Markov order 1, binary Y and A, with with several W (binary, continuous)

library(simcausal)
options(simcausal.verbose=FALSE)

t.end <- 500

#No missing values, keep it very simple
D <- DAG.empty()
D <- D +

  #Initialize first points
  node("A", t = 0, distr = "rbern",prob = 0.5) +
  node("Y", t = 0, distr = "rbern",prob = 0.5) +
  node("W1", t = 0, distr = "rbern",prob = 0.5) +
  node("W2", t = 0, distr = "runif", min = 0, max = 3) +
  node("W3", t = 0, distr = "runif", min = 0, max = 1)

# Define very simple dependence-
# A(t) should depend on W(t-1) and Y(t-1)
# Y(t) should depend on A(t) and W(t-1)
# W1(t) here depends on A(t) and Y(t), or independent.

D <- D +
  node("A", t = 1:t.end, distr = "rbern",
       prob = plogis(0.4 + 0.25*W1[t-1] - 0.45*W2[t-1] + 0.3*Y[t-1])) +
  node("Y", t = 1:t.end, distr = "rbern",
       prob = plogis(0.3 - 0.8*W1[t-1]+ 0.4*W3[t-1] + 1.2*A[t])) +
  node("W1", t = 1:t.end, distr = "rbern",
       prob = plogis(0.3 - 0.5*Y[t] + 0.6*A[t])) +
  node("W2", t = 1:t.end, distr = "runif", min = 0, max = 3) +
  node("W3", t = 1:t.end, distr = "runif", min = 0, max = 1)

lDAG <- set.DAG(D)

# take a draw of 1 from this data generating distribution
data <- sim(DAG = lDAG, n = 1, rndseed = 123)
sim_ts_s1 <- data.frame(data=t(data[-1]))

#Check out distributions of A, Y, W:
W1 <- sim_ts_s1[grep("W1_", row.names(sim_ts_s1), value = TRUE), ]
A <- sim_ts_s1[grep("A", row.names(sim_ts_s1), value = TRUE), ]
Y <- sim_ts_s1[grep("Y", row.names(sim_ts_s1), value = TRUE), ]
table(W1)
table(A)
table(Y)

#Scenario 1b:
#Markov order 1, binary Y and A, with with several W (binary, continuous)

library(simcausal)
options(simcausal.verbose=FALSE)

t.end <- 50

#No missing values, keep it very simple
D <- DAG.empty()
D <- D +

  #Initialize first points
  node("A", t = 0, distr = "rbern",prob = 0.5) +
  node("Y", t = 0, distr = "rbern",prob = 0.5) +
  node("W1", t = 0, distr = "rbern",prob = 0.5) +
  node("W2", t = 0, distr = "runif", min = 0, max = 3) +
  node("W3", t = 0, distr = "runif", min = 0, max = 1)

# Define very simple dependence-
# A(t) should depend on W(t-1) and Y(t-1)
# Y(t) should depend on A(t) and W(t-1)
# W1(t) here depends on A(t) and Y(t), or independent.

D <- D +
  node("A", t = 1:t.end, distr = "rbern",
       prob = plogis(0.4 + 0.25*W1[t-1] - 0.45*W2[t-1] + 0.3*Y[t-1])) +
  node("Y", t = 1:t.end, distr = "rbern",
       prob = plogis(0.3 - 0.8*W1[t-1]+ 0.4*W3[t-1] + 1.2*A[t])) +
  node("W1", t = 1:t.end, distr = "rbern",
       prob = plogis(0.3 - 0.5*Y[t] + 0.6*A[t])) +
  node("W2", t = 1:t.end, distr = "runif", min = 0, max = 3) +
  node("W3", t = 1:t.end, distr = "runif", min = 0, max = 1)

lDAG <- set.DAG(D)

# take a draw of 1 from this data generating distribution
data <- sim(DAG = lDAG, n = 1, rndseed = 123)
sim_ts_s1_n50 <- data.frame(data=t(data[-1]))

#Check out distributions of A, Y, W:
W1 <- sim_ts_s1_n50[grep("W1_", row.names(sim_ts_s1_n50), value = TRUE), ]
A <- sim_ts_s1_n50[grep("A", row.names(sim_ts_s1_n50), value = TRUE), ]
Y <- sim_ts_s1_n50[grep("Y", row.names(sim_ts_s1_n50), value = TRUE), ]
table(W1)
table(A)
table(Y)
