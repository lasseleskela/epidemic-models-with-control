
# SIR control functions
# Lasse Leskelä 2022-02-15
#

library(deSolve)
library(stats)


# 0 Basic path characteristics

FinalState   <- function(path) { return(tail(path,1)[1,c("S","I")]) }
InitialState <- function(path) { return(path[1,c("S","I")]) }
PeakISize    <- function(path) { return(max(path[,"I"])) }
PeakITime    <- function(path) { return(unname(path[which.max(path[,"I"]),"time"])) }


# 1 ODE models

# 1.1 SIR with no control
#
# Input:
# state = (S,I)
# par = (beta,gamma)
#
SIR <- function (t, state, par) {
  with(as.list(c(state, par)), {
    dS <- -beta*S*I
    dI <- beta*S*I - gamma*I
    list(c(dS, dI))
  })
}

# 1.1.1 Example
epar <- list(beta=0.6, gamma=0.2)
init <- list(S=0.9999, I=0.0001)
dt <- epar$gamma/10
times <- seq(0,100,by=dt)
path  <- ode(c(S=init$S, I=init$I), times, SIR, epar)
plot(path[,"time"], path[,"S"], ylim=c(0,1), type="l", col="blue")
lines(path[,"time"], path[,"I"], type="l", col="red")

# 1.2 SIR with a single constant-level lockdown
# with start time "start", duration "dur", and level "lev".
#
# Input:
# state = (S,I)
# par = (beta, gammma, start, dur, lev)
#
SIR.1LD <- function (t, state, par) {
  with(as.list(c(state, par)), {
    u <- (start<t)*(t<=(start+dur))*lev
    dS <- -beta*(1-u)*S*I
    dI <- beta*(1-u)*S*I-gamma*I
    list(c(dS, dI))
  })
}

# 1.2.1 Example
epar <- list(beta=0.6, gamma=0.2)
cpar <- list(start=15, dur=20, lev=0.75)
init <- list(S=0.9999, I=0.0001)
dt <- epar$gamma/10
times <- seq(0,100,by=dt)
path <- ode(c(S=init$S, I=init$I), times, SIR.1LD, c(epar, cpar))
plot(path[,"time"], path[,"S"], ylim=c(0,1), type="l", col="blue")
lines(path[,"time"], path[,"I"], type="l", col="red")
abline(v=cpar$start,lty=2)
abline(v=cpar$start+cpar$dur,lty=2)


# 1.3 Wait-maintain-relax control
# with critical peak prevalence b, start time t1, end time t2
#
SIR.WMR <- function(t, state, par) {
  with(as.list(c(state, par)), {
    u  <- (t>t1)*(t<t2)*(1 - 1/(1+b*beta*(t2-t)))
    dS <- -(1-u)*beta*S*I
    dI <- (1-u)*beta*S*I - gamma*I
    list(c(dS, dI))
  })
}

# 1.3.1 Example
epar <- list(beta=0.6, gamma=0.2)
cpar <- list(b=0.075, t1=17.0, t2=53.7)
init <- list(S=0.9999, I=0.0001)
dt <- epar$gamma/10
times <- seq(0,100,by=dt)
path <- ode(c(S=init$S, I=init$I), times, SIR.WMR, c(epar, cpar))
plot(path[,"time"], path[,"S"], ylim=c(0,1), type="l", col="blue")
lines(path[,"time"], path[,"I"], type="l", col="red")
abline(v=cpar$t1,lty=2)
abline(v=cpar$t2,lty=2)
abline(h=cpar$b, lty=2)


# 2. Basic SIR semigroup operator
#
# Input:
# epar = (beta, gamma) Transmission parameters (as named list)
# init = (S,I)         Initial state (as named list)
# t                    Future time instant, possibly infinite
# 
# Output:
# final = (S,I)        State after t time units (as named list)
#
# Auxiliary input:
# res = Time resolution (default 10); time slots per infectious period
#
# Implementation:
# deSolve::ode() handles transforming the initial state and the output
# state from named list to named numeric vector.
#
SIR.evolve <- function (epar, init, t, res=10) {
  if (t==0) {
    return(init)
  } else if (t==Inf) {
    rho <- epar$beta/epar$gamma
    v0 <- init$S + init$I - log(init$S)/rho
    f <- function(s) { s - log(s)/rho - v0}
    S.inf <- uniroot(f, interval=c(0,1/rho) )$root
    return(list(S=S.inf, I=0))
  } else {    
    times <- seq(0,t,by=epar$gamma/res)
    path  <- ode(c(S=init$S, I=init$I), times, SIR, epar)
    final <- list(S=tail(path,1)[,"S"], I=tail(path,1)[,"I"])
  }
}



# 3 Limiting susceptible shares

# 3.1 Final susceptible share given no control
#
# Input:
# epar = (beta, gamma) Transmission parameters (as named list)
# init = (S,I)         Initial state (as named list)
# 
# Output: Final susceptible share
#
final.s <- function (epar, init) {
  output <- SIR.evolve(epar, init, Inf)$S
}

# 3.1.1 Example
epar <- list(beta=0.6, gamma=0.2)
init <- list(S=0.9999, I=0.0001)
S.Inf <- final.s(epar, init)
print(paste0("Limiting susceptible share equals ", round(S.Inf,3)))



# 3.2 Final susceptible share given 1 lockdown.
#
# Input:
# epar = (beta, gamma)      Transmission parameters (named list)
# cpar = (start, dur, lev)  Start time, duration, and level of lockdown (named list)
# init = (S,I)              Initial state (named list)
# 
# Output: Final susceptible share
#
final.s.LD1 <- function(epar, cpar, init, res=10) {
 epar.ld <- list(beta = epar$beta*(1-cpar$lev), gamma = epar$gamma)
 x0 <- init
 x1 <- SIR.evolve(epar, x0, cpar$start, res)
 x2 <- SIR.evolve(epar.ld, x1, cpar$dur, res)
 x3 <- SIR.evolve(epar, x2, Inf)
 output <- x3$S
}

# 3.1.2 Example
epar <- list(beta=0.6, gamma=0.2)
init <- list(S=0.9999, I=0.0001)
S.Inf <- final.s.LD1(epar, cpar, init)
print(paste0("Limiting susceptible share equals ", round(S.Inf,3)))



# 4 Time to herd immunity given no control
#
# Input:
# epar = (beta, gamma)      Transmission parameters (named list)
# init = (S,I)              Initial state (named list)
# 
# Output: Time to reach herd immunity
#
TimeToHerdImmunity <- function(epar, init, res=10) {
 beta  <- epar$beta
 gamma <- epar$gamma
 s0 <- init$S
 i0 <- init$I
 rho   <- beta/gamma
 dt    <- gamma/res
 t0.max <- log(rho*s0)/(beta*i0)  # Universal upper bound of uncontrolled herd immunity time
 times <- seq(0,t0.max,by=dt)
 path0 <- ode(c(S=init$S, I=init$I), times, SIR, epar)
 t.peak0 <- PeakITime(path0)
 if (t.peak0 >= t0.max) {
   print("ERROR")
   return(-1)
 }
 return(t.peak0)
}

# 4.1.1 Example
epar <- list(beta=0.6, gamma=0.2)
init <- list(S=0.9999, I=0.0001)
t.herd <- TimeToHerdImmunity(epar, init)
print(paste0("Time to herd immunity equals ", round(t.herd,3)))



# 5 Optimal start time for minimising total incidence with
# a single lockdown of duration "dur" and level "lev".
# 
# Output:
# objective = Minimum nonsusceptible share 1-S(infty)
# minimum   = Optimal start time
#
# Implementation:
# psi0() is first used to check whether the optimal start time is 0,
# as described in Bliman and Duprez 2021.
# If not, then the R optimize() is used to search for an optimum.
# This uses the "optimise" function in R "stats" package.
#
psi0 <- function(epar, init, dur, lev, res=10) {
 S0 <- init$S
 I0 <- init$I
 dt <- epar$gamma/res
 epar.ld <- list(beta=epar$beta*(1-lev), gamma=epar$gamma)
 path <- ode(c(S=S0,I=I0), times=seq(0,dur,dt), SIR, epar.ld)
 I <- path[,"I"]
 I1 <- tail(I,1)
 integral <- I1*sum(1/I)*dt # [Eq 4; Bliman Duprez 2021]
 psi0 <- -I1/I0 - lev * epar$gamma * integral + 1 
}
Optimal.start.LD1 <- function(epar, init, dur, lev, res=10) {
 psi0.val <- psi0(epar, init, dur, lev)
 if (psi0.val >= 0) {
  val <- 1 - final.s.LD1(epar, list(start=0, dur=dur, lev=lev), init, res)
  opt <- list(minimum=0, objective=val)
 } else {
  t.herd <- TimeToHerdImmunity(epar, init, res)
  f <- function (x) { 1 - final.s.LD1(epar, list(start=x, dur=dur, lev=lev), init) }
  opt <- optimize(f, lower=0, upper=t.herd) # Minimise final size 1-S(infty)
 }
 return(opt)
}

# 5.1.1 Example
epar <- list(beta=0.6, gamma=0.2)
init <- list(S=0.9999, I=0.0001)
dur <- 20
lev <- 0.75
opt <- Optimal.start.LD1(epar, init, dur, lev)
t.opt <- opt$minimum
print(paste0("Optimal start time equals ", round(t.opt,3)))
cpar <- list(start=t.opt, dur=dur, lev=lev)
S.Inf <- final.s.LD1(epar, cpar, init)
print(paste0("Limiting susceptible for optimally started lockdown share equals ", round(S.Inf,3)))

