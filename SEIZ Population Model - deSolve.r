#load deSolve package
library(deSolve)

#SEIZ model
SEIZ <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- A - beta*S*I/N - b*S*Z/N - mu*S
    dE <- (1-p)*beta*S*I/N + (1-l)*b*S*Z/N - rho*E*I/N - eps*E - mu*E 
    dI <- p*beta*S*I/N + rho*E*I/N + eps*E - mu*I
    dZ <- l*b*S*Z/N - mu*Z
    return(list(c(dS, dE, dI, dZ)))
  })
}

#R0 <- function(parameters) with(as.list(parameters), beta/(mu+gamma))

# Set Parameters A ----------------------------------------------------------

N=600

# time at which we want the solution to be outputted
times  <- seq(0,300,by = 1)

#create list with the parameters
parameters <- c(A=30, mu=0.05, eps=0.15, beta=0.2, rho=0.2, b=0.17, l=0.9, p=0.15)

# create list with the initial conditions
initial.state <- c(S=560, E=5, I=20, Z=15)

# Results ---------------------------------------------------------------

# Solve using ode function from deSolve package
out_A <- ode(y = initial.state, times, SEIZ, parameters)

# change to data frame
out_A <- as.data.frame(out_A)

# Set Parameters B ----------------------------------------------------------

# time at which we want the solution to be outputted
times  <- seq(0,300,by = 1)

#create list with the parameters
parameters <- c(A=30, mu=0.05, eps=0.15, beta=0.2, rho=0.2, b=0.2, l=0.99, p=0.15)

# create list with the initial conditions
initial.state <- c(S=560, E=5, I=20, Z=15)

# Results ---------------------------------------------------------------

# Solve using ode function from deSolve package
out_B <- ode(y = initial.state, times, SEIZ, parameters)

# change to data frame
out_B <- as.data.frame(out_B)

# Set Parameters C ----------------------------------------------------------

# time at which we want the solution to be outputted
times  <- seq(0,300,by = 1)

#create list with the parameters
parameters <- c(A=30, mu=0.05, eps=0.15, beta=0.08, rho=0.15, b=0.08, l=0.87, p=0.15)

# create list with the initial conditions
initial.state <- c(S=560, E=5, I=20, Z=15)

# Results ---------------------------------------------------------------

# Solve using ode function from deSolve package
out_C <- ode(y = initial.state, times, SEIZ, parameters)

# change to data frame
out_C <- as.data.frame(out_C)

# Plotting A-B-C ----------------------------------------------------------------

viz_A <- par(fig=c(0,1,0.5,1),mar=c(4,4,2,5))
par(mfrow=c(1,1))
plot(S~time,data=out_A,type='l', ylim=c(0,600), xlab='Time', main=paste('SEIZ with Parameters A'), ylab='populations')
lines(E~time,data=out_A,type='l', col='green'); par(new=TRUE)
lines(I~time,data=out_A,type='l', col='blue'); par(new=TRUE)
lines(Z~time,data=out_A,type='l', col='red'); par(new=TRUE)
legend('topright', legend=c('Susceptible', 'Exposed', 'Infectious', 'Skeptical'), col=c('black','green','blue','red'), lty=1, par(viz_A), cex=0.8)

viz_B <- par(fig=c(0,1,0.5,1),mar=c(4,4,2,5))
par(mfrow=c(1,1))
plot(S~time,data=out_B,type='l', ylim=c(0,600), xlab='Time', main=paste('SEIZ with Parameters B'), ylab='populations')
lines(E~time,data=out_B,type='l', col='green'); par(new=TRUE)
lines(I~time,data=out_B,type='l', col='blue'); par(new=TRUE)
lines(Z~time,data=out_B,type='l', col='red'); par(new=TRUE)
legend('topright', legend=c('Susceptible', 'Exposed', 'Infectious', 'Skeptical'), col=c('black','green','blue','red'), lty=1, par(viz_B), cex=0.8)

viz_C <- par(fig=c(0,1,0.5,1),mar=c(4,4,2,5))
par(mfrow=c(1,1))
plot(S~time,data=out_C,type='l', ylim=c(0,600), xlab='Time', main=paste('SEIZ with Parameters C'), ylab='populations')
lines(E~time,data=out_C,type='l', col='green'); par(new=TRUE)
lines(I~time,data=out_C,type='l', col='blue'); par(new=TRUE)
lines(Z~time,data=out_C,type='l', col='red'); par(new=TRUE)
legend('topright', legend=c('Susceptible', 'Exposed', 'Infectious', 'Skeptical'), col=c('black','green','blue','red'), lty=1, par(viz_C), cex=0.8)


#SEIZ model II
SEIZ_2 <- function(time, state, parameters) {
	with(as.list(c(state, parameters)), {
		dS <- A - beta*S*I/N - b*S*Z/N - mu*S
		dE <- (1-p)*beta*S*I/N + (1-l)*b*S*Z/N - rho*E*I/N - eps*E - mu*E - phi*E - rho*E*Z/N
		dI <- p*beta*S*I/N + rho*E*I/N + eps*E - mu*I
		dZ <- l*b*S*Z/N - mu*Z + phi*E + rho*E*Z/N
		return(list(c(dS, dE, dI, dZ)))
	})
}

# Set Parameters A ----------------------------------------------------------

N=600

# time at which we want the solution to be outputted
times_2  <- seq(0,300,by = 1)

#create list with the parameters
parameters_2 <- c(A=30, mu=0.05, phi = 0.10 ,eps=0.50, beta=0.2, rho=0.15, b=0.17, l=0.9, p=0.15)

# create list with the initial conditions
initial.state_2 <- c(S=560, E=5, I=20, Z=15)

# Results ---------------------------------------------------------------

# Solve using ode function from deSolve package
out_2 <- ode(y = initial.state_2, times_2, SEIZ_2, parameters_2)

# change to data frame
out_2 <- as.data.frame(out_2)

# Other option for plotting 
viz_2 <- par(fig=c(0,1,0.5,1),mar=c(4,4,2,5))
par(mfrow=c(1,1))
plot(S~time,data=out_2,type='l', ylim=c(0,600), xlab='Time', main=paste('SEIZ with Alternate Schema'), ylab='populations')
lines(E~time,data=out_2,type='l', col='green'); par(new=TRUE)
lines(I~time,data=out_2,type='l', col='blue'); par(new=TRUE)
lines(Z~time,data=out_2,type='l', col='red'); par(new=TRUE)
legend('topright', legend=c('Susceptible', 'Exposed', 'Infectious', 'Skeptical'), col=c('black','green','blue','red'), lty=1, par(viz_2), cex=0.8)

