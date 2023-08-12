library(deSolve)


sirModel <- function(r0, i0, beta, gamma, N, times) {
    #
    # Solves the SIR model
    #param: r0, number of people vaccinated (initially)
    #param: i0, initial number of people with an infection
    #param: beta, infection rate 1/(time it takes to infect another person)
    #param: gamma, 1/(time it takes to recover)
    #param: N, number of subjects
    #param: times, steps the ODE integrator should take, we're only interested in the last value

    parameters <- c(beta = beta, gamma = gamma, N = N)

    state <- c(S = N - r0 - i0, I = i0, R = r0)

    tendency <- function(t, state, parameters) {
      with(as.list(c(state, parameters)),{
        # rate of change
        dSdt <- -beta * S * I / N
        dIdt <- beta * S * I / N - gamma * I
        dRdt <- gamma * I 
        # return the rate of change
        list(c(dSdt, dIdt, dRdt))
     }) # end with(as.list ...
    }
    
    out <- ode(y = state, times = times, func = tendency, parms = parameters)

    n = length(times)

    # R(T)
    return(out[n, 4])
}


costFunc <- function(r0, rt, rel_cost){
    #param: r0, number of people vaccinated at time 0
    #param: rt, number of people recovered at time T (the end of the period) - the number of people who have been infected (vaxxed or disease induced immunity)
    #param: rel_cost, average cost of treatment per person relative to cost of vaccination per person

    cost <- r0 + rel_cost * (rt - r0)

    return(cost)
}


generateRt <- function(r0vals, i0, beta, gamma) {
    rt <- as.numeric(lapply(r0vals, sirModel, i0, beta, gamma, N = 100, times = seq(0., 1000, by = 1)))
    return(rt)
}

# Example
# source('sir.R')
# r0vals <- seq(0, 100, by = 1)
# rtvals <- generateRt(r0vals, i0 = 5, beta = 0.12, gamma = 0.1)
# plot(r0vals, costFunc(r0vals, rtvals, rel_cost = 5))






