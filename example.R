source('sir.R')

# generate vaccination rates
r0vals <- seq(0, 100, by = 1)

# compute the number of infected individuals
rtvals <- generateRt(r0vals, i0 = 5, beta = 0.12, gamma = 0.1)

plot(r0vals, costFunc(r0vals, rtvals, rel_cost = 5))
