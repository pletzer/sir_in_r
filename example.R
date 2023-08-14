source('sir.R')

# population size
N = 100

# i0 = number of infected people (initially)
i0 = 5


# generate vaccination rates
r0vals <- seq(0, N - i0, by = 1)

# compute the number of infected individuals
rtvals <- generateRt(r0vals, i0 = i0, beta = 0.12, gamma = 0.1, N = N)

plot(r0vals, costFunc(r0vals, rtvals, rel_cost = 5))
