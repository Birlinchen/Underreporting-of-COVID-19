success <- 0:20

plot(success, dbinom(success, size=250, prob=.01),type='h')

plot(success, dbinom(success, size=250, prob=.03),type='h')
