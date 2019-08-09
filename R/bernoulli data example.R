x1 = runif(100)
beta = c(1,2)
eta = beta[1]+beta[2]*x1
y1 = rbinom(100,1,exp(eta)/(1+exp(eta)))

kette = metrohas(y1~x1, dist= "bernoulli", beta_start = c(0,0),anzahl_sim = 10000)
density(kette[500:10000,1]) # mean(0.823)
density(kette[500:10000,2]) # mean(0.84)

glm_model = glm(y1~x1, family = binomial(link="logit"))
glm_model   #intercept 0.84 #slope 1.10
