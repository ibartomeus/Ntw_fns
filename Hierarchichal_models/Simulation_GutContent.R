########################################
#Using hierarchichal models and Prey traits to estimate real predation 
#Author: Ignasi Bartomeus
#Chunks of code borrowed from Andy Royle and Richard Chandler Unmarked tutorial
#This code can also be used for plant-pollinator networks
########################################

# First let's assume gut content gives you perfect detection 
# of all links a species has in a nwteork. 

  # And let's assume that the number of links is a function of body size.
  # Create fake data
nPrey <- 20
BodySize <- runif(nPrey, 1, 3) # uniform from 1 to 3
  # Suppose that link probability increases with Body size
  # The relationship is described by an intercept of -3 and
  # a slope parameter of 2 on the logit scale
  # plogis is the inverse-logit (constrains us back to the [0-1] scale)
psi <- plogis(-3 + 2*BodySize)
  # Now we do gut content analysis on one individual and find how many of the 
  # potential preys (20 species) has been eaten by this species.
  # Remember we are in a perfect world and you detect all links with one sample!  
(z <- rbinom(nPrey, 1, psi))
  # If detection probability is 1, we can estimate link prob
  # using logistic regression
glm1 <- glm(z ~ BodySize, family=binomial)
summary(glm1)
plot(BodySize, z, xlab="BodySize", ylab="Link probability")
glm1.est <- coef(glm1)
plot(function(x) plogis(-3 + 2*x), 1, 3, add=TRUE, lwd=2)
plot(function(x) plogis(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,
     lwd=2, col="blue")
legend(1, 0.9, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=2)

#Second, we go back to the real world, where some links will be undetected

  # What if detection probability is less than one, say p=0.3?
  # Now, z can no longer be observed, ie it is latent
  # let's simulate we do gut content on 30 individuals
nIndividuals <- 30
p <- 0.3
y <- matrix(NA, nPrey, nIndividuals)
for(i in 1:nPrey) {
  y[i,] <- rbinom(nIndividuals, 1, z[i]*p)
}
  # let's see the accumulation curve of links per species sampled
cum <- c()
cum[1] <- sum(y[,1])
for(i in 2:nIndividuals){
    cum[i] <- sum(ifelse(rowSums(y[,1:i]) > 0, 1,0)) 
    }
scatter.smooth(c(1:nIndividuals), cum)
  # But gut content is expensive, and we can only afford 4 individuals
  # If so, some links will get undetected... 
  # here are the first 4 individuals
(y1 <- y[,c(1:4)])

# THE HIERARCHICAL MODELING FRAMEWORK
  # We use package unmarked. the idea is that first we estimate detection pobability 
  # and use this information to estimate the links as function of covariables.

  #First format for unmarked and summarize
#install.packages(unmarked)
library(unmarked)
umf <- unmarkedFrameOccu(y=y1, siteCovs=as.data.frame(BodySize))
summary(umf)

  # Second we Fit a model. Note: estimates are on logit scale
  # Detection covariates follow first tilde, then come prey covars
(fm.occu1 <- occu(formula= ~1 ~BodySize, data=umf)) #note several models can be AIC-compared

  # Third. We Analyze results
  # Explore the detection prob, we can back-transform using:
(beta1 <- coef(fm.occu1))
backTransform(fm.occu1, type="det") # detection prob estimate with SE

  # And undertand linkage rules (estimates of covariates)
plot(function(x) plogis(beta1[1] + beta1[2]*x), 1, 3,
     xlab="BodySize", ylab="Expected link probability", col = "blue", lwd=2)
plot(function(x) plogis(-3 + 2*x), 1, 3, add=TRUE, lwd=2) #expected data

  # And finnally we can do predictions for the actual values of BodySize, 
  # (or even other values!)
newdat <- data.frame(BodySize=BodySize)
temp <- predict(fm.occu1, type="state", newdata=newdat)
  # And we compare all data visually
data <- as.data.frame(cbind(y1, BodySize, rowSums(y1), z, temp))
data

# If you think this may be useful, please use that approach.
# I am happy to colaborate, and otherwise "here is my paper, cite me maybe"
# Bartomeus I. (2013) Understanding linkage rules in plant-pollinator networks by using 
# hierarchical models that incorporate pollinator detectability and plant traits. 
# Sumitted to PLoSONE.


