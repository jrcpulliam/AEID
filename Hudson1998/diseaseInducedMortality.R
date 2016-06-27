## Disease-induced mortality in the Hudson model

## per capita (grouse) disease-induced mortality as a function
## infection of prevalence
diMort <- function(prev,alpha=3*10^-4,k=1){
  alpha*(1+prev*(k+1)/k)
}

par(bty='L',lwd=3)
curve(diMort(x),xlab='Prevalence',ylab='Per capita disease induced mortality rate')
curve(diMort(x,k=1.8),add=T,lty=2)
curve(diMort(x,k=0.5),add=T,lty=3)


barplot(dnbinom(0:20,1,mu=10))
barplot(dnbinom(0:20,1.8,mu=10))
barplot(dnbinom(0:20,0.5,mu=10))
