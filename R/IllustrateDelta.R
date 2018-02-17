
# This script was used to generate Figure 4 in the paper
# Morten Morup and Lars K. Hansen, "Archetypal
# Analysis for Machine Learning and Data Mining", submitted NeuroComputing
# 2011
illustrate_delta <- function() {
  N=1000     # Number of observations
  tresh=0.8  # Level of truncation of simplex
  DD=3       # Dimensionality of simplex

  # Generate synethetic data
  XC=matrix(c(cos(0), cos(2*pi/3), cos(2*pi/3*2), sin(0),sin(2*pi/3), sin(2*pi/3*2)), nrow=2, byrow=T)
  S=-log(matrix(runif(DD*N), nrow=DD))
  S=t(t(S)/colSums(S))
  IJ=which(S>tresh, arr.ind = T)
  I=IJ[,1]; J=IJ[,2]
  S <- S[,-J]
  NN=ncol(S)
  # Add noise with standard deviation sigma
  sigma=0.0
  X=XC%*%S + sigma*matrix(rnorm(2*NN), nrow=2)

  # Plot the generated data
  #figure;
  #hold on;
  plot(X[1,],X[2,])

  # Estimate the AA/PCH model using PCHA.m
  noc=3  # Number of components
  ms2=1  # Widht of line in generated plot
  delta=c(0, 0.25, 0.5) # values of \delta

  colors=c('black','red','green') # Mark in seperate colors the 3 different PCHA solutions
  for(k in 1:3) {
    ret = PCHA(X,noc,1:ncol(X),1:ncol(X),delta[k])
    lines(ret$XC[1, c(1,2,3,1)], ret$XC[2,c(1,2,3,1)],col=colors[k]) #,'linewidth',1.5*ms2,'LineStyle','-')
    #axis off;
    #axis equal;
  }
  legend('topright',
         legend=c('observations','delta=0','delta=0.25','delta=0.5'),
         col=colors, lty=1)

}
