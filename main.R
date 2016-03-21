#Given a function, calculates the
#Fourier series up to n terms

#Clear previous variables
rm(list=ls())

#For complex integration
library(elliptic)

#The function to approximate
fun = function(x){
    return(exp(-abs(x)))
}

#Upper bound (lower bound is -bound)
bound = pi

#a_0
a_0 = 1/(2*bound)*integrate(Vectorize(fun),-bound,bound)$value
#Half-interval
#a_0 = 1/(bound)*integrate(Vectorize(fun),0,bound)$value

#The number of repitions
n = 30 

#Vectors for the coefficients
a = rep(0,n)
b = rep(0,n)
c = rep(0,n)

#Calculate the coefficients
for(k in 1:n){
    a[k]=1/bound * integrate(function(x){return(fun(x)*cos(k*pi*x/bound))},-bound,bound)$value
    b[k]=1/bound * integrate(function(x){return(fun(x)*sin(k*pi*x/bound))},-bound,bound)$value
    
    #For half interval
    #a[k]=2/bound * integrate(function(x){return(fun(x)*cos(k*pi*x/bound))},0,bound)$value
    #b[k]=2/bound * integrate(function(x){return(fun(x)*sin(k*pi*x/bound))},0,bound)$value
    
    #Complex coefficients
    c[k] = 1/(2*bound) * myintegrate(function(x){return(fun(x)*exp(complex(imaginary=-k*x)))},-bound,bound)
}

#Fourier expansion
s_n = function(x){
    s = a_0
    for(i in 1:n)
        s=s+a[i]*cos(i*pi*x/bound)+b[i]*sin(i*pi*x/bound)
    return(s)
}

# http://stackoverflow.com/questions/7144118/how-to-save-a-plot-as-image-on-the-disk

#x-coordinates
mult = 2 #How wide should the plot be?
xc = seq(-mult*bound,mult*bound,by=0.01)

#plot fun
plot(xc,lapply(xc,fun),type="l",xlab="x",ylab="y",main="f(x)")
grid()
axis(1, pos=0, labels=F)
axis(2, pos=0, labels=F)
# https://stat.ethz.ch/pipermail/r-help/2007-September/141870.html
