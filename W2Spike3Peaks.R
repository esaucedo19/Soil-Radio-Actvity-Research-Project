library(plotrix)
Q2Spike <-read.table("oxford-w2-soil-spike-stop-5-4-2023.txt",header=FALSE)
x<-c(1:8192)
plotCI(x,Q2Spike[,1],sqrt(Q2Spike[,1]),col = "purple",xlim = c(5800,6800),
       ylim = c(300,4500),main = "Oxford W2 Spiked Soil",xlab = "Channel", ylab = "Counts")
df1 <- data.frame(x=x[5800:6400],y=Q2Spike[5800:6400,1])

fitsSpike1 <- nls(y~a*exp(-(x-c)^2/(2*b^2))+d*x+e,df1,start = list(a=5500,b=100,c=6270,d =-0.14,e = 250),
                    trace = T, weights = 1/df1[,2])

df2 <- data.frame(x=x[6400:6900],y=Q2Spike[6400:6900,1])

fitsSpike2 <- nls(y~f*exp(-(x-h)^2/(2*g^2))+j*exp(-(x-l)^2/(2*k^2))+d*x+e,df2
                 ,start = list(d =-0.14,e = 250,f=800,g=80,h=6550,j=800,k=80,l=6725),
                 trace = T, weights = 1/df2[,2])

fitsSpike1;summary(fitsSpike1)
fitsSpike2;summary(fitsSpike2)
d=-0.43896;e=3344.38993;f=223.44521;g=72.49046;h=6566.65984
j=257.40546;k=56.44585;l=6742.15508 

#Gaussian Graphed Separate
curve(f*exp(-(x-h)^2/(2*g^2))+j*exp(-(x-l)^2/(2*k^2))+d*x+e,from = 6400, to = 6900, 
      col = "pink",add = T, lwd = 5)

a=3513.07807;b=87.92918;c=6179.37917;d=-0.41565;e=2977.74193
curve(a*exp(-(x-c)^2/(2*b^2))+d*x+e,from = 5800, to = 6400, 
      col = "pink",add = T, lwd = 5)


#Graphed Together
df3 <- data.frame(x=x[5900:6900],y=Q2Spike[5900:6900,1])
fitsSpikeAll <- nls(y~a*exp(-(x-c)^2/(2*b^2))+f*exp(-(x-h)^2/(2*g^2))+j*exp(-(x-l)^2/(2*k^2))+d*x+e
                    ,df3,start = list(a=5500,b=81,c=6275,d =-0.045,e =394,f=800,
                                      g=82,h=6550,j=800,k=80,l=6750), trace = T,
                    weights = 1/df3[,2])

fitsSpikeAll;summary(fitsSpikeAll)
a=3555.4132;b=86.4410;c=6181.1168;d=-0.7499
e=5033.0913;f=293.9218;g=72.3330;h=6507.8125
j=638.0083;k=160.7200;l=6742.6250

curve(a*exp(-(x-c)^2/(2*b^2))+f*exp(-(x-h)^2/(2*g^2))+j*exp(-(x-l)^2/(2*k^2))+d*x+e,
      from = 5900, to = 6900, 
      col = "red",add = T, lwd = 5)

nls.control(printEval = T)
