library(plotrix)

#Background
Q2Bck <-read.table("oxford-soil-stop-4-27-2023.txt",header=FALSE)
x<-c(1:8192)
plotCI(x,Q2Bck[,1],sqrt(Q2Bck[,1]),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Bck[,1],sqrt(Q2Bck[,1]),xlim = c(6000,6350),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Bck[,1],sqrt(Q2Bck[,1]),col = "purple",xlim = c(5900,6400),ylim = c(100,600),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
df <- data.frame(x=x[5900:6400],y=Q2Bck[5900:6400,1])
fitsQ2Bck <- nls(y~a*exp(-(x-c)^2/(2*b^2))+d*x+e,df,start = list(a=550,b=100,c=6215,d =-0.14,e = 250),trace = T, weights = 1/df[,2])
fitsQ2Bck;summary(fitsQ2Bck)
a=2.980e+02 ;b=8.497e+01;c=6.187e+03;d=-2.663e-02;e=2.960e+02  
plotCI(x,Q2Bck[,1],sqrt(Q2Bck[,1]),col = "purple",xlim = c(5900,6400),ylim = c(100,600),main = "Oxford W2 Raw Soil Counts",xlab = "Channel", ylab = "Counts")
curve(a*exp(-(x-c)^2/(2*b^2))+d*x+e,from = 5900, to = 6400, col = "pink",add = T, lwd = 5)
curve(d*x+e,from = 5900, to = 6400, col = "red",add = T, lwd = 5)

#RawSoil
Q2Raw <-read.table("oxford-w2-soil-spike-stop-5-4-2023.txt",header=FALSE)
x<-c(1:8192)
plotCI(x,Q2Raw[,1],sqrt(Q2Raw[,1]),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Raw[,1],sqrt(Q2Raw[,1]),xlim = c(6000,6350),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Raw[,1],sqrt(Q2Raw[,1]),col = "purple",xlim = c(5900,6400),ylim = c(1000,5000),main = "Oxford Background Counts",xlab = "Channel", ylab = "Counts")
df <- data.frame(x=x[5900:6400],y=Q2Raw[5900:6400,1])
fitsQ2Raw <- nls(y~a*exp(-(x-c)^2/(2*b^2))+d*x+e,df,start = list(a=5000,b=100,c=6175,d =-0.14,e = 250),trace = T, weights = 1/df[,2])
fitsQ2Raw;summary(fitsQ2Bck)
a=2.980e+02 ;b=8.497e+01;c=6.187e+03;d=-2.663e-02;e=2.960e+02  
plotCI(x,Q2Raw[,1],sqrt(Q2Raw[,1]),col = "red",xlim = c(5900,6400),ylim = c(1000,5000),main = "Oxford W2 Raw Soil Counts",xlab = "Channel", ylab = "Counts")
curve(a*exp(-(x-c)^2/(2*b^2))+d*x+e,from = 5900, to = 6400, col = "blue",add = T, lwd = 5)
curve(d*x+e,from = 5900, to = 6400, col = "black",add = T, lwd = 5)

#SpikedSoil
Q2Spike <-read.table("oxford-w2-soil-spike-stop-5-4-2023.txt",header=FALSE)
x<-c(1:8192)
plotCI(x,Q2Spike[,1],sqrt(Q2Spike[,1]),main = "Oxford Spiked Soil",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Spike[,1],sqrt(Q2Spike[,1]),xlim = c(6000,6350),main = "Oxford Spiked Soil",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Spike[,1],sqrt(Q2Spike[,1]),col = "purple",xlim = c(5800,6800),
       ylim = c(300,4500),main = "Oxford W2 Spiked Soil",xlab = "Channel", ylab = "Counts")
plotCI(x,Q2Spike[,1],sqrt(Q2Spike[,1]),col = "purple",xlim = c(5900,6400),
       ylim = c(500,5000),main = "Oxford W2 Spiked Soil",xlab = "Channel", ylab = "Counts")
df <- data.frame(x=x[5800:6800],y=Q2Spike[5800:6800,1])

fitsQ2Spike <- nls(y~a*exp(-(x-c)^2/(2*b^2))+d*x+e+f*exp(-(x-h)^2/(2*g^2))
                   +(i*exp(-(x-k)^2/(2*j^2))),df,
                   start = list(a=4000,b=90,c=6215,d =-0.14,e = 250,f=800,g=90,
                               h=6555,i=800,j=6650,k=90),trace = T, weights = 1/df[,2])

fitsQ2Spike;summary(fitsQ2Spike)
fitsQ2S2 <- nls(y~i*exp(-(x-j)^2/(2*k^2)),df,
                start = list(i=800,j=6650,k=90),trace = T,weights = 1/df[,2])
fitsQ2Spike;summary(fitsQ2S2)
a=3528.76702 ;b=86.94500;c=6178.74116;d=-0.42153;e=3014.92323;f=459.68271;
g=100.74981;h=6572.98139
i=3079.140;j=6172.168;k= 133.073
curve(a*exp(-(x-c)^2/(2*b^2))+d*x+e+f*exp(-(x-h)^2/(2*g^2)),from = 5800, to = 6800, 
      col = "pink",add = T, lwd = 5)
curve(i*exp(-(x-j)^2/(2*k^2)),from = 6700, to = 6800, col = "pink",add = T,lwd = 5)
curve(d*x+e,from = 5800, to = 6800, col = "black",add = T, lwd = 5)

