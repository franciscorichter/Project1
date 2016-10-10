library("grid")
l=0.2 #Length of the square
grid.newpage()
gr <- rectGrob(width=l, height=l, name="gr") #Basic Square
pts <- data.frame(level=1, x=0.5, y=0.1, alfa=0) #Centers of the squares
for (i in 2:10) #10=Deep of the fractal. Feel free to change it
{
  df<-pts[pts$level==i-1,]
  for (j in 1:nrow(df))
  {
    pts <- rbind(pts,
                 c(i,
                   df[j,]$x-2*l*((1/sqrt(2))^(i-1))*sin(df[j,]$alfa+pi/4)-0.5*l*((1/sqrt(2))^(i-2))*sin(df[j,]$alfa+pi/4-3*pi/4),
                   df[j,]$y+2*l*((1/sqrt(2))^(i-1))*cos(df[j,]$alfa+pi/4)+0.5*l*((1/sqrt(2))^(i-2))*cos(df[j,]$alfa+pi/4-3*pi/4),
                   df[j,]$alfa+pi/4))
    pts <- rbind(pts,
                 c(i,
                   df[j,]$x-2*l*((1/sqrt(2))^(i-1))*sin(df[j,]$alfa-pi/4)-0.5*l*((1/sqrt(2))^(i-2))*sin(df[j,]$alfa-pi/4+3*pi/4),
                   df[j,]$y+2*l*((1/sqrt(2))^(i-1))*cos(df[j,]$alfa-pi/4)+0.5*l*((1/sqrt(2))^(i-2))*cos(df[j,]$alfa-pi/4+3*pi/4),
                   df[j,]$alfa-pi/4))
  }
}
for (i in 1:nrow(pts))
{
  grid.draw(editGrob(gr, vp=viewport(x=pts[i,]$x, y=pts[i,]$y, w=((1/sqrt(2))^(pts[i,]$level-1)), h=((1/sqrt(2))^(pts[i,]$level-1)), angle=pts[i,]$alfa*180/pi),
                     gp=gpar(col=0, lty="solid", fill=rgb(139*(nrow(pts)-i)/(nrow(pts)-1),
                                                          (186*i+69*nrow(pts)-255)/(nrow(pts)-1),
                                                          19*(nrow(pts)-i)/(nrow(pts)-1),
                                                          alpha= (-110*i+200*nrow(pts)-90)/(nrow(pts)-1), max=255))))
}


library(ggplot2)
library(numDeriv)
library(RColorBrewer)
library(gridExtra)
## Polynom: choose only one or try yourself
f  <- function (z) {z^3-1}        #Blurry 1
#f  <- function (z) {z^4+z-1}     #Blurry 2
f  <- function (z) {z^5+z^3+z-1} #Blurry 3
z <- outer(seq(-2, 2, by = 0.01),1i*seq(-2, 2, by = 0.01),'+')
for (k in 1:5) z <- z-f(z)/matrix(grad(f, z), nrow=nrow(z))
## Supressing texts, titles, ticks, background and legend.
opt <- theme(legend.position="none",
             panel.background = element_blank(),
             axis.ticks=element_blank(),
             axis.title=element_blank(),
             axis.text =element_blank())
z <- data.frame(expand.grid(x=seq(ncol(z)), y=seq(nrow(z))), z=as.vector(exp(-Mod(f(z)))))
# Create plots. Choose a palette with display.brewer.all()
p1 <- ggplot(z, aes(x=x, y=y, color=z)) + geom_tile() + scale_colour_gradientn(colours=brewer.pal(8, "Paired")) + opt
p2 <- ggplot(z, aes(x=x, y=y, color=z)) + geom_tile() + scale_colour_gradientn(colours=brewer.pal(7, "Paired")) + opt
p3 <- ggplot(z, aes(x=x, y=y, color=z)) + geom_tile() + scale_colour_gradientn(colours=brewer.pal(6, "Paired")) + opt
p4 <- ggplot(z, aes(x=x, y=y, color=z)) + geom_tile() + scale_colour_gradientn(colours=brewer.pal(5, "Paired")) + opt
# Arrange four plots in a 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol=2)


depth <- 9
angle<-30 #Between branches division
L <- 0.90 #Decreasing rate of branches by depth
nstars <- 300 #Number of stars to draw
mstars <- matrix(runif(2*nstars), ncol=2)
branches <- rbind(c(1,0,0,abs(jitter(0)),1,jitter(5, amount = 5)), data.frame())
colnames(branches) <- c("depth", "x1", "y1", "x2", "y2", "inertia")
for(i in 1:depth)
{
  df <- branches[branches$depth==i,]
  for(j in 1:nrow(df))
  {
    branches <- rbind(branches, c(df[j,1]+1, df[j,4], df[j,5], df[j,4]+L^(2*i+1)*sin(pi*(df[j,6]+angle)/180), df[j,5]+L^(2*i+1)*cos(pi*(df[j,6]+angle)/180), df[j,6]+angle+jitter(10, amount = 8)))
    branches <- rbind(branches, c(df[j,1]+1, df[j,4], df[j,5], df[j,4]+L^(2*i+1)*sin(pi*(df[j,6]-angle)/180), df[j,5]+L^(2*i+1)*cos(pi*(df[j,6]-angle)/180), df[j,6]-angle+jitter(10, amount = 8)))
  }
}
nodes <- rbind(as.matrix(branches[,2:3]), as.matrix(branches[,4:5]))
png("image.png", width = 1200, height = 600)
plot.new()
par(mai = rep(0, 4), bg = "gray12")
plot(nodes, type="n", xlim=c(-7, 3), ylim=c(0, 5))
for (i in 1:nrow(mstars))
{
  points(x=10*mstars[i,1]-7, y=5*mstars[i,2], col = "blue4", cex=.7, pch=16)
  points(x=10*mstars[i,1]-7, y=5*mstars[i,2], col = "blue",  cex=.3, pch=16)
  points(x=10*mstars[i,1]-7, y=5*mstars[i,2], col = "white", cex=.1, pch=16)
}
# The moon
points(x=-5, y=3.5, cex=40, pch=16, col="lightyellow")
# The tree
for (i in 1:nrow(branches)) {lines(x=branches[i,c(2,4)], y=branches[i,c(3,5)], col = paste("gray", as.character(sample(seq(from=50, to=round(50+5*branches[i,1]), by=1), 1)), sep = ""), lwd=(65/(1+3*branches[i,1])))}
rm(branches)
dev.off()


u=0.918 #Parameter between 0 and 1
n=2000 #Number of particles
m=40 #Number of iterations
ikeda=data.frame(it=1,x1=runif(n, min = -40, max = 40), y1=runif(n, min = -40, max = 40))
ikeda$x2=1+u*(ikeda$x1*cos(0.4-6/(1+ikeda$x1^2+ikeda$y1^2))-ikeda$y1*sin(0.4-6/(1+ikeda$x1^2+ikeda$y1^2)))
ikeda$y2=  u*(ikeda$x1*sin(0.4-6/(1+ikeda$x1^2+ikeda$y1^2))+ikeda$y1*cos(0.4-6/(1+ikeda$x1^2+ikeda$y1^2)))
for (k in 1:m)
{
  df=as.data.frame(cbind(rep(k+1,n),
                         ikeda[ikeda$it==k,]$x2,
                         ikeda[ikeda$it==k,]$y2,
                         1+u*(ikeda[ikeda$it==k,]$x2*cos(0.4-6/(1+ikeda[ikeda$it==k,]$x2^2+ikeda[ikeda$it==k,]$y2^2))-ikeda[ikeda$it==k,]$y2*sin(0.4-6/(1+ikeda[ikeda$it==k,]$x2^2+ikeda[ikeda$it==k,]$y2^2))),
                         u*(ikeda[ikeda$it==k,]$x2*sin(0.4-6/(1+ikeda[ikeda$it==k,]$x2^2+ikeda[ikeda$it==k,]$y2^2))+ikeda[ikeda$it==k,]$y2*cos(0.4-6/(1+ikeda[ikeda$it==k,]$x2^2+ikeda[ikeda$it==k,]$y2^2)))))
  names(df)=names(ikeda)
  ikeda=rbind(df, ikeda)
}
plot.new()
par(mai = rep(0, 4), bg = "gray12")
plot(c(0,0),type="n", xlim=c(-35, 35), ylim=c(-35,35))
apply(ikeda, 1, function(x) lines(x=c(x[2],x[4]), y=c(x[3],x[5]), col = paste("gray", as.character(min(round(jitter(x[1]*80/(m-1)+(20*m-100)/(m-1), amount=5)), 100)), sep = ""), lwd=0.1))
