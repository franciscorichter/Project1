library(P1)
s = phyl2(seed=runif(1,1,100000))
par(mfrow=(c(1,2))) # 1 row but 2 plot panels
plot(s$newick,direction = 'upwards',show.tip.label=FALSE,edge.width=2,edge.color = 'purple',type="fan")
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex,direction = 'upwards',show.tip.label=FALSE,edge.width=2,edge.color = 'purple',type="fan")
axisPhylo()


M=data.frame()
names(M) = c('lambda','beta','mu')



#WIP

B = cbind(BBB,llik_obt,llik_real)
data.viewer()
deducer()


itt=100
tt=5
Vals = vector(mode='list',length = itt)
for (k in 1:itt){
  #print(k)
  s1 = phyl2(tt=tt, lambda0=0.82, mu0=0.47, K=Inf, seed=runif(1,1,10000))
  Vals[[k]]=data.frame(time=cumsum(s1$t),n=s1$n)
}

#plot(approx(Vals[[6]]$time,Vals[[6]]$n))
#approx(Vals[[6]]$time,Vals[[6]]$n,xout = seq(0,15,by=0.1))

AA = matrix(nrow=length(seq(0,tt,by=0.1)),ncol=itt)
for (i in 1:itt){
  AA[,i]=approx(Vals[[i]]$time,Vals[[i]]$n,xout = seq(0,tt,by=0.1))$y
}

plot((seq(0,tt,by=0.1)-tt),rowMeans(AA, na.rm = TRUE, dims = 1))


vils = data.frame(time=cumsum(bt),n=c(2:length(bt),length(bt)))
A2 = approx(vils$time,vils$n,xout = seq(0,tt,by=0.1))$y
A3 = approx(dendroica, 2:(length(dendroica)+1),xout = seq(0,5,by=0.1))$y
a1 = data.frame(time=(seq(0,5,by=0.1)-5),val=rowMeans(AA, na.rm = TRUE, dims = 1))
a2 = data.frame(time=(seq(0,5,by=0.1)-5),val=A3,mu ='reconstructed phylogeny')
smp = rbind(a1,a2)
ggplot(smp,aes(x=time, y=val))+geom_smooth()+ylab('number of lineages (log)')+xlab('time (Myr)')+ scale_y_log10()+ggtitle('Dendroica')
plot(s$newick)

data.frame(time=(seq(0,5,by=0.1)-5),val=rowMeans(AA, na.rm = TRUE, dims = 1),mu='complete phylogeny')

itt=10000
Vals = vector(mode='list',length = itt)
for (k in 1:itt){
  #print(k)
  s1 = phyl2(tt=64.95, lambda0=0.16, mu0=0.03, K=45, seed=runif(1,1,10000))
  Vals[[k]]=data.frame(time=cumsum(s1$t),n=s1$n)
}

#plot(approx(Vals[[6]]$time,Vals[[6]]$n))
#approx(Vals[[6]]$time,Vals[[6]]$n,xout = seq(0,15,by=0.1))

AA = matrix(nrow=length(seq(0,5,by=0.1)),ncol=itt)
for (i in 1:itt){
  AA[,i]=approx(Vals[[i]]$time,Vals[[i]]$n,xout = seq(0,5,by=0.1))$y
}

plot((seq(0,5,by=0.1)-5),rowMeans(AA, na.rm = TRUE, dims = 1),log="y")

vils = data.frame(time=cumsum(bt),n=c(2:length(bt),length(bt)))
A2 = approx(vils$time,vils$n,xout = seq(0,tt,by=0.1))$y

smp = cbind(data.frame(time=(seq(0,tt,by=0.1)-tt),val=rowMeans(AA, na.rm = TRUE, dims = 1),mu='complete phylogeny'))
ggplot(smp,aes(x=time, y=val))+geom_smooth()+ylab('number of lineages (log)')+xlab('time (Myr)')+ scale_y_log10()+ggtitle('Foraminifera')



nmnight = mnight
nmnight[,2] = (nmnight[,1]-nmnight[,3])/(nmnight[,2])

xtable(head(rbind(c(4,30,1),nmnight),n=7))

nm = M
nm[,2] = (nm[,1]-nm[,3])/(nm[,2])

xtable(head(rbind(c(4,30,1),nm),n=9))
