#### newick trees for communication
#simplest tree
obs = '((B:1,A:1):1,C:2);'
#library(ape)
obs = read.tree(text=obs)
plot(obs,show.tip.label = F)
# complete/incomplete trees
comp = '((B:1,(A:0.6,D:0.3):0.4):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
comp= read.tree(text = comp)
plot(comp,show.tip.label = F)
#library(phytools)
et <- fancyTree(comp, type = "droptip", tip = getExtinct(comp), cex = 0.7)
## reconstruction process
par(mfrow=c(4,1))
rec0 = '((B:1,A:1):1,C:2);'
rec0 = read.tree(text = rec0)
plot(rec0,show.tip.label = FALSE,type='cladogram')
rec1 = '((B:1,A:1):1,(F:1,C:1.5):0.5);'
rec1 = read.tree(text = rec1)
plot(rec1,show.tip.label = F,type = 'cladogram')
plot(rec1,show.tip.label = F,edge.color = c("black","black","black","black","darkgreen","black"),edge.width = 2,edge.lty = c(rep(1,4),4,1))
rec2 = '((B:1,A:1):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
rec2 = read.tree(text = rec2)
plot(rec2)
rec3 = '((B:1,(A:0.6,D:0.3):0.4):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
rec3 = read.tree(text = rec3)
plot(rec3)

## code from stackoverflow
library(jpeg)
logo <- readJPEG("Downloads/Symbol1.jpg")
logo2 <- as.raster(logo)
r <- nrow(logo2)/ncol(logo2) # aspect ratio
s <- 0.4 # symbol size

# display plot to obtain its size
plot(rec1, edge.color = c("black","black","black","black","darkgreen","black"),edge.width = 2, edge.lty = c(rep(1,4),4,1))
lims <- par("usr") # plot area size
file_r <- (lims[2]-lims[1]) / (lims[4]-lims[3]) # aspect ratio for the file
file_s <- 480   # file size

# save tree with added symbol
png("tree_logo.png", height=file_s, width=file_s*file_r)
plot(rec1, show.tip.label = F,
     edge.color = c("black","black","black","black","darkgreen","black"),
     edge.width = 2, edge.lty = c(rep(1,4),4,1))
rasterImage(logo2, 1.6, 2.8, 1.6+s/r, 2.8+s)

# add axis
axisPhylo()
mtext(expression(Delta*italic("t")["i"]), side = 1, line = 3)
dev.off()

