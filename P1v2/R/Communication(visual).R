#### newick trees for communication
#simplest tree
obs = '((B:1,A:1):1,C:2);'
#library(ape)
obs = read.tree(text=obs)
plot(obs)
# complete/incomplete trees
comp = '((B:1,(A:0.6,D:0.3):0.4):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
comp= read.tree(text = comp)
plot(comp)
#library(phytools)
et <- fancyTree(comp, type = "droptip", tip = getExtinct(comp), cex = 0.7)
## reconstruction process
par(mfrow=c(4,1))
rec0 = '((B:1,A:1):1,C:2);'
rec0 = read.tree(text = rec0)
plot(rec0)
rec1 = '((B:1,A:1):1,(F:1,C:1.5):0.5);'
rec1 = read.tree(text = rec1)
plot(rec1)
rec2 = '((B:1,A:1):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
rec2 = read.tree(text = rec2)
plot(rec2)
rec3 = '((B:1,(A:0.6,D:0.3):0.4):1,((E:0.6,F:0.8):0.2,C:1.5):0.5);'
rec3 = read.tree(text = rec3)
plot(rec3)

