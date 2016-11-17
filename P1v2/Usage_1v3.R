st = sim_phyl(ct=2,seed=8)
plot(st$newick)
st2 = drop.fossil(st$newick)
plot(st2)
st2 = phylo2p(st2)
# Recontruction and estimations
rec = rec_
