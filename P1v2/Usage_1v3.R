st = sim_phyl()
plot(st$newick)
st2 = drop.fossil(st$newick)
plot(st2)
st2 = phylo2p(st2)
# Recontruction and estimations
rec = rec_
