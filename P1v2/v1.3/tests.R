#Tests

##1. test to check update_tree
st = sim_phyl(ct=2,seed=8)
plot(st$newick)
st$t
st$br
#op1
t_spe = 0.8
t_ext = 1.1
update_tree(wt=st$t,t_spe=t_spe,t_ext = t_ext,E=st$E,n=st$n) # works properly
#op2
t_spe = 0.8
t_ext = 0.85
update_tree(wt=st$t,t_spe=t_spe,t_ext = t_ext,E=st$E,n=st$n) # works properly
#op3
t_spe = 1.1
t_ext = 1.9
update_tree(wt=st$t,t_spe=t_spe,t_ext = t_ext,E=st$E,n=st$n) # works properly
#op4
t_spe = 1.1
t_ext = 1.2
update_tree(wt=st$t,t_spe=t_spe,t_ext = t_ext,E=st$E,n=st$n) # works properly
#op5
t_spe = 1.9
t_ext = 1.95
update_tree(wt=st$t,t_spe=t_spe,t_ext = t_ext,E=st$E,n=st$n) # works properly

#2. test to check phylo2p
