\name{reconst_tree}
\alias{reconst_tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reconstructed phylogenetic tree
%%  ~~function to do ... ~~
}
\description{Given an incomplete phylogenetic tree and parameters of the model, this function is able to generate simulated reconstructed trees.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reconst_tree(bt, pars, model = "dd", tt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bt}{Branching times of the incomplete phylogenetic tree
%%     ~~Describe \code{bt} here~~
}
  \item{pars}{Parameters of the model
%%     ~~Describe \code{pars} here~~
}
  \item{model}{Model for diversification rates (only diversity-dependence model available so far)
%%     ~~Describe \code{model} here~~
}
  \item{tt}{Crown time
%%     ~~Describe \code{tt} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
s <- phyl2()  # simulate phylogenetic tree
p <- mle_dd(n=s$n,E=s$E,t=s$t) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)

p2 <- mle_dd(n=s2$n,E=s2$E,t=s2$t) # estimate parameters of incomplete tree
c(p2$lambda,(p2$lambda-p2$mu)/p2$beta,p2$mu) #we can see no nice estimation and mu=0 always

rec_tree =reconst_tree(bt=s2$t,pars=c(0.8,0.0175,0.1),tt=15) #simulate a random recostructed tree with given parameters

p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
c(p3$lambda,(p3$lambda-p3$mu)/p3$beta,p3$mu) # estimation of parameters of the reconstructed tree

}

