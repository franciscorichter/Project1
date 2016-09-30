\name{phyl2}
\alias{phyl2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to generate random tree, with Newick output.
%%  ~~function to do ... ~~
}
\description{Function to generate random trees, given the parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
 phyl2()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{Crown time}
  \item{lambda,mu,etc}{parameters for the simulation}
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
s <- phyl2(tt=15)
p <- mle_dd(n=s$n,E=s$E,t=s$t)
c(p$lambda,p$beta,p$mu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
