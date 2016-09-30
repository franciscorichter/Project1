\name{llik}
\alias{llik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Likelihood function under the diversity-depencence model
%%  ~~function to do ... ~~
}
\description{ It computes the likelihood function of the dd model...
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
llik(b,n,E,t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{ vector with parameters lambda,beta and mu.
%%     ~~Describe \code{x} here~~
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
s <- dd_sim(c(0.8,0.1,40),15)
llik(c(0.8,0.1,40),n=s$n,s$E,s$t)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

