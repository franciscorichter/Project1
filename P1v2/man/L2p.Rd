\name{L2p}
\alias{L2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to transforme an L table to an p object
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
L2p(L,ct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{L table
%%     ~~Describe \code{x} here~~
}
  \item{ct}{crown time
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
s3 <- dd_sim(c(0.8,0.1,40),15)
s1 <-Ltop1(L=s3$L)
p2 <- mle_dd(n=s1$n,E=s1$E,t=s1$t)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

