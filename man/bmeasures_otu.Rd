\name{bmeasures_otu}
\alias{bmeasures_otu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The quantities of the OTUs table
}
\description{
To generate the quantities of the Operational Taxonomic Units (OTU) table
}
\usage{
bmeasures_otu(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     The binary row vector which is represented as a row matrix
}
  \item{y}{
     The binary row vector which is represented as a row matrix
}
}

\details{
The expression of two binary instances x and y can be expressed by the OTUs table in 2x2 contingency table,
where:
\code{a} is the number of features where the value of both x and y are 1 (positive matches),
\code{b} and \code{c} are the number of features where the value of x is 0 and y is 1 and vice versa, respectively (absence mismatches), and \code{d} is the number of features where the values of both x and y are 0 (negative matches).

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
# Generates sample data
row <- 30
col <- 10
sample <- matrix(sample(0:1, row * col, replace = TRUE), row, col)

ob1 <- t(as.matrix(sample[1,]))
ob2 <- t(as.matrix(sample[2,]))

# Defined the quantities of the OTUs table
otu <- bmeasures_otu(ob1, ob2)

# Show the result
otu
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
