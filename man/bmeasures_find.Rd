\name{bmeasures_find}
\alias{bmeasures_find}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding an appropriate equation using ROC analysis
}
\description{
To select the most suitable similarity and dissimilarity equations using ROC analysis
}
\usage{
bmeasures_find(inFile, setSeed=0, numSample=20, numRange=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inFile}{
      Input file ("filename.csv""). It is a two-dimensional matrix with rows represent samples and columns represent features/variables. The first row indicates the column title and the last column indicates the class of samples.
}
  \item{setSeed}{
      Set the seed of R's random number generator, which is useful for creating simulations or random objects that can be reproduced. It can be utilized by assigning setSeed=1. The default value is 0.
}
  \item{setSample}{
      The number of generated samples used to calculate the mean of area under the ROC curve (AUC).
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
This function requires ROCR package (https://rocr.bioinf.mpi-sb.mpg.de).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Finding an appropriate similarity/dissimilarity equation
# The result is also saved as "bmeasures_filename.csv"

output <- bmeasures_find(inFile="immune_cancer.csv")

# set the setSeed parameter (setSeed=1) to obtain reproducible results
output <- bmeasures_find(inFile="immune_cancer.csv", setSeed=1)

# define the number of AUC iterations used for calculating meanAUC
output <- bmeasures_find(inFile="immune_cancer.csv", setSeed=1, numSample=50)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
