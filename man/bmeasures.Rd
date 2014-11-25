\name{bmeasures}
\alias{bmeasures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Binary Similarity/Dissimilarity Coefficient
}
\description{
To calculate the binary similarity or dissimilarity (distance) coefficient between two binary vectors
}
\usage{
bmeasures(x, y, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     The first binary vector (represented as matrix and row vector)
}
  \item{y}{
     The second binary vector (represented as matrix and row vector)
}
  \item{method}{
     The equation that can be used to calculate the binary similarity/dissimilarity coefficient. The List of binary
     similarity/dissimilarity measures are as follows: 
     \tabular{ll}{
     \code{"eq_01"} Jaccard similarity,                 \tab \code{"eq_41"} Kulczynski-2 similarity,\cr 
     \code{"eq_02"} Dice similarity,                    \tab \code{"eq_42"} Driver & Kroeber similarity,\cr          
     \code{"eq_03"} Czekanowski similarity,             \tab \code{"eq_43"} Johnson similarity,\cr
     \code{"eq_04"} 3W-Jaccard similarity,              \tab \code{"eq_44"} Dennis similarity,\cr
     \code{"eq_05"} Nei & Li similarity,                \tab \code{"eq_45"} Simpson similarity,\cr
     \code{"eq_06"} Sokal & Sneath-1 similarity,        \tab \code{"eq_46"} Braun & Banquet similarity,\cr
     \code{"eq_07"} Sokal & Michener similarity,        \tab \code{"eq_47"} Fager & McGowan similarity,\cr
     \code{"eq_08"} Sokal & Sneath-2 similarity,        \tab \code{"eq_48"} Forbes-2 similarity,\cr
     \code{"eq_09"} Roger & Tanimoto similarity,        \tab \code{"eq_49"} Sokal & Sneath-4 similarity,\cr
     \code{"eq_10"} Faith similarity,                   \tab \code{"eq_50"} Gower similarity,\cr
     \code{"eq_11"} Gower & Legendre similarity,        \tab \code{"eq_51"} Pearson-1 similarity,\cr
     \code{"eq_12"} Intersection similarity,            \tab \code{"eq_52"} Pearson-2 similarity,\cr
     \code{"eq_13"} Inner product similarity,           \tab \code{"eq_53"} Pearson-3 similarity,\cr
     \code{"eq_14"} Russell & Rao similarity,           \tab \code{"eq_54"} Pearson & Heron-1 similarity,\cr
     \code{"eq_15"} Hamming distance,                   \tab \code{"eq_55"} Pearson & Heron-2 similarity,\cr
     \code{"eq_16"} Euclid distance,                    \tab \code{"eq_56"} Sokal & Sneath-3 similarity,\cr
     \code{"eq_17"} Squared-euclid distance,            \tab \code{"eq_57"} Sokal & Sneath-5 similarity,\cr 
     \code{"eq_18"} Canberra distance,                  \tab \code{"eq_58"} Cole similarity,\cr
     \code{"eq_19"} Manhattan distance,                 \tab \code{"eq_59"} Stiles similarity,\cr
     \code{"eq_20"} Mean-Manhattan distance,            \tab \code{"eq_60"} Ochiai-2 similarity,\cr
     \code{"eq_21"} Cityblock distance,                 \tab \code{"eq_61"} Yuleq similarity,\cr
     \code{"eq_22"} Minkowski distance,                 \tab \code{"eq_62"} Yuleq distance,\cr
     \code{"eq_23"} Vari distance,                      \tab \code{"eq_63"} Yulew similarity,\cr
     \code{"eq_24"} Size Difference distance,           \tab \code{"eq_64"} Kulczynski-1 similarity,\cr
     \code{"eq_25"} Shape Difference distance,          \tab \code{"eq_65"} Tanimoto similarity,\cr 
     \code{"eq_26"} Pattern Difference distance,        \tab \code{"eq_66"} Disperson similarity,\cr
     \code{"eq_27"} Lance & Williams distance,          \tab \code{"eq_67"} Hamann similarity,\cr
     \code{"eq_28"} Bray & Curtis distance,             \tab \code{"eq_68"} Michael similarity,\cr
     \code{"eq_29"} Hellinger distance,                 \tab \code{"eq_69"} Goodman&Kruskal similarity,\cr
     \code{"eq_30"} Chord distance,                     \tab \code{"eq_70"} Anderberg similarity,\cr
     \code{"eq_31"} Cosine similarity,                  \tab \code{"eq_71"} Baroni-Urbani&Buser-1 similarity,\cr
     \code{"eq_32"} Gilbert & Wells similarity,         \tab \code{"eq_72"} Baroni-Urbani&Buser-2 similarity,\cr
     \code{"eq_33"} Ochiai-1 similarity,                \tab \code{"eq_73"} Peirce similarity,\cr
     \code{"eq_34"} Forbes-1 similarity,                 \tab \code{"eq_74"} Eyraud similarity,\cr
     \code{"eq_35"} Fossum similarity,                  \tab \code{"eq_75"} Tarantula similarity,\cr
     \code{"eq_36"} Sorgenfrei similarity,              \tab \code{"eq_76"} Ample similarity,\cr
     \code{"eq_37"} Mountford similarity,               \tab \code{"eq_77"} Derived Rusell-Rao similarity,\cr
     \code{"eq_38"} Otsuka similarity,                  \tab \code{"eq_78"} Derived Jaccard similarity,\cr
     \code{"eq_39"} Mc Connaughey similarity,           \tab \code{"eq_79"} Variant of Correlation similarity.\cr
     \code{"eq_40"} Tarwid similarity,\cr
     }
}
}
\details{
%%List of 79 binary similarity/dissimilarity measures\cr

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
b_coeff <- bmeasures(ob1, ob2, method="eq_01")


# Binary coeffients between object samples
dsample <- dim(sample)
for (data1 in 1:(dsample[1]-1))
{
  for(data2 in (data1+1):dsample[1])
  {
    obj1 <- t(as.matrix(sample[data1,]))
    obj2 <- t(as.matrix(sample[data2,]))
    
    cat(data1, " - ", data2, "\n")
    cat(obj1,"\n")
    cat(obj2,"\n")
    
    result <- bmeasures(obj1, obj2, method="eq_01")
    
    id <- paste(data1,"-",data2)
    id <- gsub(" ","",id)
    
    # result: a, b, c, d, similarity/dissimilarity coefficient
    result <- c(id=id, result)
    
    # binding the results as output
    if ((data1==1) & (data2==(data1+1)))
    {
      output <- result
    } else
    {
      output <- rbind(output, result)
    }
  }
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
