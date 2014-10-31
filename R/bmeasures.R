#' The OTUs table  
#' 
#' Generates the OTUs table between 2 binary vectors
#' @param x The first binary vector (represented as matrix and row vector)
#' @param y The second binary vector (represented as matrix and row vector)
#' @return The four quantities a,b,c,d in the OTUs table are defined as follows: a is the number of features where the value of both x and y are 1 (positive matches), b and c are the number of features where the value of x is 0 and y is 1 and vice versa, respectively (absence mismatches) and d is the number of features where the values of both x and y are 0 (negative matches).
#' @export 
bmeasures_otu <- function (x,y)
{
  dx <- dim(x)
  dy <- dim(y)
  
  if (dx[1] != dy[1])
  {
    cat("Inputs have different number of rows \n")
  } else if (dx[2] != dy[2])
  {
    cat("Inputs have different number of columns \n")
  } else
  {
    cat("")
  }
  
  # bind row vectors
  input <- as.matrix(rbind(x,y))
  
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  for (i in 1:dx[2])
  {
    k <- input[1,i]
    l <- input[2,i]
    
    if ((k==1) & (l==1))
    {
      a <- a + 1
    } else if ((k==0) & (l==1))
    {
      b <- b + 1
    } else if ((k==1) & (l==0))
    {
      c <- c + 1
    } else # ((k==0) & (l==0))
    {
      d <- d + 1
    }
  }
  
  otu_out <- c(a=a,b=b,c=c,d=d)
  remove (a,b,c,d)
  return (otu_out)
}


#' The Binary Similarity/Dissimilarity Coefficient  
#' 
#' Calculate the binary similarity/dissimilarity coefficient between 2 binary vectors
#' @param x The first binary vector (represented as matrix and row vector)
#' @param y The second binary vector (represented as matrix and row vector)
#' @param method The equation ID e.g. eq_01
#' @return The binary similarity/dissimilarity coeffifient
#' @export 
bmeasures <- function (x,y,method)
{
  otuTable <- bmeasures_otu(x,y)
  a <- otuTable[1]
  b <- otuTable[2]
  c <- otuTable[3]
  d <- otuTable[4]
  n <- a+b+c+d;
  
  switch (method,
  eq_01 = {
#     cat("equation 1 \n")
    coef <- a/(a+b+c)
  },
  eq_02 = {
#     cat("equation 2 \n")
    coef <- a/(2*a+b+c)
  },
  eq_03 = {
#     cat("equation 3 \n")
    coef <- (2*a)/(2*a+b+c)
  },
  eq_04 = {
#     cat("equation 4 \n")
    coef <- (3*a)/(3*a+b+c)
  },
  eq_05 = {
#     cat("equation 5 \n")
    coef <- (2*a)/((a+b)+(a+c))
  },
  eq_06 = {
#     cat("equation 6 \n")
    coef <- a/(a+(2*b)+(2*c))
  },
  eq_07 = {
#     cat("equation 7 \n")
    coef <- (a+d)/n
  },
  eq_08 = {
#     cat("equation 8 \n")
    coef <- (2*(a+d))/((2*a)+b+c+(2*d))
  },
  eq_09 = {
#     cat("equation 9 \n")
    coef <- (a+d)/(a+(2*(b+c))+d)
  },
  eq_10 = {
#     cat("equation 10 \n")
    coef <- (a+(0.5*d))/n
  },
  eq_11 = {
#     cat("equation 11 \n")
    coef <- (a+d)/(a+(0.5*(b+c))+d)
  },
  eq_12 = {
#     cat("equation 12 \n")
    coef <- a
  },
  eq_13 = {
#     cat("equation 13 \n")
    coef <- a+d
  },
  eq_14 = {
#     cat("equation 14 \n")
    coef <- a/n
  },
  eq_15 = {
#     cat("equation 15 \n")
    coef <- b+c
  },
  eq_16 = {
#     cat("equation 16 \n")
    coef <- sqrt(b+c)
  },
  eq_17 = {
#     cat("equation 17 \n")
    coef <- sqrt((b+c)^2)
  },
  eq_18 = {
#     cat("equation 18 \n")
    coef <- sqrt((b+c)^2)
  },
  eq_19 = {
#     cat("equation 19 \n")
    coef <- b+c
  },
  eq_20 = {
#     cat("equation 20 \n")
    coef <- (b+c)/n
  },
  eq_21 = {
#     cat("equation 21 \n")
    coef <- b+c
  },
  eq_22 = {
#     cat("equation 22 \n")
    coef <- b+c
  },
  eq_23 = {
#     cat("equation 23 \n")
    coef <- (b+c)/(4*n)
  },
  eq_24 = {
#     cat("equation 24 \n")
    coef <- ((b+c)^2)/(n^2)
  },
  eq_25 = {
#     cat("equation 25 \n")
    coef <- ((n*(b+c))-((b-c)^2))/(n^2)
  },
  eq_26 = {
#     cat("equation 26")
    coef <- (4*b*c)/(n^2)
  },
  eq_27 = {
#     cat("equation 27 \n")
    coef <- (b+c)/(2*a+b+c)
  },
  eq_28 = {
#     cat("equation 28 \n")
    coef <- (b+c)/(2*a+b+c)
  },
  eq_29 = {
#     cat("equation 29 \n")
    coef <- 2 * (sqrt(1-(a/sqrt((a+b)*(a+c)))))
  },
  eq_30 = {
#     cat("equation 30 \n")
    coef <- (sqrt(2*(1-(a/(sqrt((a+b)*(a+c)))))))
  },
  eq_31 = {
#     cat("equation 31 \n")
    coef <- a/(sqrt((a+b)*(a+c))^2)
  },
  eq_32 = {
#     cat("equation 32 \n")
    coef <- log(a) - log(n) - log((a+b)/n) - log((a+c)/n)
  },
  eq_33 = {
#     cat("equation 33 \n")
    coef <- a/(sqrt((a+b)*(a+c)))
  },
  eq_34 = {
#     cat("equation 34 \n")
    coef <- (n*a)/((a+b)*(a+c))
  },
  eq_35 = {
#     cat("equation 35 \n")
    coef <- (n*((a-0.5)^2))/((a+b)*(a+c))
  },
  eq_36 = {
#     cat("equation 36 \n")
    coef <-(a^2)/((a+b)*(a+c))
  }, 
  eq_37 = {
#     cat("equation 37 \n")
    coef <- a/(0.5*(a*b + a*c)+(b*c))
  },
  eq_38 = {
#     cat("equation 38 \n")
    coef <- a/(((a+b)*(a+c))^0.5)
  },
  eq_39 = {
#     cat("equation 39 \n")
    coef <- ((a^2)-(b*c))/((a+b)*(a+c))
  },
  eq_40 = {
#     cat("equation 40 \n")
    coef <-((n*a) - (a+b)*(a+c))/((n*a) + (a+b)*(a+c))
  },
  eq_41 = {
#     cat("equation 41 \n")
    coef <- ((a/2)*(2*a+b+c))/((a+b)*(a+c))
  },
  eq_42 = {
#     cat("equation 42 \n")
    coef <- (a/2)*((1/(a+b))+(1/(a+c)))
  },
  eq_43 = {
#     cat("equation 43 \n")
    coef <- (a/(a+b))+(a/(a+c))
  },
  eq_44 = {
#     cat("equation 44 \n")
    coef <- ((a*d)-(b*c))/(sqrt(n*(a+b)*(a+c)))
  },
  eq_45 = {
#     cat("equation 45 \n")
    coef <-a/(min((a+b),(a+c)))
  },
  eq_46 = {
#     cat("equation 46 \n")
    coef <- a/(max((a+b),(a+c)))
  },
  eq_47 = {
#     cat("equation 47 \n")
    coef <- (a/sqrt((a+b)*(a+c)))-(max((a+b),(a+c))/2)
  },
  eq_48 = {
#     cat("equation 48 \n")
    coef <- ((n*a)-((a+b)*(a+c)))/((n*min(a+b,a+c))-(a+b)*(a+c))
  },
  eq_49 = {
#     cat("equation 49 \n")
    coef <- 0.25 * ((a/(a+b))+(a/(a+c))+(d/(b+d))+(d/(c+d)))
  },
  eq_50 = {
#     cat("equation 50 \n")
    coef <- (a+d)/(sqrt((a+b)*(a+c)*(b+d)*(c+d)))
  },
  eq_51 = {
#     cat("equation 51 \n")
    x2 <- (n*(((a*d)-(b*c))^2))/((a+b)*(a+c)*(c+d)*(b+d))
    coef <- x2
  },
  eq_52 = {
#     cat("equation 52 \n")
    x2 <- (n*(((a*d)-(b*c))^2))/((a+b)*(a+c)*(c+d)*(b+d))
    coef <- sqrt(x2/(n+x2))
  },
  eq_53 = {
#     cat("equation 53 \n")
    p <- ((a*d) - (b*c))/(sqrt((a+b)*(a+c)*(b+d)*(c+d)))
    coef <- sqrt(p/(n+p))
  },
  eq_54 = {
#     cat("equation 54 \n")
    p <- ((a*d) - (b*c))/(sqrt((a+b)*(a+c)*(b+d)*(c+d)))
    coef <- p
  },
  eq_55 = {
#     cat("equation 55 \n")
    coef <- cos((pi*sqrt(b*c))/(sqrt(a*d)+sqrt(b*c)))
  },
  eq_56 = {
#     cat("equation 56 \n")
    coef <- (a+d)/(b+c)
  },
  eq_57 = {
#     cat("equation 57 \n")
    coef <- (a*d)/(((a+b)*(a+c)*(b+d)*(c+d))^0.5)
  },
  eq_58 = {
#     cat("equation 58 \n")
    coef <- (sqrt(2) * (a*d - b*c)) / (sqrt( ((a*d - b*c)^2) - (a+b)*(a+c)*(b+d)*(c+d)))
  },
  eq_59 = {
#     cat("equation 59 \n")
    coef <- log10( (n*(((abs(a*d - b*c)) - n/2)^2)) / ((a+b)*(a+c)*(b+d)*(c+d)))
  },
  eq_60 = {
#     cat("equation 60 \n")
    coef <- (a*d)/(sqrt((a+b)*(a+c)*(b+d)*(c+d)))
  },
  eq_61 = {
#     cat("equation 61 \n")
    coef <- (a*d - b*c)/(a*d + b*c)
  },
  eq_62 = {
#     cat("equation 62 \n")
    coef <-((2*b*c)/(a*d + b*c))
  },
  eq_63 = {
#     cat("equation 63 \n")
    coef <- (sqrt(a*d) - sqrt(b*c))/(sqrt(a*d) + sqrt(b*c))
  },
  eq_64 = {
#     cat("equation 64 \n")
    coef <- a/(b+c)
  },
  eq_65 = {
#     cat("equation 65 \n")
    coef <- a/((a+b)+(a+c)-a)
  },
  eq_66 = {
#     cat("equation 66 \n")
    coef <- (a*d - b*c)/(n^2)
  },
  eq_67 = {
#     cat("equation 67 \n")
    coef <- ((a+d)-(b+c))/n
  },
  eq_68 = {
#     cat("equation 68 \n")
    coef <- (4*(a*d - b*c))/(((a+d)^2)+((b+c)^2))
  },
  eq_69 = {
#     cat("equation 69 \n")
    sig <- max(a,b) + max(c,d) + max(a,c) + max(b,d)
    sigt <- max(a+c,b+d) + max(a+b,c+d)
    coef <-(sig-sigt)/(2*n - sigt)
  },
  eq_70 = {
#     cat("equation 70 \n")
    sig <- max(a,b) + max(c,d) + max(a,c) + max(b,d)
    sigt <- max(a+c,b+d) + max(a+b,c+d)
    coef <- (sig - sigt)/(2*n)
  },
  eq_71 = {
#     cat("equation 71 \n")
    coef <- (sqrt(a*d)+a)/(sqrt(a*d)+a+b+c)
  },
  eq_72 = {
#     cat("equation 72 \n")
    coef <- (sqrt(a*d)+a-(b+c))/(sqrt(a*d)+a+b+c)
  },
  eq_73 = {
#     cat("equation 73 \n")
    coef <- (a*b + b*c)/((a*b)+(2*b*c)+(c*d))
  },
  eq_74 = {
#     cat("equation 74 \n")
    coef <- ((n^2) * (n*a - (a+b)*(a+c))) / ((a+b)*(a+c)*(b+d)*(c+d))
  },
  eq_75 = {
#     cat("equation 75 \n")
    coef <- (a*(c+d))/(c*(a+b))
  },
  eq_76 = {
#     cat("equation 76 \n")
    coef <- abs((a*(c+d))/(c*(a+b)))
  },
  eq_77 = {
#     cat("equation 77 \n")
    coef <- log(1+a)/log(1+n)
  },
  eq_78 = {
#     cat("equation 78 \n")
    coef <- log(1+a)/log(1+a+b+c)
  },
  eq_79 = {
#     cat("equation 79 \n")
    coef <-(log(1+a*d)-log(1+b*c))/log(1+(n^2)/4)
  },
  {
    cat("No desired equation. Please check it again.")
  }
  )
  
  result <- c(a,b,c,d,coef)
  result <- t(as.matrix(result))
  colnames(result) <- c("a","b","c","d","coef")
  return(result)
}

