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


#' The Binary Similarity/Dissimilarity Coefficient between two vectors 
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
            coef <- a/(sqrt((a+b)*(a+c)))
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


#' Finding a suitable binary similarity and dissimilarity measures  
#' 
#' @param inFile Input file as a matrix with rows represent samples and columns represent features/variables
#' @param setSeed Parameter to use the seed of random generator
#' @param numSample The number of samples to calculate the mean of AUC
#' @return The list of equations ordered by meanAUC
#' @export 
bmeasures_find <- function (inFile, setSeed=0, numSample=20){
  # read input file
  dataIn <- read.csv(inFile, header = TRUE)
  dataIn <- as.matrix(dataIn)
  ndataIn <- nrow(dataIn)
  cdataIn <- ncol(dataIn)
  
  maxPair <- (ndataIn * (ndataIn - 1))/2
  nFeat <- cdataIn-1
  
  # selected equation
  redEq <- c("eq_01", "eq_02", "eq_04", "eq_06", "eq_07", "eq_08", "eq_09", "eq_10", "eq_12", "eq_15", "eq_16", "eq_24", "eq_25", "eq_26", "eq_27", "eq_29", "eq_31", "eq_34", "eq_35", "eq_36", "eq_39", "eq_40", "eq_44", "eq_45", "eq_46", "eq_47", "eq_48", "eq_49", "eq_50", "eq_51", "eq_52", "eq_54", "eq_55", "eq_57", "eq_59", "eq_61", "eq_62", "eq_63", "eq_66", "eq_68", "eq_71", "eq_74", "eq_77", "eq_78", "eq_79")
  
  # temporary out matrix
  mat <- matrix(0,maxPair,length(redEq)+3)   # cols: Sim/Dissim Eq. 1-45, class1, class2, match/mismatch
  tmpID <- matrix(0,maxPair,2)
  
  # calculate similarity/dissimilarity measures using selected equations
  # for each equation; for all sample pairs
  cat("\nMeasuring the similarity/dissimilarity measures between samples: \n")
  
  id <- 1
  for(i in 1:(ndataIn-1)){
    cat("   Samples:", i,"of", ndataIn, "\n")
    
    for (j in (i+1):ndataIn){
      
      x <- t(as.matrix(dataIn[i,1:nFeat]))
      y <- t(as.matrix(dataIn[j,1:nFeat]))
      
      # calculate the otu table
      otuTable <- bmeasures_otu(x,y)
      a <- otuTable[1]
      b <- otuTable[2]
      c <- otuTable[3]
      d <- otuTable[4]
      
      # match and mismatch class
      classx <- dataIn[i,cdataIn]
      classy <- dataIn[j,cdataIn]
      
      mat[id,length(redEq)+1] <- classx
      mat[id,length(redEq)+2] <- classy
      
      if(classx == classy)
        mat[id,length(redEq)+3] <- 1
      
      # calculate similarity/dissimilarity measures using all equations
      for(eq in 1:length(redEq)){
        mat[id,eq] <- bmeasures_coef(a,b,c,d,method=redEq[eq])
      }      
      
      id <- id+1
    }
  }
  
  colnames(mat) <- c(redEq, "classx", "classy", "matchxy")
  
  # remove column "classx" and "classy"
  mat <- mat[,c(-46,-47)]         # character
  
  # eliminate equations with NA/Infinite values
  mat[is.infinite(mat)] <- NA
  mat <- mat[,!is.na(colSums(mat))]       #remove col with infinite and NA values
  
  label <- colnames(mat)
  z <- length(label)
  dissim <- c("eq_15", "eq_16", "eq_17", "eq_18", "eq_19", "eq_20", "eq_21", "eq_22", "eq_23", "eq_24", "eq_25", "eq_26", "eq_27", "eq_28", "eq_29", "eq_30", "eq_62")
  
  # classify into match and mismatch classes
  oneJ <- subset(mat, mat[,z] == 1)
  zeroJ <- subset(mat, mat[,z] == 0)
  
  
  # ROC Analysis -----------------------------------------------------------------------------------
  cat("\nThe ROC Analysis:\n")
  library(ROCR)
  roneJ <- nrow(oneJ)
  rzeroJ <- nrow(zeroJ)
  
  # num. of samples in the match class are smaller than num. of samples in the mismatch class
  if(roneJ <= rzeroJ){
    data_match <- matrix(1,roneJ,2)       # Sim coeff values, label (1/0)
    data_mismatch <- matrix(0,roneJ,2)
  } else {
    data_match <- matrix(1,rzeroJ,2)
    data_mismatch <- matrix(0,rzeroJ,2)
  }
  eq_auc <- matrix(0,(z-1),(numSample+3))              # eqIDs, eqName, type (S/D), auc of 20 samples
  
  for (eq in 1:(z-1)){
    cat("   Processing:", label[eq], "(", eq, "/", (z-1), ")\n")
    
    maxSim1J <- max(oneJ[,eq])
    minSim1J <- min(oneJ[,eq])
    
    maxSim0J <- max(zeroJ[,eq])
    minSim0J <- min(zeroJ[,eq])
    
    maxSimJ <- max(maxSim1J, maxSim0J)
    minSimJ <- min(minSim1J, minSim0J)
    
    # normalize using (X-min)/(max-min)
    oneJ[,eq] <- as.numeric(as.matrix((oneJ[,eq]-minSimJ)/(maxSimJ-minSimJ)))
    zeroJ[,eq] <- as.numeric(as.matrix((zeroJ[,eq]-minSimJ)/(maxSimJ-minSimJ)))
    
    # sim/dissim measures
    # S = 1 - D^2 
    mysim <- label[eq] %in% dissim
    if (mysim == "TRUE"){
      type <- "D"
      
      # S = 1 - D^2
      oneJ[,eq] <- 1 - (oneJ[,eq]^2)
      zeroJ[,eq] <- 1- (zeroJ[,eq]^2)
    } else {
      type <- "S"
    }
    
    
    # ROCR package (AUC) ------------------------   
    # Update Sept 5, 2014 : same size (balanced match and mismatch class) and iteration numSample 
    
    if(roneJ <= rzeroJ){
      data_match[,1] <- oneJ[,eq]
    } else {
      data_mismatch[,1] <- zeroJ[,eq]
    }
    
    for (zsample in 1:numSample)
    {
      # cat("   Sample:", zsample, "\n")
      
      if(setSeed==1){
        set.seed(zsample*10)
      }
      
      if(roneJ <= rzeroJ){
        # select sample from zeroJtemp randomly (having the same dimension with match class)    
        zero_sample <- sample(1:rzeroJ, roneJ, replace=FALSE)     # generate random index
        zeroJtemp <- zeroJ[zero_sample,eq]
        
        data_mismatch[,1] <- zeroJtemp
      } else {
        # select sample from oneJtemp randomly (having the same dimension with mismatch class)    
        one_sample <- sample(1:roneJ, rzeroJ, replace=FALSE)     # generate random index
        oneJtemp <- oneJ[one_sample,eq]
        
        data_match[,1] <- oneJtemp
      }
      
      
      rocr <- rbind(data_match, data_mismatch)
      colnames(rocr) <- c("simCoefVal", "labels")
      
      pred <- prediction(rocr[,1], rocr[,2])
      perf <- performance(pred,"tpr","fpr")
      
      # AUC analysis   
      auc <- performance(pred,"auc")
      auc_label <- unlist (auc@y.name)
      auc_values <- round(as.numeric(unlist(auc@y.values)), digits=4)
      
      if(zsample > 1){
        eq_auc[eq,(zsample+3)] <- auc_values
      } else {
        eq_auc[eq,1] <- label[eq]
        eq_auc[eq,2] <- bmeasures_eqname(label[eq])
        eq_auc[eq,3] <- type
        eq_auc[eq,4] <- auc_values  
      }
    }
    # end zsample
    
  }
  # end ROC analysis
  
  # Mean of AUC
  tmp_eqAuc <- as.matrix((eq_auc[,4:ncol(eq_auc)]))
  class(tmp_eqAuc) <- "numeric"
  meanAUC <- as.matrix(rowMeans(tmp_eqAuc))
  
  eq_auc <- cbind(eq_auc, meanAUC)
  iter <- sprintf("iter%02d",1:numSample)
  colnames(eq_auc) <- c("Eq.IDs","Equations", "Sim/Dissim",iter,"meanAUC")
  
  # order by meanAUC descending
  eq_auc <- eq_auc[order(-meanAUC),]
  
  #View(eq_auc)
  
  # View top-n recommended equations
  cat("\nTop-10 recommended equations:\n")
  print(eq_auc[1:10,c(1,2,ncol(eq_auc))])
  
  # save output
  fout <- paste("bmeasures_", inFile)
  fout <- gsub(" ", "", fout)
  write.csv(eq_auc, file = fout, row.names = FALSE)
  
  return (eq_auc)
}


#' The Binary Similarity/Dissimilarity Coefficient using the quantities of OTUs table  
#' 
#' @param a The number of features where the value of both x and y are 1 (positive matches)
#' @param b The number of features where the value of x is 0 and y is 1 (absence mismatches)
#' @param c The number of features where the value of x is 1 and y is 0 (absence mismatches)
#' @param d he number of features where the values of both x and y are 0 (negative matches)
#' @return The binary similarity/dissimilarity coeffifient
#' @export 
bmeasures_coef <- function (a, b, c, d, method)
{
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
            coef <- a/(sqrt((a+b)*(a+c)))
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
  
  #   result <- c(a,b,c,d,coef)
  #   result <- t(as.matrix(result))
  #   colnames(result) <- c("a","b","c","d","coef")
  return(coef)
}


#' Get equation's name from its ID
#' 
#' @param eqId Equation's ID, i.e. "eq_01", "eq_02", etc
#' @return eqName name of similarity/dissimilarity measures, i.e. "Jaccard similarity", "Dice-2 similarity", etc
#' @export
bmeasures_eqname <- function (eqId){
  
  switch (eqId,
          eq_01 = {
            eqName <- c("Jaccard similarity")
          },
          eq_02 = {
            eqName <- c("Dice-2 similarity")
          },
          eq_03 = {
            eqName <- c("Dice-1/Czekanowski similarity")
          },
          eq_04 = {
            eqName <- c("3W-Jaccard similarity")
          },
          eq_05 = {
            eqName <- c("Nei & Li similarity")
          },
          eq_06 = {
            eqName <- c("Sokal & Sneath-1 similarity")
          },
          eq_07 = {
            eqName <- c("Sokal & Michener similarity")
          },
          eq_08 = {
            eqName <- c("Sokal & Sneath-2 similarity")
          },
          eq_09 = {
            eqName <- c("Roger & Tanimoto similarity")
          },
          eq_10 = {
            eqName <- c("Faith similarity")
          },
          eq_11 = {
            eqName <- c("Gower & Legendre similarity")
          },
          eq_12 = {
            eqName <- c("Intersection similarity")
          },
          eq_13 = {
            eqName <- c("Inner product similarity")
          },
          eq_14 = {
            eqName <- c("Russell & Rao similarity")
          },
          eq_15 = {
            eqName <- c("Hamming distance")
          },
          eq_16 = {
            eqName <- c("Euclid distance")
          },
          eq_17 = {
            eqName <- c("Squared-euclid distance")
          },
          eq_18 = {
            eqName <- c("Canberra distance")
          },
          eq_19 = {
            eqName <- c("Manhattan distance")
          },
          eq_20 = {
            eqName <- c("Mean-Manhattan distance")
          },
          eq_21 = {
            eqName <- c("Cityblock distance")
          },
          eq_22 = {
            eqName <- c("Minkowski distance")
          },
          eq_23 = {
            eqName <- c("Vari distance")
          },
          eq_24 = {
            eqName <- c("Size Difference distance")
          },
          eq_25 = {
            eqName <- c("Shape Difference distance")
          },
          eq_26 = {
            eqName <- c("Pattern Difference distance")
          },
          eq_27 = {
            eqName <- c("Lance & Williams distance")
          },
          eq_28 = {
            eqName <- c("Bray & Curtis distance")
          },
          eq_29 = {
            eqName <- c("Hellinger distance")
          },
          eq_30 = {
            eqName <- c("Chord distance")
          },
          eq_31 = {
            eqName <- c("Cosine similarity")
          },
          eq_32 = {
            eqName <- c("Gilbert & Wells similarity")
          },
          eq_33 = {
            eqName <- c("Ochiai-1 similarity")
          },
          eq_34 = {
            eqName <- c("Forbes-1 similarity")
          },
          eq_35 = {
            eqName <- c("Fossum similarity")
          },
          eq_36 = {
            eqName <- c("Sorgenfrei similarity")
          }, 
          eq_37 = {
            eqName <- c("Mountford similarity")
          },
          eq_38 = {
            eqName <- c("Otsuka similarity")
          },
          eq_39 = {
            eqName <- c("Mc Connaughey similarity")
          },
          eq_40 = {
            eqName <- c("Tarwid similarity")
          },
          eq_41 = {
            eqName <- c("Kulczynski-2 similarity")
          },
          eq_42 = {
            eqName <- c("Driver & Kroeber similarity")
          },
          eq_43 = {
            eqName <- c("Johnson similarity")
          },
          eq_44 = {
            eqName <- c("Dennis similarity")
          },
          eq_45 = {
            eqName <- c("Simpson similarity")
          },
          eq_46 = {
            eqName <- c("Braun & Banquet similarity")
          },
          eq_47 = {
            eqName <- c("Fager & McGowan similarity")
          },
          eq_48 = {
            eqName <- c("Forbes-2 similarity")
          },
          eq_49 = {
            eqName <- c("Sokal & Sneath-4 similarity")
          },
          eq_50 = {
            eqName <- c("Gower similarity")
          },
          eq_51 = {
            eqName <- c("Pearson-1 similarity")
          },
          eq_52 = {
            eqName <- c("Pearson-2 similarity")
          },
          eq_53 = {
            eqName <- c("Pearson-3 similarity")
          },
          eq_54 = {
            eqName <- c("Pearson & Heron-1 similarity")
          },
          eq_55 = {
            eqName <- c("Pearson & Heron-2 similarity")
          },
          eq_56 = {
            eqName <- c("Sokal & Sneath-3 similarity")
          },
          eq_57 = {
            eqName <- c("Sokal & Sneath-5 similarity")
          },
          eq_58 = {
            eqName <- c("Cole similarity")
          },
          eq_59 = {
            eqName <- c("Stiles similarity")
          },
          eq_60 = {
            eqName <- c("Ochiai-2 similarity")
          },
          eq_61 = {
            eqName <- c("Yuleq similarity")
          },
          eq_62 = {
            eqName <- c("Yuleq distance")
          },
          eq_63 = {
            eqName <- c("Yulew similarity")
          },
          eq_64 = {
            eqName <- c("Kulczynski-1 similarity")
          },
          eq_65 = {
            eqName <- c("Tanimoto similarity")
          },
          eq_66 = {
            eqName <- c("Disperson similarity")
          },
          eq_67 = {
            eqName <- c("Hamann similarity")
          },
          eq_68 = {
            eqName <- c("Michael similarity")
          },
          eq_69 = {
            eqName <- c("Goodman&Kruskal similarity")
          },
          eq_70 = {
            eqName <- c("Anderberg similarity")
          },
          eq_71 = {
            eqName <- c("Baroni-Urbani&Buser-1 similarity")
          },
          eq_72 = {
            eqName <- c("Baroni-Urbani&Buser-2 similarity")
          },
          eq_73 = {
            eqName <- c("Peirce similarity")
          },
          eq_74 = {
            eqName <- c("Eyraud similarity")
          },
          eq_75 = {
            eqName <- c("Tarantula similarity")
          },
          eq_76 = {
            eqName <- c("Ample similarity")
          },
          eq_77 = {
            eqName <- c("Derived Rusell-Rao similarity")
          },
          eq_78 = {
            eqName <- c("Derived Jaccard similarity")
          },
          eq_79 = {
            eqName <- c("Variant of Correlation similarity")
          },
          {
            cat("Unoknown.")
          }
  )
  
  
  return (eqName)
}
