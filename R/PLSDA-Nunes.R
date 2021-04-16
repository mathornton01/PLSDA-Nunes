# Author: Micah Thornton
# Date: 04-16-2021
# Description: PLSDA-Nunes.R 
#   This is an implementation of the Cleiton Nunes PLS Regression and DA 
#   procedure meant to hold consistently the results of the MatLab implementation
#   of the PLSDA procedure hosted on MatLab central, which was used in a prior 
#   study by Lu, Das, et al. 



#' Partial Least Squares Modeling 
#' 
#' This function is a wrapper that will allow for generaic modeling applying the
#' traditional Partial Least Squares Procedure where the response is either a 
#' numerical response, or an outcome class (which can have as many classes as
#' possible). 
#' @param x The matrix of independent variables, or stimuli assocated with the 
#' responses in y.  The matrix should have a number of rows equal to the number 
#' of observations, and a number of columns equal to the number of variables. 
#' @param y The matrix of dependent variables, or responses assocaited with the 
#' stimuli in x.  The matrix should have a number of rowa equal to the number of 
#' observations, and a number of columns equal to the number of variables. 
#' @param vl The number of latent variables to include in the model, this is 
#' essentially the number of times that the iterated least squares procedure is 
#' performed and the regression coefficients captured. 
#' @param da a boolean indicator, which indicates whether the output variables 
#' in this case should be handeled as classes or as continuous variables for 
#' regression true means use discriminant analysis, and treat the outcomes/responses
#' as class variables. 
#' @return A list of several different values containing the results of the PLSDA
#' procedure. 
#' @details They are: 
#' \itemize{
#' \item{"x"}{The original independent matrix}
#' \item{"y"}{The original dependent matrix}
#' \item{"type"}{The type of analysis discriminant analysis (da) or regression (re)}
#' \item{"VLvar"}{cumulative variance (%) explained by model for X and Y}
#' \item{"Ypc"}{Predicted variables or predicted classes}
#' \item{"Tm"}{x-scores}
#' \item{"P"}{x-loadings}
#' \item{"W"}{x-weights}
#' \item{"U"}{y-scores}
#' \item{"Q"}{y-loadings}
#' \item{"B"}{Regression vectors (betas) for the PLSR/DA}
#' }
#' @examples 
#' Prepare the data from acetylene.mat in matlab for testing. 
#' x1 <- c(rep(1300,6),rep(1200,6),rep(1100,4))
#' x2 <- c(7.5,9,11,13.5,17,23,5.3,7.5,11,13.5,17,23,5.3,7.5,11,17); 
#' x3 <- c(0.012,0.012,0.0115,0.013,0.0135,0.012,0.04,0.038,0.0320,0.0260,
#'        0.0340,0.0410,0.0840,0.0980,0.0920,0.0860);
#' y <- c(49, 50.2,50.5,48.5,47.5,44.5,28,31.5,34.5,35,38,38.5,15,17,20.5,29.5)
#' ybin <- y < mean(y);
#' ythree <- acetylene <- data.frame(x1=x1,x2=x2,x3=x3,y=y,ybin=ybin);
#' pls_model <- PLS(acetylene[,1:3],acetylene[,5],2,T);
#' @export
PLS <- function(x,y,vl,da){
  # x <- factor matrix, n observations by p variables 
  # y <- response matrix, n observations, by k variables
  # vl <- number of latent variables to model 
  # da <- boolean indicating whether discriminant analysis or not
  yo <- y # Don't modify original y
  
  
  if (da){ # If discriminant analysis
    if (is.vector(y) != 1) { # Cannot handle mutli class discriminant analysis
      errorCondition("Not currently capable of Multiclassification discriminant 
                     analysis"); 
    }
    else {
      u <- unique(y) # Get the classes from y 
      t <- length(u) # The number of unique classes in y
      sc <- length(y) # The number of observations 
      yc <- matrix(0,nrow=sc,ncol=t) # For storing unit vectors indicating class
      for (i in 1:sc){
        for (j in 1:t){
          if (y[i] == u[j]){
            yc[i,j] <- 1; 
          }
        }
      }
    }
    y <- yc;
  }
  
  plsr_res <- PLSR(x,y,vl); # Run the actual PLSR code
  #DEBUG: print(plsr_res); # Checks if proper calculation completed 
  ypc <- as.matrix(cbind(Int=rep(1,nrow(x)),x)) %*% plsr_res$B; 
  
  if (da){
    nc <- apply(abs(ypc-1), 1, function(x) {which.min(x)}); 
    d <- apply(abs(ypc-1), 1, function(x) {min(x)}); 
    nnc <- nc; 
    for (i in 1:sc){
      for (j in 1:t){
        if (nc[i] == j){
          nnc[i] <- u[j]; 
        }
      }
    }
    ypc <- nnc;
  }
  s <- nrow(as.matrix(ypc)); 
  v <- ncol(as.matrix(ypc));
  rmsec <- matrix(0,1,v); 
  suc <- c(); 
  ypcm <- as.matrix(ypc);
  yom <- as.matrix(yo);
  for (inv in 1:v){
    if (da){
      suc <- c(suc,1-(sum((ypcm[,inv]-yom[,inv]) != 0)/nrow(as.matrix(ypc)))*100);
    } else {
      rmsec[1,inv] <- sqrt(sum((ypcm[,inv]-yom[,inv])^2)/nrow(as.matrix(ypc)));
    } 
  }
  if (da){
    return(list(X=x,Y=y,class=u,type='da',VLvar=plsr_res$var_LV, YPC = ypc, Tm=plsr_res$Tm, 
                P = plsr_res$P, W = plsr_res$W, U = plsr_res$U, Q= plsr_res$Q, B=plsr_res$B))
  } else {
    return(list(X=x,Y=y,type='re',VLvar=plsr_res$var_LV, YPC = ypc, Tm=plsr_res$Tm, 
                P = plsr_res$P, W = plsr_res$W, U = plsr_res$U, Q= plsr_res$Q, B=plsr_res$B))
  }
}

PLSR <- function(x,y,lv){
  n <- nrow(x) # Number of observations 
  m <- ncol(x) # Number of factor variables 
  l <- nrow(y) # Number of observations, same as n
  k <- ncol(y) # Number of response variables (in DA number of levels)
  P <- matrix(0, nrow = m, ncol = lv) # The X-Loading Matrix. num factors by latent variables
  Q <- matrix(0, nrow = k, ncol = lv) # The Y-Loading Matrix, num variables by latent variables 
  Tm <- matrix(0, nrow = n, ncol = lv) # X-Score Matrix (projection of the X's)
  U <- matrix(0, nrow = n, ncol = lv) # Y-Score Matrix (projection of the Y's)
  W <- matrix(0, nrow = m, ncol = lv) # UNKNOWN number of variables by latent variables
  V <- matrix(0, nrow = m, ncol = lv) # UNKNOWN number of variable by latent variables
  cv <- t(x) %*% y; # First Step
  
  # Iterate through for all latent variables desired 
  for (i in 1:lv){
    # The singular value decomposition of CV is taken, here I will need to determine 
    # what the correspondant return arguments are. 
    svd.ret <- svd(cv);
    ai <- svd.ret$u[,1]; # These values differ slightly from those computed in MatLab due to SVD (weary of the 2 here - may need to be last)
    bi <- svd.ret$d[1];
    ci <- svd.ret$v[,1]; # These values differ slightly from those computed in MatLab due to SVD procedure
    ti <- as.matrix(x) %*% as.matrix(ai); # Compute the current observations scores
    normti <- norm(ti,"2"); # Compute the norm, and normalize the score vector  
    ti <- ti/normti; # normalilze the score vector 
    P[,i] <- t(as.matrix(x)) %*% as.matrix(ti); # Get the X-Loadings for component i
    qi <- bi*ci/normti; # get the Y-Loadings for component 
    Q[,i] <- qi; 
    Tm[,i] <- ti; 
    U[,i] <- y %*% qi; 
    W[,i] <- ai/normti; 
    vi <- P[,i]; 
    
    for (re in 1:2){
      if(i==1){next;}
      for (j in 1:(i-1)){
        vj <- V[,j]; 
        vi <- vi-as.numeric(t(as.matrix(vi)) %*% vj)*vj
      }
    }
    
    vi <- vi/norm(vi,"2");
    V[,i] <- vi; 
    cv <- cv - vi %*% (t(vi) %*% cv); 
    Vi <- V[,1:i];
    if (i==1){cv <- cv - kronecker(Vi, (t(as.matrix(Vi)) %*% cv));}
    else {cv <- cv - Vi %*% (t(as.matrix(Vi)) %*% cv);
    }
    
  }
  
  for (i in 1:lv){
    ui <- U[,i]; 
    for (re in 1:2){
      if (i == 1){next;}
      for (j in 1:(i-1)){
        tj <- Tm[,j]; 
        ui <- ui - as.numeric(t(as.matrix(ui)) %*% as.matrix(tj))*tj
      }
    }
    U[,i] <- ui; 
  }
  
  B <- W %*% t(as.matrix(Q)); 
  B <- rbind((colMeans(y)-colMeans(x)%*%B),B);
  var_LV <- rbind(colSums(P^2)/sum(colSums(x^2)),colSums(Q^2)/sum(colSums(y^2))); 
  var_LV <- apply(t(var_LV*100),2,cumsum); 
  return(list(Tm=Tm,P=P,U=U,Q=Q,W=W,B=B,var_LV=var_LV))
}

runTestPLSR <- function(){
  # Prepare the data from acetylene.mat in matlab for testing. 
  x1 <- c(rep(1300,6),rep(1200,6),rep(1100,4))
  x2 <- c(7.5,9,11,13.5,17,23,5.3,7.5,11,13.5,17,23,5.3,7.5,11,17); 
  x3 <- c(0.012,0.012,0.0115,0.013,0.0135,0.012,0.04,0.038,0.0320,0.0260,
          0.0340,0.0410,0.0840,0.0980,0.0920,0.0860);
  y <- c(49, 50.2,50.5,48.5,47.5,44.5,28,31.5,34.5,35,38,38.5,15,17,20.5,29.5)
  ybin <- y < mean(y);
  ythree <- 
    acetylene <- data.frame(x1=x1,x2=x2,x3=x3,y=y,ybin=ybin);
  
  # Run and return the test results 
  pls_model <- PLS(acetylene[,1:3],acetylene[,5],2,T);
  
  return(pls_model);
}

pls_model <- runTestPLSR();

print(pls_model)