#' @title Parallel Safe Bayesian Additive Regression Trees
#'
#' @description A parallelized implementation of safeBART
#' @param lambda A real number between 0 and 1 that determines the splitting probability in the prior (which is used as the importance sampler of tree models). Quadrianto and Ghahramani (2015) recommend a value less than 0.5 .
#' @param num_trees The number of trees to be sampled.
#' @param seed The seed for random number generation.
#' @param num_cats The number of possible values for the outcome variable.
#' @param y The training data vector of outcomes.
#' @param original_datamat The original training data. Currently all variables must be continuous. The training data does not need to be transformed before being entered to this function.
#' @param alpha_parameters Vector of prior parameters.
#' @param beta_par The power to which the likelihood is to be raised. For BMA, set beta_par=1.
#' @param original_datamat The original test data. This matrix must have the same number of columns (variables) as the training data. Currently all variables must be continuous. The test data does not need to be transformed before being entered to this function.
#' @param ncores The number of cores to be used in parallelization.
#' @param valid_trees If equal to 1, restrict splits so that they describe feasible/valid partitions. e.g. can't have a rule x1<0.75 as the splitting rule for the left child node of a parent node with the splitting rule x1<0.5
#' @param tree_prior 1 = BART prior, 2= spike-and-tree, otherwise default prior by Novi and Quandrianto
#' @param imp_sampler Importance sampler for trees. 1 = BART prior, 2= spike-and-tree, otherwise default prior by Novi and Quandrianto
#' @param alpha_BART The alpha parameter for the standard BART prior.
#' @param beta_BART The beta parameter for the standard BART prior.
#' @param fast_approx If equal to 1, use an approximate BIC weighted average and do not invert matrices for each model (should also use SVD).
#' @return A matrix of probabilities with the number of rows equl to the number of test observations and the number of columns equal to the number of possible outcome categories.
#' @useDynLib safeBart, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @examples
#' beta_par <- 0.5
#'
#' N <- 100
#' p<- 5
#' set.seed(100)
#'
#' epsilon <- rnorm(N)
#'
#' xcov <- matrix(runif(N*p), nrow=N)
#'
#' y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon
#' # <- rep(1,N) + epsilon
#'
#' epsilontest <- rnorm(N)
#'
#' xcovtest <- matrix(runif(N*p), nrow=N)
#' ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+5*xcovtest[,5]+epsilontest
#' #ytest <- rep(1,N) + epsilontest
#'
#'
#'
#'
#' Num_split_vars <- 10
#'
#' lambda <- 0.45
#' Num_models <- 10000
#' num_trees1 <- 5
#'
#' seed1 <- 42
#' ncores <- 7
#'
#'
#'
#' examplepreds1 <- safeBart_parallel(seed1,
#'                                    y, xcov,xcovtest,
#'                                    lambda=0.45,
#'                                    num_models=Num_models,
#'                                    num_trees=num_trees1,
#'                                    beta_par=beta_par,
#'                                    ncores=ncores,
#'                                    outsamppreds=1,
#'                                    nu=3,
#'                                    a=3,
#'                                    sigquant=0.9,
#'                                    valid_trees=1)
#'
#' cbind(examplepreds1,ytest )
#' @export

approx_ProbitBARTIS <- function(seed,
                              y,
                              original_datamat,test_datamat,
                              lambda=0.45,
                              num_models=1000,
                              num_trees=5,
                              beta_par=1,
                              ncores=1,
                              outsamppreds=1,
                              nu=3,
                              a=3,
                              sigquant=0.9,
                              valid_trees=1,
                              tree_prior=0,
                              imp_sampler=0,
                              alpha_BART=0.95,
                              beta_BART=2,
                              s_t_hyperprior=1,
                              p_s_t=0.5,
                              a_s_t=1,
                              b_s_t=3,
                              lambda_poisson=10,
                              fast_approx=0){

  if(ncores>num_models ) stop("ncores > num_models")


  if(is.factor(y)) {
    if(length(levels(y)) != 2) stop("y is a factor with number of levels != 2")
    #binary = TRUE
    #  y.train = as.numeric(y.train)-1
    #stop("Response must be a numeric vector")
  } else {
    if(!((length(unique(y)) == 2) & (max(y) == 1) & (min(y) == 0))) {
      stop("y.train is a numeric with number of unique elements != 2")
      #cat('NOTE: assumming numeric response is binary\n')
      #binary = TRUE
      #stop("Response must be a numeric vector")
    }
  }


  Zlatent.train <- ifelse(y==1,3.1,-3.1)


  sigma=sd(Zlatent.train)/(max(Zlatent.train)-min(Zlatent.train))
  qchi = qchisq(1.0-sigquant,nu,1,0);
  lambdaBART = (sigma*sigma*qchi)/nu;

  if(is.vector(original_datamat) | is.factor(original_datamat)| is.data.frame(original_datamat)) original_datamat = as.matrix(original_datamat)
  if(is.vector(test_datamat) | is.factor(test_datamat)| is.data.frame(test_datamat)) test_datamat = as.matrix(test_datamat)

  if((!is.matrix(original_datamat))) stop("argument x.train must be a double matrix")
  if((!is.matrix(test_datamat)) ) stop("argument x.test must be a double matrix")

  if(nrow(original_datamat) != length(Zlatent.train)) stop("number of rows in x.train must equal length of Zlatent.train.train")
  if((ncol(test_datamat)!=ncol(original_datamat))) stop("input x.test must have the same number of columns as x.train")

  sBARToutput=sBART_onefunc_parallel(lambda,
                                     num_models,
                                     num_trees,
                                     seed,
                                     Zlatent.train,
                                     original_datamat,
                                     beta_par,
                                     test_datamat,
                                     ncores,
                                     outsamppreds,
                                     nu,
                                     a,
                                     lambdaBART,
                                     valid_trees,
                                     tree_prior,
                                     imp_sampler,
                                     alpha_BART,
                                     beta_BART,
                                     s_t_hyperprior,
                                     p_s_t,
                                     a_s_t,
                                     b_s_t,
                                     lambda_poisson,
                                     fast_approx)


  #names(sBARToutput) <- c("Predictions_latent")

  ret <- list()
  ret[[1]] <- sBARToutput

  ret[[2]] <- pnorm(sBARToutput)

  names(ret) <- c("Predictions_latent", "pred_probs")

  ret

  #sBARToutput

}
