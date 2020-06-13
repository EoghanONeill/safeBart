#' @title Variable inclusion probabilities as defined by Linero (2018)
#'
#' @description This measure defines the posterior inclusion probability of a variable as the model-probability weighted sum of indicator variables for whether the variable was used in any splitting rules in any of the trees in the sum-of-tree model.
#' @param object A BCF-IS object obtained using the train_BCF_IS or train_BCF_IS_no_output function.
#' @export
#' @return A vector of posterior inclusion probabilities. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input BCF-IS object.
#' @examples
#' #set the seed
#' set.seed(100)
#' #simulate some data
#' N <- 100
#' p<- 100
#' epsilon <- rnorm(N)
#' xcov <- matrix(runif(N*p), nrow=N)
#' y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon
#' epsilontest <- rnorm(N)
#' xcovtest <- matrix(runif(N*p), nrow=N)
#' ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+
#'   5*xcovtest[,5]+epsilontest
#'
#' #Train the object
#' bart_bma_example <- bartBMA(x.train = xcov,y.train=y,x.test=xcovtest,zero_split = 1,
#'                             only_max_num_trees = 1,split_rule_node = 0)
#' #Obtain the variable importances
#' varIncProb(bart_bma_example)

varIncProb_bcfis<-function(object){
  #object will be BCF-IS object.

  num_to_add_mu = 0
  num_to_add_tau = 0

  if(object$include_pi2 == 0) {
    #include_pi2 = 0
    # control
    num_to_add_mu = num_to_add_mu + object$num_ps
  }
  if(object$include_pi2 == 1) {
    #include_pi2 = 1
    #moderate

    num_to_add_tau = num_to_add_tau + object$num_ps

  }
  if(object$include_pi2 == 2) {
    #include_pi2 = 2
    # both

    num_to_add_mu = num_to_add_mu + object$num_ps
    num_to_add_tau = num_to_add_tau + object$num_ps

  }
  if(object$include_pi2 == 4) {
    #include_pi2 = 4
    # none

  }

  imp_vars2_mu=get_weighted_var_imp(num_vars=object$numvars + num_to_add_mu,
                                    BIC=object$model_probs,
                                    sum_trees=object$sumoftrees_mu)
  #res<-apply((imp_vars2[[3]]>0)*imp_vars2[[1]],2,sum)
  res_mu <- t(imp_vars2_mu[[3]]>0)%*%imp_vars2_mu[[1]]

  imp_vars2_tau=get_weighted_var_imp(num_vars=object$numvars + num_to_add_tau,
                                     BIC=object$model_probs,
                                     sum_trees=object$sumoftrees_tau)
  #res<-apply((imp_vars2[[3]]>0)*imp_vars2[[1]],2,sum)
  res_tau <- t(imp_vars2_tau[[3]]>0)%*%imp_vars2_tau[[1]]


  ret<-list()
  length(ret)<-2
  ret[[1]] <- res_mu
  ret[[2]] <- res_tau

  class(ret)<-"varPIP_list.bcfIS"
  ret

}
