# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;

//######################################################################################################################//


//' @title ECDF transformation of the training data
//'
//' @description Quadrianto and Ghahramani (2015) reccomend the use of the probability intergral transform to transform the continuous input features. The code is edited from https://github.com/dmbates/ecdfExample
//' @param originaldata Training data matrix
// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
NumericMatrix cpptrans_cdf(NumericMatrix originaldata){
  NumericMatrix transformedData(originaldata.nrow(), originaldata.ncol());
  for(int i=0; i<originaldata.ncol();i++){
    NumericVector samp= originaldata(_,i);
    NumericVector sv(clone(samp));
    std::sort(sv.begin(), sv.end());
    double nobs = samp.size();
    NumericVector ans(nobs);
    for (int k = 0; k < samp.size(); ++k)
      ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
    //NumericVector ansnum = ans;
    transformedData(_,i) = (ans+1)/nobs;
  }
  return transformedData;

}

//######################################################################################################################//

//' @title ECDF transformation of the test data
//'
//' @description Quadrianto and Ghahramani (2015) reccomend the use of the probability intergral transform to transform the continuous input features. The code is edited from https://github.com/dmbates/ecdfExample
//' @param originaldata Training data matrix
//' @param testdata Test data matrix
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix cpptrans_cdf_test(NumericMatrix originaldata, NumericMatrix testdata){
  NumericMatrix transformedData(testdata.nrow(), testdata.ncol());
  for(int i=0; i<testdata.ncol();i++){
    NumericVector samp= testdata(_,i);
    NumericVector svtest = originaldata(_,i);
    NumericVector sv(clone(svtest));
    std::sort(sv.begin(), sv.end());
    double nobs = samp.size();
    NumericVector ans(nobs);
    double nobsref = svtest.size();
    for (int k = 0; k < samp.size(); ++k){
      ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
    }
    //NumericVector ansnum = ans;
    transformedData(_,i) = (ans)/nobsref;
  }
  return transformedData;

}
//######################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));

  return(y_scaled);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector get_original(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);

  return(original_y);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_nodes(NumericMatrix tree_table){
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

  //arma::vec colmat=arma_tree.col(4);
  //arma::uvec term_nodes=arma::find(colmat==-1);

  //arma::vec colmat=arma_tree.col(2);
  //arma::uvec term_nodes=arma::find(colmat==0);

  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==0);

  term_nodes=term_nodes+1;

  return(wrap(term_nodes));
}

//######################################################################################################################//

#include <math.h>       /* tgamma */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_treepreds(NumericVector original_y, int num_cats, NumericVector alpha_pars,
                   NumericMatrix originaldata, //NumericMatrix test_data,
                   NumericMatrix treetable//, NumericMatrix tree_data
) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.

  //test_data is a nxp matrix with the same variable names as the training data the model was built on

  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values

  //term_node_means is a vector storing the terminal node mean values
  arma::vec orig_y_arma= as<arma::vec>(original_y);
  arma::vec alpha_pars_arma= as<arma::vec>(alpha_pars);

  double lik_prod=1;
  double alph_prod=1;
  for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
    alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
  }
  double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
  double alph_term=gam_alph_sum/alph_prod;

  arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
  arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


  //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

  //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

  NumericVector terminal_nodes=find_term_nodes(treetable);
  //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  //NumericVector tree_predictions;

  //now for each internal node find the observations that belong to the terminal nodes

  //NumericVector predictions(test_data.nrow());
  //List term_obs(terminal_nodes.size());

  if(terminal_nodes.size()==1){
    //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    //predictions=rep(nodemean,test_data.nrow());
    //Rcout << "Line 67 .\n";

    //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
    //term_obs[0]= temp_obsvec;
    double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

    double num_prod=1;
    double num_sum=0;
    //Rcout << "Line 129.\n";

    for(int k=0; k<num_cats; k++){
      //assuming categories of y are from 1 to num_cats
      arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
      double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
      arma_tree_table(0,5+k)= m_plus_alph/denom_temp ;

      //for likelihood calculation
      num_prod=num_prod*tgamma(m_plus_alph);
      num_sum=num_sum +m_plus_alph ;
    }

    lik_prod= alph_term*num_prod/tgamma(num_sum);

  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      //arma::mat subdata=testd;
      int curr_term=terminal_nodes[i];

      int row_index;
      int term_node=terminal_nodes[i];
      //Rcout << "Line 152.\n";


      //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
      //Why should the ro index be different for a right daughter?
      //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
      row_index=0;

      // if(curr_term % 2==0){
      //   //term node is left daughter
      //   row_index=terminal_nodes[i];
      // }else{
      //   //term node is right daughter
      //   row_index=terminal_nodes[i]-1;
      // }




      //save the left and right node data into arma uvec

      //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
      //arma::vec left_nodes=arma_tree.col(0);
      //arma::vec right_nodes=arma_tree.col(1);

      arma::vec left_nodes=arma_tree_table.col(0);
      arma::vec right_nodes=arma_tree_table.col(1);



      arma::mat node_split_mat;
      node_split_mat.set_size(0,3);
      //Rcout << "Line 182. i = " << i << " .\n";

      while(row_index!=1){
        //for each terminal node work backwards and see if the parent node was a left or right node
        //append split info to a matrix
        int rd=0;
        arma::uvec parent_node=arma::find(left_nodes == term_node);

        if(parent_node.size()==0){
          parent_node=arma::find(right_nodes == term_node);
          rd=1;
        }

        //want to cout parent node and append to node_split_mat

        node_split_mat.insert_rows(0,1);

        //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
        //node_split_mat(0,0)=treetable(parent_node[0],2);
        //node_split_mat(0,1)=treetable(parent_node[0],3);

        //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
        //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

        node_split_mat(0,0)=arma_tree_table(parent_node[0],2);
        node_split_mat(0,1)=arma_tree_table(parent_node[0],3);

        node_split_mat(0,2)=rd;
        row_index=parent_node[0]+1;
        term_node=parent_node[0]+1;
      }

      //once we have the split info, loop through rows and find the subset indexes for that terminal node!
      //then fill in the predicted value for that tree
      //double prediction = tree_data(term_node,5);
      arma::uvec pred_indices;
      int split= node_split_mat(0,0)-1;

      //Rcout << "Line 224.\n";
      //Rcout << "split = " << split << ".\n";
      //arma::vec tempvec = testd.col(split);
      arma::vec tempvec = arma_orig_data.col(split);
      //Rcout << "Line 227.\n";


      double temp_split = node_split_mat(0,1);

      if(node_split_mat(0,2)==0){
        pred_indices = arma::find(tempvec <= temp_split);
      }else{
        pred_indices = arma::find(tempvec > temp_split);
      }
      //Rcout << "Line 236.\n";

      arma::uvec temp_pred_indices;

      //arma::vec data_subset = testd.col(split);
      arma::vec data_subset = arma_orig_data.col(split);

      data_subset=data_subset.elem(pred_indices);

      //now loop through each row of node_split_mat
      int n=node_split_mat.n_rows;
      //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
      //Rcout << "Line 248.\n";

      for(int j=1;j<n;j++){
        int curr_sv=node_split_mat(j,0);
        double split_p = node_split_mat(j,1);

        //data_subset = testd.col(curr_sv-1);
        //Rcout << "Line 255.\n";
        //Rcout << "curr_sv = " << curr_sv << ".\n";
        data_subset = arma_orig_data.col(curr_sv-1);
        //Rcout << "Line 258.\n";

        data_subset=data_subset.elem(pred_indices);

        if(node_split_mat(j,2)==0){
          //split is to the left
          temp_pred_indices=arma::find(data_subset <= split_p);
        }else{
          //split is to the right
          temp_pred_indices=arma::find(data_subset > split_p);
        }
        pred_indices=pred_indices.elem(temp_pred_indices);

        if(pred_indices.size()==0){
          continue;
        }

      }
      //Rcout << "Line 199. i = " << i <<  ".\n";

      //double nodemean=tree_data(terminal_nodes[i]-1,5);
      //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
      //predictions[predind]= nodemean;
      //term_obs[i]=predind;

      double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
      //Rcout << "Line 207. predind = " << predind <<  ".\n";
      //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
      // << "Line 207. term_node = " << term_node <<  ".\n";

      double num_prod=1;
      double num_sum=0;

      for(int k=0; k<num_cats; k++){
        //assuming categories of y are from 1 to num_cats
        arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
        double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);

        arma_tree_table(curr_term-1,5+k)= m_plus_alph/denom_temp ;

        num_prod=num_prod*tgamma(m_plus_alph);
        num_sum=num_sum +m_plus_alph ;
      }


      lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);
      //Rcout << "Line 297.\n";


    }
    //Rcout << "Line 301.\n";

  }
  //List ret(1);
  //ret[0] = term_obs;

  //ret[0] = terminal_nodes;
  //ret[1] = term_obs;
  //ret[2] = predictions;
  //return(term_obs);
  //Rcout << "Line 309";

  //return(wrap(arma_tree_table));

  List ret(2);
  ret[0]=wrap(arma_tree_table);
  ret[1]=lik_prod;

  return(ret);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
//' @title For a set of trees, obtain tree matrices with predictions, and obtain model weights
//' @export
// [[Rcpp::export]]
List get_treelist(NumericVector original_y, int num_cats, NumericVector alpha_pars,
                  double beta_pow,
                  NumericMatrix originaldata, //NumericMatrix test_data,
                  List treetable_list//, NumericMatrix tree_data
){

  //List overall_term_nodes_trees(overall_sum_trees.size());
  //List overall_term_obs_trees(overall_sum_trees.size());
  //List overall_predictions(overall_sum_trees.size());

  List overall_treetables(treetable_list.size());
  NumericVector overall_liks(treetable_list.size());

  for(int i=0;i<treetable_list.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = treetable_list[i];

    //NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      //List sum_tree=treetable_list[i];

      //save all info in list of list format the same as the trees.
      //List term_nodes_trees(sum_tree.size());
      //List term_obs_trees(sum_tree.size());
      //NumericMatrix predictions(num_obs,sum_tree.size());

      // for(int k=0;k<sum_tree.size();k++){
      //   NumericMatrix tree_table=sum_tree[k];
      //   List tree_info=get_termobs_test_data(test_data, tree_table) ;
      //   //NumericVector term_nodes=tree_info[0];
      //   //term_nodes_trees[k]=term_nodes;
      //   term_obs_trees[k]=tree_info;
      //   //umericVector term_preds=tree_info[2];
      //   //predictions(_,k)=term_preds;
      // }


      List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
                                           originaldata,
                                           treetable_list[i]  );

      overall_treetables[i]= treepred_output[0];
      double templik = as<double>(treepred_output[1]);
      overall_liks[i]= pow(templik,beta_pow);



      //overall_term_nodes_trees[i]=term_nodes_trees;
      //overall_term_obs_trees[i]= term_obs_trees;
      //overall_predictions[i]=predictions;
    }else{
      // NumericMatrix sum_tree=overall_sum_trees[i];
      // List tree_info=get_termobs_test_data(test_data, sum_tree) ;
      // //overall_term_nodes_trees[i]=tree_info[0];
      // List term_obs_trees(1);
      // term_obs_trees[0]=tree_info ;
      // //NumericVector term_preds=tree_info[2];
      // //NumericVector predictions=term_preds;
      // overall_term_obs_trees[i]= term_obs_trees;
      // //overall_predictions[i]=predictions;
      //

      //overall_treetables[i]=get_treepreds(original_y, num_cats, alpha_pars,
      //                                    originaldata,
      //                                    treetable_list[i]  );


      List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
                                           originaldata,
                                           treetable_list[i]  );

      overall_treetables[i]= treepred_output[0];
      double templik = as<double>(treepred_output[1]);
      overall_liks[i]= pow(templik,beta_pow);

    }
  }
  //List ret(1);
  //ret[0]=overall_term_nodes_trees;
  //ret[0]=overall_term_obs_trees;
  //ret[2]=overall_predictions;
  //return(overall_term_obs_trees);

  //return(overall_treetables);

  overall_liks=overall_liks/sum(overall_liks);

  List ret(2);
  ret[0]=overall_treetables;
  ret[1]=overall_liks;
  return(ret);



}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat get_test_probs(NumericVector weights, int num_cats,
                         NumericMatrix testdata, //NumericMatrix test_data,
                         NumericMatrix treetable//, NumericMatrix tree_data
) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.

  //test_data is a nxp matrix with the same variable names as the training data the model was built on

  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values

  //term_node_means is a vector storing the terminal node mean values
  // arma::vec orig_y_arma= as<arma::vec>(original_y);
  // arma::vec alpha_pars_arma= as<arma::vec>(alpha_pars);
  //
  // double lik_prod=1;
  // double alph_prod=1;
  // for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
  //   alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
  // }
  // double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
  // double alph_term=gam_alph_sum/alph_prod;

  arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
  arma::mat arma_test_data(testdata.begin(), testdata.nrow(), testdata.ncol(), false);


  //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

  //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

  NumericVector terminal_nodes=find_term_nodes(treetable);
  //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  //NumericVector tree_predictions;

  //now for each internal node find the observations that belong to the terminal nodes

  //NumericVector predictions(test_data.nrow());

  arma::mat pred_mat(testdata.nrow(),num_cats);
  //arma::vec filled_in(testdata.nrow());


  //List term_obs(terminal_nodes.size());
  if(terminal_nodes.size()==1){

    //Rcout << "Line 422. \n";


    pred_mat=repmat(arma_tree_table(0,arma::span(5,5+num_cats-1)),testdata.nrow(),1);


    //Rcout << "Line 424. \n";


    // for(int k=0; k<num_cats; k++){
    // pred_mat(_,k)=rep(treetable(0,5+k),testdata.nrow());
    // }
    //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    //predictions=rep(nodemean,test_data.nrow());
    //Rcout << "Line 67 .\n";

    //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
    //term_obs[0]= temp_obsvec;
    // double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);
    //
    // double num_prod=1;
    // double num_sum=0;

    // for(int k=0; k<num_cats; k++){
    //   //assuming categories of y are from 1 to num_cats
    //   arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
    //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
    //   arma_tree_table(0,5+k)= m_plus_alph/denom_temp ;
    //
    //   //for likelihood calculation
    //   num_prod=num_prod*tgamma(m_plus_alph);
    //   num_sum=num_sum +m_plus_alph ;
    // }
    //
    // lik_prod= alph_term*num_prod/tgamma(num_sum);
    //
  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      //arma::mat subdata=testd;
      int curr_term=terminal_nodes[i];

      int row_index;
      int term_node=terminal_nodes[i];


      //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
      //Why should the ro index be different for a right daughter?
      //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
      row_index=0;

      // if(curr_term % 2==0){
      //   //term node is left daughter
      //   row_index=terminal_nodes[i];
      // }else{
      //   //term node is right daughter
      //   row_index=terminal_nodes[i]-1;
      // }








      //save the left and right node data into arma uvec

      //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
      //arma::vec left_nodes=arma_tree.col(0);
      //arma::vec right_nodes=arma_tree.col(1);

      arma::vec left_nodes=arma_tree_table.col(0);
      arma::vec right_nodes=arma_tree_table.col(1);



      arma::mat node_split_mat;
      node_split_mat.set_size(0,3);
      //Rcout << "Line 124. i = " << i << " .\n";

      while(row_index!=1){
        //for each terminal node work backwards and see if the parent node was a left or right node
        //append split info to a matrix
        int rd=0;
        arma::uvec parent_node=arma::find(left_nodes == term_node);

        if(parent_node.size()==0){
          parent_node=arma::find(right_nodes == term_node);
          rd=1;
        }

        //want to cout parent node and append to node_split_mat

        node_split_mat.insert_rows(0,1);

        //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
        //node_split_mat(0,0)=treetable(parent_node[0],2);
        //node_split_mat(0,1)=treetable(parent_node[0],3);

        //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
        //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

        node_split_mat(0,0)=arma_tree_table(parent_node[0],2);
        node_split_mat(0,1)=arma_tree_table(parent_node[0],3);

        node_split_mat(0,2)=rd;
        row_index=parent_node[0]+1;
        term_node=parent_node[0]+1;
      }

      //once we have the split info, loop through rows and find the subset indexes for that terminal node!
      //then fill in the predicted value for that tree
      //double prediction = tree_data(term_node,5);
      arma::uvec pred_indices;
      int split= node_split_mat(0,0)-1;

      //arma::vec tempvec = testd.col(split);
      arma::vec tempvec = arma_test_data.col(split);


      double temp_split = node_split_mat(0,1);

      if(node_split_mat(0,2)==0){
        pred_indices = arma::find(tempvec <= temp_split);
      }else{
        pred_indices = arma::find(tempvec > temp_split);
      }

      arma::uvec temp_pred_indices;

      //arma::vec data_subset = testd.col(split);
      arma::vec data_subset = arma_test_data.col(split);

      data_subset=data_subset.elem(pred_indices);

      //now loop through each row of node_split_mat
      int n=node_split_mat.n_rows;
      //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
      //Rcout << "Line 174. node_split_mat= " << node_split_mat << ". n = " << n << ".\n";


      for(int j=1;j<n;j++){
        int curr_sv=node_split_mat(j,0);
        double split_p = node_split_mat(j,1);

        //data_subset = testd.col(curr_sv-1);
        data_subset = arma_test_data.col(curr_sv-1);

        data_subset=data_subset.elem(pred_indices);

        if(node_split_mat(j,2)==0){
          //split is to the left
          temp_pred_indices=arma::find(data_subset <= split_p);
        }else{
          //split is to the right
          temp_pred_indices=arma::find(data_subset > split_p);
        }
        pred_indices=pred_indices.elem(temp_pred_indices);

        if(pred_indices.size()==0){
          continue;
        }

      }
      //Rcout << "Line 199. i = " << i <<  ".\n";

      //double nodemean=tree_data(terminal_nodes[i]-1,5);
      //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
      //predictions[predind]= nodemean;
      //term_obs[i]=predind;

      //Rcout << "Line 635. \n";
      //Rcout << "pred_indices = " << pred_indices << ".\n";

      //pred_mat.rows(pred_indices)=arma::repmat(arma_tree_table(curr_term-1,arma::span(5,5+num_cats-1)),pred_indices.n_elem,1);
      pred_mat.each_row(pred_indices)=arma_tree_table(curr_term-1,arma::span(5,4+num_cats));



      //Rcout << "Line 588. \n";

      // for(int k=0; k<num_cats; k++){
      //   pred_mat(predind,k)=rep(treetable(curr_term-1,5+k),predind.size());
      // }


      // double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
      // //Rcout << "Line 207. predind = " << predind <<  ".\n";
      // //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
      // // << "Line 207. term_node = " << term_node <<  ".\n";
      //
      // double num_prod=1;
      // double num_sum=0;
      //
      // for(int k=0; k<num_cats; k++){
      //   //assuming categories of y are from 1 to num_cats
      //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
      //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
      //
      //   arma_tree_table(curr_term-1,5+k)= m_plus_alph/denom_temp ;
      //
      //   num_prod=num_prod*tgamma(m_plus_alph);
      //   num_sum=num_sum +m_plus_alph ;
      // }


      //lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);


    }
  }
  //List ret(1);
  //ret[0] = term_obs;

  //ret[0] = terminal_nodes;
  //ret[1] = term_obs;
  //ret[2] = predictions;
  //return(term_obs);

  //return(wrap(arma_tree_table));

  //List ret(2);
  //ret[0]=wrap(arma_tree_table);
  //ret[1]=lik_prod;

  return(pred_mat);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Given tree tables and model weights, obtain predicted probabilities for test data.
//' @export
// [[Rcpp::export]]

NumericMatrix get_test_prob_overall(NumericVector weights, int num_cats,
                                    NumericMatrix testdata, //NumericMatrix test_data,
                                    List treetable_list//, NumericMatrix tree_data
){

  //List overall_term_nodes_trees(overall_sum_trees.size());
  //List overall_term_obs_trees(overall_sum_trees.size());
  //List overall_predictions(overall_sum_trees.size());

  List overall_treetables(treetable_list.size());
  NumericVector overall_liks(treetable_list.size());


  arma::mat pred_mat_overall=arma::zeros<arma::mat>(testdata.nrow(),num_cats);


  for(int i=0;i<treetable_list.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = treetable_list[i];

    //NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      //List sum_tree=treetable_list[i];

      //save all info in list of list format the same as the trees.
      //List term_nodes_trees(sum_tree.size());
      //List term_obs_trees(sum_tree.size());
      //NumericMatrix predictions(num_obs,sum_tree.size());

      // for(int k=0;k<sum_tree.size();k++){
      //   NumericMatrix tree_table=sum_tree[k];
      //   List tree_info=get_termobs_test_data(test_data, tree_table) ;
      //   //NumericVector term_nodes=tree_info[0];
      //   //term_nodes_trees[k]=term_nodes;
      //   term_obs_trees[k]=tree_info;
      //   //umericVector term_preds=tree_info[2];
      //   //predictions(_,k)=term_preds;
      // }

      //Rcout << "Line 682. i== " << i << ". \n";

      arma::mat treeprob_output = get_test_probs(weights, num_cats,
                                                 testdata,
                                                 treetable_list[i]  );

      //Rcout << "Line 688. i== " << i << ". \n";

      double weighttemp = weights[i];
      //Rcout << "Line 691. i== " << i << ". \n";

      pred_mat_overall = pred_mat_overall + weighttemp*treeprob_output;
      //Rcout << "Line 694. i== " << i << ". \n";


      //overall_treetables[i]= treepred_output[0];
      //double templik = as<double>(treepred_output[1]);
      //overall_liks[i]= pow(templik,beta_pow);



      //overall_term_nodes_trees[i]=term_nodes_trees;
      //overall_term_obs_trees[i]= term_obs_trees;
      //overall_predictions[i]=predictions;
    }else{
      // NumericMatrix sum_tree=overall_sum_trees[i];
      // List tree_info=get_termobs_test_data(test_data, sum_tree) ;
      // //overall_term_nodes_trees[i]=tree_info[0];
      // List term_obs_trees(1);
      // term_obs_trees[0]=tree_info ;
      // //NumericVector term_preds=tree_info[2];
      // //NumericVector predictions=term_preds;
      // overall_term_obs_trees[i]= term_obs_trees;
      // //overall_predictions[i]=predictions;
      //

      //overall_treetables[i]=get_treepreds(original_y, num_cats, alpha_pars,
      //                                    originaldata,
      //                                    treetable_list[i]  );


      // List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
      //                                      originaldata,
      //                                      treetable_list[i]  );
      //
      // overall_treetables[i]= treepred_output[0];
      // double templik = as<double>(treepred_output[1]);
      // overall_liks[i]= pow(templik,beta_pow);


      //Rcout << "Line 732. i== " << i << ". \n";

      arma::mat treeprob_output = get_test_probs(weights, num_cats,
                                                 testdata,
                                                 treetable_list[i]  );

      //Rcout << "Line 738. i== " << i << ". \n";

      double weighttemp = weights[i];
      //Rcout << "Line 741. i== " << i << ". \n";
      //Rcout << "treeprob_output.n_rows" << treeprob_output.n_rows << ".\n";
      //Rcout << "treeprob_output.n_cols" << treeprob_output.n_cols << ".\n";


      pred_mat_overall = pred_mat_overall + weighttemp*treeprob_output;
      //Rcout << "Line 744. i== " << i << ". \n";
      //Rcout << "pred_mat_overall " << pred_mat_overall << ". \n";

    }
  }
  //List ret(1);
  //ret[0]=overall_term_nodes_trees;
  //ret[0]=overall_term_obs_trees;
  //ret[2]=overall_predictions;
  //return(overall_term_obs_trees);

  //return(overall_treetables);

  // overall_liks=overall_liks/sum(overall_liks);
  //
  // List ret(2);
  // ret[0]=overall_treetables;
  // ret[1]=overall_liks;
  // return(ret);

  return(wrap(pred_mat_overall));

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::depends(dqrng)]]
//// [[Rcpp::depends(BH)]]
//// [[Rcpp::depends(dqrng, BH, RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

//#include <dqrng.h>

//#include <boost/random/binomial_distribution.hpp>
//using binomial = boost::random::binomial_distribution<int>;
//' @title Draw a set of trees from the prior.
//' @export
// [[Rcpp::export]]
List draw_trees(double lambda, int num_trees, int seed, int num_split_vars, int num_cats ){

  //dqrng::dqRNGkind("Xoroshiro128+");
  //dqrng::dqset_seed(IntegerVector::create(seed));

  //use following with binomial?
  //dqrng::xoshiro256plus rng(seed);

  //std::vector<int> lambdavec = {lambda, 1-lambda};

  //typedef boost::mt19937 RNGType;
  //boost::random::uniform_int_distribution<> sample_splitvardist(1,num_split_vars);
  //boost::variate_generator< RNGType, boost::uniform_int<> >  sample_splitvars(rng, sample_splitvardist);

  //boost::random::uniform_real_distribution<double> b_unifdist(0,1);
  //boost::variate_generator< RNGType, boost::uniform_real<> >  b_unif_point(rng, b_unifdist);



  // std::random_device device;
  // std::mt19937 gen(device());

  //possibly use seed?
  //// std::mt19937 gen(seed);


  // std::bernoulli_distribution coin_flip(lambda);

  // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
  //std::uniform_real_distribution<> dis_cont_unif(0, 1);



  List table_list(num_trees);




  for(int j=0; j<num_trees;j++){

    //If parallelizing, define the distributinos before this loop
    //and use lrng and the following two lines
    //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
    //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


    //NumericVector treenodes_bin(0);
    //arma::uvec treenodes_bin(0);

    std::vector<int> treenodes_bin;


    int count_terminals = 0;
    int count_internals = 0;

    //int count_treebuild = 0;

    while(count_internals > (count_terminals -1)){

      //Also consider standard library and random header
      // std::random_device device;
      // std::mt19937 gen(device());
      // std::bernoulli_distribution coin_flip(lambda);
      // bool outcome = coin_flip(gen);


      //int tempdraw = coin_flip(gen);

      //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


      //int tempdraw = Rcpp::rbinom(1,lambda,1);
      int tempdraw = R::rbinom(1,lambda);
      treenodes_bin.push_back(tempdraw);

      //Rcout << "tempdraw = " << tempdraw << ".\n" ;

      //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;
      //need to update rng if use boost?
      //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));
      if(tempdraw==1){
        count_internals=count_internals+1;
      }else{
        count_terminals=count_terminals+1;
      }

    }//end of while loop creating parent vector treenodes_bin

    //Consider making this an armadillo vector
    //IntegerVector split_var_vec(treenodes_bin.size());
    //arma::uvec split_var_vec(treenodes_bin.size());
    std::vector<int> split_var_vec(treenodes_bin.size());

    //loop drawing splitting variables
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_var_vec[i] = -1;
      }else{
        // also consider the standard library function uniform_int_distribution
        // might need random header
        // This uses the Mersenne twister

        //Three lines below should probably be outside all the loops
        // std::random_device rd;
        // std::mt19937 engine(rd());
        // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
        //
        // split_var_vec[i] <- distsampvar(engine);

        // split_var_vec[i] <- distsampvar(gen);


        //consider using boost
        //might need to update rng
        //split_var_vec[i] <- sample_splitvars(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

        //not sure if this returns an integer or a vector?
        //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
        //could try
        split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
        //could also try RcppArmadillo::rmultinom

      }

    }// end of for-loop drawing split variables


    //Consider making this an armadillo vector
    //NumericVector split_point_vec(treenodes_bin.size());
    //arma::vec split_point_vec(treenodes_bin.size());
    std::vector<double> split_point_vec(treenodes_bin.size());


    //loop drawing splitting points
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_point_vec[i] = -1;
      }else{


        //////////////////////////////////////////////////////////
        //following function not reccommended
        //split_point_vec[i] = std::rand();
        //////////////////////////////////////////////////////////
        ////Standard library:
        ////This should probably be outside all the loops
        ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
        ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
        ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

        // split_point_vec[i] = dis_cont_unif(gen);

        //////////////////////////////////////////////////////////
        //from armadillo
        split_point_vec[i] = arma::randu();

        //////////////////////////////////////////////////////////
        //probably not adviseable for paralelization
        //From Rcpp
        split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

        //////////////////////////////////////////////////////////
        //consider using boost
        //might need to update rng
        //split_point_vec[i] <- b_unif_point(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

        //not sure if this returns an integer or a vector?





      }

    }// end of for-loop drawing split points


    //Create tree table matrix

    //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

    //Rcout << "Line 1037. \n";
    //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

    //initialize with zeros. Not sure if this is necessary
    arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),5+num_cats);
    //Rcout << "Line 1040. \n";


    //tree_table1(_,2) = wrap(split_var_vec);
    //tree_table1(_,3) = wrap(split_point_vec);
    //tree_table1(_,4) = wrap(treenodes_bin);

    //It might be more efficient to make everything an armadillo object initially
    // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
    arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
    arma::colvec split_point_vec_arma(split_point_vec);
    arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


    //Rcout << "Line 1054. \n";

    tree_table1.col(2) = split_var_vec_arma;
    tree_table1.col(3) = split_point_vec_arma;
    tree_table1.col(4) = treenodes_bin_arma;


    //Rcout << "Line 1061. j = " << j << ". \n";



    // Now start filling in left daughter and right daughter columns
    std::vector<int> rd_spaces;
    int prev_node = -1;

    for(unsigned int i=0; i<treenodes_bin.size();i++){
      //Rcout << "Line 1061. i = " << i << ". \n";
      if(prev_node==0){
        //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
        //Rcout << "Line 1073. j = " << j << ". \n";

        tree_table1(rd_spaces.back(), 1)=i+1;
        //Rcout << "Line 1076. j = " << j << ". \n";

        rd_spaces.pop_back();
      }
      if(treenodes_bin[i]==1){
        //Rcout << "Line 1081. j = " << j << ". \n";

        tree_table1(i,0) = i+2;
        rd_spaces.push_back(i);
        prev_node = 1;
        //Rcout << "Line 185. j = " << j << ". \n";

      }else{                  // These 2 lines unnecessary if begin with matrix of zeros
        //Rcout << "Line 1089. j = " << j << ". \n";
        tree_table1(i,0)=0 ;
        tree_table1(i,1) = 0 ;
        prev_node = 0;
        //Rcout << "Line 1093. j = " << j << ". \n";

      }
    }//
    //Rcout << "Line 1097. j = " << j << ". \n";

    table_list[j]=wrap(tree_table1);
    //Rcout << "Line 1100. j = " << j << ". \n";

  }//end of loop over all trees

  return(table_list);
}//end of function definition

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Safe-Bayesian Random Forest. Initial test function.
//' @export
// [[Rcpp::export]]

NumericMatrix sBayesRF(double lambda, int num_trees,
                       int seed, int num_cats,
                       NumericVector y, NumericMatrix original_datamat,
                       NumericVector alpha_parameters, double beta_par,
                       NumericMatrix test_datamat){

  int num_split_vars= original_datamat.ncol();
  NumericMatrix Data_transformed = cpptrans_cdf(original_datamat);
  NumericMatrix testdat_trans = cpptrans_cdf_test(original_datamat,test_datamat);

  //Rcout << "Line 1134 . \n";
  List table_list = draw_trees(lambda, num_trees, seed, num_split_vars, num_cats );
  //Rcout << "Line 1136 . \n";


  List tree_list_output = get_treelist(y, num_cats, alpha_parameters, beta_par,
                                       Data_transformed,
                                       table_list  );
  //Rcout << "Line 1141 . \n";

  NumericMatrix probmat = get_test_prob_overall(tree_list_output[1],num_cats,
                                                testdat_trans,
                                                tree_list_output[0]);
  //Rcout << "Line 1146 . \n";

  return(probmat);

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Safe-Bayesian Random Forest
//'
//' @description An implementation of the Safe-Bayesian Random Forest described by Quadrianto and Ghahramani (2015)
//' @param lambda A real number between 0 and 1 that determines the splitting probability in the prior (which is used as the importance sampler of tree models). Quadrianto and Ghahramani (2015) recommend a value less than 0.5 .
//' @param num_trees The number of trees to be sampled.
//' @param seed The seed for random number generation.
//' @param num_cats The number of possible values for the outcome variable.
//' @param y The training data vector of outcomes. This must be a vector of integers between 1 and num_cats.
//' @param original_datamat The original training data. Currently all variables must be continuous. The training data does not need to be transformed before being entered to this function.
//' @param alpha_parameters Vector of prior parameters.
//' @param beta_par The power to which the likelihood is to be raised. For BMA, set beta_par=1.
//' @param original_datamat The original test data. This matrix must have the same number of columns (variables) as the training data. Currently all variables must be continuous. The test data does not need to be transformed before being entered to this function.
//' @export
// [[Rcpp::export]]

NumericMatrix sBayesRF_onefunc(double lambda, int num_trees,
                               int seed, int num_cats,
                               NumericVector y, NumericMatrix original_datamat,
                               NumericVector alpha_parameters, double beta_par,
                               NumericMatrix test_datamat){

  int num_split_vars= original_datamat.ncol();


  ///////////////////////
  //NumericMatrix Data_transformed = cpptrans_cdf(original_datamat);
  NumericMatrix Data_transformed(original_datamat.nrow(), original_datamat.ncol());
  for(int i=0; i<original_datamat.ncol();i++){
    NumericVector samp= original_datamat(_,i);
    NumericVector sv(clone(samp));
    std::sort(sv.begin(), sv.end());
    double nobs = samp.size();
    NumericVector ans(nobs);
    for (int k = 0; k < samp.size(); ++k)
      ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
    //NumericVector ansnum = ans;
    Data_transformed(_,i) = (ans+1)/nobs;
  }



  /////////////////////////////////////
  //NumericMatrix testdat_trans = cpptrans_cdf_test(original_datamat,test_datamat);
  NumericMatrix testdat_trans(test_datamat.nrow(), test_datamat.ncol());
  for(int i=0; i<test_datamat.ncol();i++){
    NumericVector samp= test_datamat(_,i);
    NumericVector svtest = original_datamat(_,i);
    NumericVector sv(clone(svtest));
    std::sort(sv.begin(), sv.end());
    double nobs = samp.size();
    NumericVector ans(nobs);
    double nobsref = svtest.size();
    for (int k = 0; k < samp.size(); ++k){
      ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
    }
    //NumericVector ansnum = ans;
    testdat_trans(_,i) = (ans)/nobsref;
  }



  /////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////////////
  //List table_list = draw_trees(lambda, num_trees, seed, num_split_vars, num_cats );



  //dqrng::dqRNGkind("Xoroshiro128+");
  //dqrng::dqset_seed(IntegerVector::create(seed));

  //use following with binomial?
  //dqrng::xoshiro256plus rng(seed);

  //std::vector<int> lambdavec = {lambda, 1-lambda};

  //typedef boost::mt19937 RNGType;
  //boost::random::uniform_int_distribution<> sample_splitvardist(1,num_split_vars);
  //boost::variate_generator< RNGType, boost::uniform_int<> >  sample_splitvars(rng, sample_splitvardist);

  //boost::random::uniform_real_distribution<double> b_unifdist(0,1);
  //boost::variate_generator< RNGType, boost::uniform_real<> >  b_unif_point(rng, b_unifdist);



  // std::random_device device;
  // std::mt19937 gen(device());

  //possibly use seed?
  //// std::mt19937 gen(seed);


  // std::bernoulli_distribution coin_flip(lambda);

  // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
  //std::uniform_real_distribution<> dis_cont_unif(0, 1);


  arma::vec orig_y_arma= as<arma::vec>(y);
  arma::vec alpha_pars_arma= as<arma::vec>(alpha_parameters);

  arma::mat arma_orig_data(Data_transformed.begin(), Data_transformed.nrow(), Data_transformed.ncol(), false);
  arma::mat arma_test_data(testdat_trans.begin(), testdat_trans.nrow(), testdat_trans.ncol(), false);


  arma::mat pred_mat_overall=arma::zeros<arma::mat>(test_datamat.nrow(),num_cats);


  //List overall_treetables(num_trees);
  NumericVector overall_liks(num_trees);


  //overall_treetables[i]= wrap(tree_table1);
  //double templik = as<double>(treepred_output[1]);
  //overall_liks[i]= pow(lik_prod,beta_pow);



  for(int j=0; j<num_trees;j++){

    //If parallelizing, define the distributinos before this loop
    //and use lrng and the following two lines
    //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
    //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


    //NumericVector treenodes_bin(0);
    //arma::uvec treenodes_bin(0);

    std::vector<int> treenodes_bin;


    int count_terminals = 0;
    int count_internals = 0;

    //int count_treebuild = 0;

    while(count_internals > (count_terminals -1)){

      //Also consider standard library and random header
      // std::random_device device;
      // std::mt19937 gen(device());
      // std::bernoulli_distribution coin_flip(lambda);
      // bool outcome = coin_flip(gen);


      //int tempdraw = coin_flip(gen);

      //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


      //int tempdraw = Rcpp::rbinom(1,lambda,1);
      int tempdraw = R::rbinom(1,lambda);
      treenodes_bin.push_back(tempdraw);

      //Rcout << "tempdraw = " << tempdraw << ".\n" ;

      //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;
      //need to update rng if use boost?
      //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));
      if(tempdraw==1){
        count_internals=count_internals+1;
      }else{
        count_terminals=count_terminals+1;
      }

    }//end of while loop creating parent vector treenodes_bin

    //Consider making this an armadillo vector
    //IntegerVector split_var_vec(treenodes_bin.size());
    //arma::uvec split_var_vec(treenodes_bin.size());
    std::vector<int> split_var_vec(treenodes_bin.size());

    //loop drawing splitting variables
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_var_vec[i] = -1;
      }else{
        // also consider the standard library function uniform_int_distribution
        // might need random header
        // This uses the Mersenne twister

        //Three lines below should probably be outside all the loops
        // std::random_device rd;
        // std::mt19937 engine(rd());
        // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
        //
        // split_var_vec[i] <- distsampvar(engine);

        // split_var_vec[i] <- distsampvar(gen);


        //consider using boost
        //might need to update rng
        //split_var_vec[i] <- sample_splitvars(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

        //not sure if this returns an integer or a vector?
        //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
        //could try
        split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
        //could also try RcppArmadillo::rmultinom

      }

    }// end of for-loop drawing split variables


    //Consider making this an armadillo vector
    //NumericVector split_point_vec(treenodes_bin.size());
    //arma::vec split_point_vec(treenodes_bin.size());
    std::vector<double> split_point_vec(treenodes_bin.size());


    //loop drawing splitting points
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_point_vec[i] = -1;
      }else{


        //////////////////////////////////////////////////////////
        //following function not reccommended
        //split_point_vec[i] = std::rand();
        //////////////////////////////////////////////////////////
        ////Standard library:
        ////This should probably be outside all the loops
        ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
        ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
        ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

        // split_point_vec[i] = dis_cont_unif(gen);

        //////////////////////////////////////////////////////////
        //from armadillo
        //split_point_vec[i] = arma::randu();

        //////////////////////////////////////////////////////////
        //probably not adviseable for paralelization
        //From Rcpp
        split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

        //////////////////////////////////////////////////////////
        //consider using boost
        //might need to update rng
        //split_point_vec[i] <- b_unif_point(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

        //not sure if this returns an integer or a vector?





      }

    }// end of for-loop drawing split points


    //Create tree table matrix

    //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

    //Rcout << "Line 1037. \n";
    //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

    //initialize with zeros. Not sure if this is necessary
    arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),5+num_cats);
    //Rcout << "Line 1040. \n";


    //tree_table1(_,2) = wrap(split_var_vec);
    //tree_table1(_,3) = wrap(split_point_vec);
    //tree_table1(_,4) = wrap(treenodes_bin);

    //It might be more efficient to make everything an armadillo object initially
    // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
    arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
    arma::colvec split_point_vec_arma(split_point_vec);
    arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


    //Rcout << "Line 1054. \n";

    tree_table1.col(2) = split_var_vec_arma;
    tree_table1.col(3) = split_point_vec_arma;
    tree_table1.col(4) = treenodes_bin_arma;


    //Rcout << "Line 1061. j = " << j << ". \n";



    // Now start filling in left daughter and right daughter columns
    std::vector<int> rd_spaces;
    int prev_node = -1;

    for(unsigned int i=0; i<treenodes_bin.size();i++){
      //Rcout << "Line 1061. i = " << i << ". \n";
      if(prev_node==0){
        //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
        //Rcout << "Line 1073. j = " << j << ". \n";

        tree_table1(rd_spaces.back(), 1)=i+1;
        //Rcout << "Line 1076. j = " << j << ". \n";

        rd_spaces.pop_back();
      }
      if(treenodes_bin[i]==1){
        //Rcout << "Line 1081. j = " << j << ". \n";

        tree_table1(i,0) = i+2;
        rd_spaces.push_back(i);
        prev_node = 1;
        //Rcout << "Line 185. j = " << j << ". \n";

      }else{                  // These 2 lines unnecessary if begin with matrix of zeros
        //Rcout << "Line 1089. j = " << j << ". \n";
        tree_table1(i,0)=0 ;
        tree_table1(i,1) = 0 ;
        prev_node = 0;
        //Rcout << "Line 1093. j = " << j << ". \n";

      }
    }//
    //Rcout << "Line 1097. j = " << j << ". \n";





    //List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
    //                                     originaldata,
    //                                     treetable_list[i]  );


    //use armadillo object tree_table1

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////



    double lik_prod=1;
    double alph_prod=1;
    for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
      alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
    }
    double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
    double alph_term=gam_alph_sum/alph_prod;

    //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
    //arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


    //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
    //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

    //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

    //NumericVector terminal_nodes=find_term_nodes(treetable);

    //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

    //arma::vec colmat=arma_tree.col(4);
    //arma::uvec term_nodes=arma::find(colmat==-1);

    //arma::vec colmat=arma_tree.col(2);
    //arma::uvec term_nodes=arma::find(colmat==0);

    arma::vec colmat=tree_table1.col(4);
    arma::uvec term_nodes=arma::find(colmat==0);

    term_nodes=term_nodes+1;

    NumericVector terminal_nodes= wrap(term_nodes);




    //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
    //NumericVector tree_predictions;

    //now for each internal node find the observations that belong to the terminal nodes

    //NumericVector predictions(test_data.nrow());
    //List term_obs(terminal_nodes.size());
    if(terminal_nodes.size()==1){
      //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
      //predictions=rep(nodemean,test_data.nrow());
      //Rcout << "Line 67 .\n";

      //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
      //term_obs[0]= temp_obsvec;
      double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

      double num_prod=1;
      double num_sum=0;
      //Rcout << "Line 129.\n";

      for(int k=0; k<num_cats; k++){
        //assuming categories of y are from 1 to num_cats
        arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
        double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        tree_table1(0,5+k)= m_plus_alph/denom_temp ;

        //for likelihood calculation
        num_prod=num_prod*tgamma(m_plus_alph);
        num_sum=num_sum +m_plus_alph ;
      }

      lik_prod= alph_term*num_prod/tgamma(num_sum);

    }
    else{
      for(int i=0;i<terminal_nodes.size();i++){
        //arma::mat subdata=testd;
        int curr_term=terminal_nodes[i];

        int row_index;
        int term_node=terminal_nodes[i];
        //Rcout << "Line 152.\n";


        //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
        //Why should the ro index be different for a right daughter?
        //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
        row_index=0;

        // if(curr_term % 2==0){
        //   //term node is left daughter
        //   row_index=terminal_nodes[i];
        // }else{
        //   //term node is right daughter
        //   row_index=terminal_nodes[i]-1;
        // }




        //save the left and right node data into arma uvec

        //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
        //arma::vec left_nodes=arma_tree.col(0);
        //arma::vec right_nodes=arma_tree.col(1);

        arma::vec left_nodes=tree_table1.col(0);
        arma::vec right_nodes=tree_table1.col(1);



        arma::mat node_split_mat;
        node_split_mat.set_size(0,3);
        //Rcout << "Line 182. i = " << i << " .\n";

        while(row_index!=1){
          //for each terminal node work backwards and see if the parent node was a left or right node
          //append split info to a matrix
          int rd=0;
          arma::uvec parent_node=arma::find(left_nodes == term_node);

          if(parent_node.size()==0){
            parent_node=arma::find(right_nodes == term_node);
            rd=1;
          }

          //want to cout parent node and append to node_split_mat

          node_split_mat.insert_rows(0,1);

          //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
          //node_split_mat(0,0)=treetable(parent_node[0],2);
          //node_split_mat(0,1)=treetable(parent_node[0],3);

          //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
          //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

          node_split_mat(0,0)=tree_table1(parent_node[0],2);
          node_split_mat(0,1)=tree_table1(parent_node[0],3);

          node_split_mat(0,2)=rd;
          row_index=parent_node[0]+1;
          term_node=parent_node[0]+1;
        }

        //once we have the split info, loop through rows and find the subset indexes for that terminal node!
        //then fill in the predicted value for that tree
        //double prediction = tree_data(term_node,5);
        arma::uvec pred_indices;
        int split= node_split_mat(0,0)-1;

        //Rcout << "Line 224.\n";
        //Rcout << "split = " << split << ".\n";
        //arma::vec tempvec = testd.col(split);
        arma::vec tempvec = arma_orig_data.col(split);
        //Rcout << "Line 227.\n";


        double temp_split = node_split_mat(0,1);

        if(node_split_mat(0,2)==0){
          pred_indices = arma::find(tempvec <= temp_split);
        }else{
          pred_indices = arma::find(tempvec > temp_split);
        }
        //Rcout << "Line 236.\n";

        arma::uvec temp_pred_indices;

        //arma::vec data_subset = testd.col(split);
        arma::vec data_subset = arma_orig_data.col(split);

        data_subset=data_subset.elem(pred_indices);

        //now loop through each row of node_split_mat
        int n=node_split_mat.n_rows;
        //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
        //Rcout << "Line 248.\n";

        for(int j=1;j<n;j++){
          int curr_sv=node_split_mat(j,0);
          double split_p = node_split_mat(j,1);

          //data_subset = testd.col(curr_sv-1);
          //Rcout << "Line 255.\n";
          //Rcout << "curr_sv = " << curr_sv << ".\n";
          data_subset = arma_orig_data.col(curr_sv-1);
          //Rcout << "Line 258.\n";

          data_subset=data_subset.elem(pred_indices);

          if(node_split_mat(j,2)==0){
            //split is to the left
            temp_pred_indices=arma::find(data_subset <= split_p);
          }else{
            //split is to the right
            temp_pred_indices=arma::find(data_subset > split_p);
          }
          pred_indices=pred_indices.elem(temp_pred_indices);

          if(pred_indices.size()==0){
            continue;
          }

        }
        //Rcout << "Line 199. i = " << i <<  ".\n";

        //double nodemean=tree_data(terminal_nodes[i]-1,5);
        //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
        //predictions[predind]= nodemean;
        //term_obs[i]=predind;

        double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
        //Rcout << "Line 207. predind = " << predind <<  ".\n";
        //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
        // << "Line 207. term_node = " << term_node <<  ".\n";

        double num_prod=1;
        double num_sum=0;

        for(int k=0; k<num_cats; k++){
          //assuming categories of y are from 1 to num_cats
          arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
          double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);

          tree_table1(curr_term-1,5+k)= m_plus_alph/denom_temp ;

          num_prod=num_prod*tgamma(m_plus_alph);
          num_sum=num_sum +m_plus_alph ;
        }


        lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);
        //Rcout << "Line 297.\n";


      }
      //Rcout << "Line 301.\n";

    }
    //List ret(1);
    //ret[0] = term_obs;

    //ret[0] = terminal_nodes;
    //ret[1] = term_obs;
    //ret[2] = predictions;
    //return(term_obs);
    //Rcout << "Line 309";

    //return(wrap(arma_tree_table));

    //List ret(2);
    //ret[0]=wrap(arma_tree_table);
    //ret[1]=lik_prod;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////






    //overall_treetables[j]= wrap(tree_table1);


    //double templik = as<double>(treepred_output[1]);

    double templik = pow(lik_prod,beta_par);
    overall_liks[j]= templik;






    //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
    //arma::mat arma_test_data(testdata.begin(), testdata.nrow(), testdata.ncol(), false);


    //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
    //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

    //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

    //NumericVector terminal_nodes=find_term_nodes(treetable);
    //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
    //NumericVector tree_predictions;

    //now for each internal node find the observations that belong to the terminal nodes

    //NumericVector predictions(test_data.nrow());

    arma::mat pred_mat(test_datamat.nrow(),num_cats);
    //arma::vec filled_in(testdata.nrow());


    //List term_obs(terminal_nodes.size());
    if(terminal_nodes.size()==1){

      //Rcout << "Line 422. \n";


      pred_mat=repmat(tree_table1(0,arma::span(5,5+num_cats-1)),test_datamat.nrow(),1);


      //Rcout << "Line 424. \n";


      // for(int k=0; k<num_cats; k++){
      // pred_mat(_,k)=rep(treetable(0,5+k),testdata.nrow());
      // }
      //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
      //predictions=rep(nodemean,test_data.nrow());
      //Rcout << "Line 67 .\n";

      //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
      //term_obs[0]= temp_obsvec;
      // double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);
      //
      // double num_prod=1;
      // double num_sum=0;

      // for(int k=0; k<num_cats; k++){
      //   //assuming categories of y are from 1 to num_cats
      //   arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
      //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
      //   arma_tree_table(0,5+k)= m_plus_alph/denom_temp ;
      //
      //   //for likelihood calculation
      //   num_prod=num_prod*tgamma(m_plus_alph);
      //   num_sum=num_sum +m_plus_alph ;
      // }
      //
      // lik_prod= alph_term*num_prod/tgamma(num_sum);
      //
    }
    else{
      for(int i=0;i<terminal_nodes.size();i++){
        //arma::mat subdata=testd;
        int curr_term=terminal_nodes[i];

        int row_index;
        int term_node=terminal_nodes[i];


        //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
        //Why should the ro index be different for a right daughter?
        //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
        row_index=0;

        // if(curr_term % 2==0){
        //   //term node is left daughter
        //   row_index=terminal_nodes[i];
        // }else{
        //   //term node is right daughter
        //   row_index=terminal_nodes[i]-1;
        // }








        //save the left and right node data into arma uvec

        //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
        //arma::vec left_nodes=arma_tree.col(0);
        //arma::vec right_nodes=arma_tree.col(1);

        arma::vec left_nodes=tree_table1.col(0);
        arma::vec right_nodes=tree_table1.col(1);



        arma::mat node_split_mat;
        node_split_mat.set_size(0,3);
        //Rcout << "Line 124. i = " << i << " .\n";

        while(row_index!=1){
          //for each terminal node work backwards and see if the parent node was a left or right node
          //append split info to a matrix
          int rd=0;
          arma::uvec parent_node=arma::find(left_nodes == term_node);

          if(parent_node.size()==0){
            parent_node=arma::find(right_nodes == term_node);
            rd=1;
          }

          //want to cout parent node and append to node_split_mat

          node_split_mat.insert_rows(0,1);

          //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
          //node_split_mat(0,0)=treetable(parent_node[0],2);
          //node_split_mat(0,1)=treetable(parent_node[0],3);

          //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
          //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

          node_split_mat(0,0)=tree_table1(parent_node[0],2);
          node_split_mat(0,1)=tree_table1(parent_node[0],3);

          node_split_mat(0,2)=rd;
          row_index=parent_node[0]+1;
          term_node=parent_node[0]+1;
        }

        //once we have the split info, loop through rows and find the subset indexes for that terminal node!
        //then fill in the predicted value for that tree
        //double prediction = tree_data(term_node,5);
        arma::uvec pred_indices;
        int split= node_split_mat(0,0)-1;

        //arma::vec tempvec = testd.col(split);
        arma::vec tempvec = arma_test_data.col(split);


        double temp_split = node_split_mat(0,1);

        if(node_split_mat(0,2)==0){
          pred_indices = arma::find(tempvec <= temp_split);
        }else{
          pred_indices = arma::find(tempvec > temp_split);
        }

        arma::uvec temp_pred_indices;

        //arma::vec data_subset = testd.col(split);
        arma::vec data_subset = arma_test_data.col(split);

        data_subset=data_subset.elem(pred_indices);

        //now loop through each row of node_split_mat
        int n=node_split_mat.n_rows;
        //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
        //Rcout << "Line 174. node_split_mat= " << node_split_mat << ". n = " << n << ".\n";


        for(int j=1;j<n;j++){
          int curr_sv=node_split_mat(j,0);
          double split_p = node_split_mat(j,1);

          //data_subset = testd.col(curr_sv-1);
          data_subset = arma_test_data.col(curr_sv-1);

          data_subset=data_subset.elem(pred_indices);

          if(node_split_mat(j,2)==0){
            //split is to the left
            temp_pred_indices=arma::find(data_subset <= split_p);
          }else{
            //split is to the right
            temp_pred_indices=arma::find(data_subset > split_p);
          }
          pred_indices=pred_indices.elem(temp_pred_indices);

          if(pred_indices.size()==0){
            continue;
          }

        }
        //Rcout << "Line 199. i = " << i <<  ".\n";

        //double nodemean=tree_data(terminal_nodes[i]-1,5);
        //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
        //predictions[predind]= nodemean;
        //term_obs[i]=predind;

        //Rcout << "Line 635. \n";
        //Rcout << "pred_indices = " << pred_indices << ".\n";

        //pred_mat.rows(pred_indices)=arma::repmat(arma_tree_table(curr_term-1,arma::span(5,5+num_cats-1)),pred_indices.n_elem,1);
        pred_mat.each_row(pred_indices)=tree_table1(curr_term-1,arma::span(5,4+num_cats));



        //Rcout << "Line 588. \n";

        // for(int k=0; k<num_cats; k++){
        //   pred_mat(predind,k)=rep(treetable(curr_term-1,5+k),predind.size());
        // }


        // double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
        // //Rcout << "Line 207. predind = " << predind <<  ".\n";
        // //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
        // // << "Line 207. term_node = " << term_node <<  ".\n";
        //
        // double num_prod=1;
        // double num_sum=0;
        //
        // for(int k=0; k<num_cats; k++){
        //   //assuming categories of y are from 1 to num_cats
        //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
        //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        //
        //   arma_tree_table(curr_term-1,5+k)= m_plus_alph/denom_temp ;
        //
        //   num_prod=num_prod*tgamma(m_plus_alph);
        //   num_sum=num_sum +m_plus_alph ;
        // }


        //lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);


      }
    }






    //THIS SHOULD BE DIFFERENT IF THE CODE IS TO BE PARALLELIZED
    //EACH THREAD SHOULD OUTPUT ITS OWN MATRIX AND SUM OF LIKELIHOODS
    //THEN ADD THE MATRICES TOGETHER AND DIVIDE BY THE TOTAL SUM OF LIKELIHOODS
    //OR JUST SAVE ALL MATRICES TO ONE LIST

    pred_mat_overall = pred_mat_overall + templik*pred_mat;




    //arma::mat treeprob_output = get_test_probs(weights, num_cats,
    //                                           testdata,
    //                                           treetable_list[i]  );

    //Rcout << "Line 688. i== " << i << ". \n";

    //double weighttemp = weights[i];
    //Rcout << "Line 691. i== " << i << ". \n";

    //pred_mat_overall = pred_mat_overall + weighttemp*treeprob_output;



  }//end of loop over all trees




  ///////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////



  double sumlik_total= sum(overall_liks);
  pred_mat_overall=pred_mat_overall*(1/sumlik_total);
  //Rcout << "Line 1141 . \n";
  //Rcout << "Line 1146 . \n";

  return(wrap(pred_mat_overall));

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Safe-Bayesian Random Forest in C++
//'
//' @description An implementation of the Safe-Bayesian Random Forest described by Quadrianto and Ghahramani (2015)
//' @param lambda A real number between 0 and 1 that determines the splitting probability in the prior (which is used as the importance sampler of tree models). Quadrianto and Ghahramani (2015) recommend a value less than 0.5 .
//' @param num_trees The number of trees to be sampled.
//' @param seed The seed for random number generation.
//' @param num_cats The number of possible values for the outcome variable.
//' @param y The training data vector of outcomes. This must be a vector of integers between 1 and num_cats.
//' @param original_datamat The original training data. Currently all variables must be continuous. The training data does not need to be transformed before being entered to this function.
//' @param alpha_parameters Vector of prior parameters.
//' @param beta_par The power to which the likelihood is to be raised. For BMA, set beta_par=1.
//' @param original_datamat The original test data. This matrix must have the same number of columns (variables) as the training data. Currently all variables must be continuous. The test data does not need to be transformed before being entered to this function.
//' @export
// [[Rcpp::export]]

NumericMatrix sBayesRF_onefunc_arma(double lambda, int num_trees,
                                    int seed, int num_cats,
                                    NumericVector y, NumericMatrix original_datamat,
                                    NumericVector alpha_parameters, double beta_par,
                                    NumericMatrix test_datamat){

  int num_split_vars= original_datamat.ncol();
  arma::mat data_arma= as<arma::mat>(original_datamat);
  arma::mat testdata_arma= as<arma::mat>(test_datamat);
  arma::vec orig_y_arma= as<arma::vec>(y);
  arma::vec alpha_pars_arma= as<arma::vec>(alpha_parameters);

  ///////////////////////
  //NumericMatrix Data_transformed = cpptrans_cdf(original_datamat);
  // NumericMatrix Data_transformed(original_datamat.nrow(), original_datamat.ncol());
  // for(int i=0; i<original_datamat.ncol();i++){
  //   NumericVector samp= original_datamat(_,i);
  //   NumericVector sv(clone(samp));
  //   std::sort(sv.begin(), sv.end());
  //   double nobs = samp.size();
  //   NumericVector ans(nobs);
  //   for (int k = 0; k < samp.size(); ++k)
  //     ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
  //   //NumericVector ansnum = ans;
  //   Data_transformed(_,i) = (ans+1)/nobs;
  // }



  //arma::mat arma_orig_data(Data_transformed.begin(), Data_transformed.nrow(), Data_transformed.ncol(), false);



  //NumericMatrix transformedData(originaldata.nrow(), originaldata.ncol());

  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  arma::mat arma_orig_data(data_arma.n_rows,data_arma.n_cols);
  for(unsigned int k=0; k<data_arma.n_cols;k++){
    arma::vec samp= data_arma.col(k);
    arma::vec sv=arma::sort(samp);
    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      while (sv(j) < ssampi && j < sv.size()) ++j;
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }
    arma_orig_data.col(k)=(ans+1)/nobs;
  }





  /////////////////////////////////////
  // NumericMatrix testdat_trans = cpptrans_cdf_test(original_datamat,test_datamat);
  // //NumericMatrix testdat_trans(test_datamat.nrow(), test_datamat.ncol());
  // for(int i=0; i<test_datamat.ncol();i++){
  //   NumericVector samp= test_datamat(_,i);
  //   NumericVector svtest = original_datamat(_,i);
  //   NumericVector sv(clone(svtest));
  //   std::sort(sv.begin(), sv.end());
  //   double nobs = samp.size();
  //   NumericVector ans(nobs);
  //   double nobsref = svtest.size();
  //   for (int k = 0; k < samp.size(); ++k){
  //     ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
  //   }
  //   //NumericVector ansnum = ans;
  //   testdat_trans(_,i) = (ans)/nobsref;
  // }





  //NumericMatrix transformedData(originaldata.nrow(), originaldata.ncol());
  //arma::mat data_arma= as<arma::mat>(originaldata);

  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  arma::mat arma_test_data(testdata_arma.n_rows,testdata_arma.n_cols);
  for(unsigned int k=0; k<data_arma.n_cols;k++){
    arma::vec ref= data_arma.col(k);
    arma::vec samp= testdata_arma.col(k);

    arma::vec sv=arma::sort(samp);
    arma::vec sref=arma::sort(ref);

    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    double nobsref = ref.n_elem;

    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      if(j+1>sref.size()){
      }else{
        while (sref(j) < ssampi && j < sref.size()){
          ++j;
          if(j==sref.size()) break;
        }
      }
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }

    arma_test_data.col(k)=(ans)/nobsref;

  }







  /////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////////////
  //List table_list = draw_trees(lambda, num_trees, seed, num_split_vars, num_cats );



  //dqrng::dqRNGkind("Xoroshiro128+");
  //dqrng::dqset_seed(IntegerVector::create(seed));

  //use following with binomial?
  //dqrng::xoshiro256plus rng(seed);

  //std::vector<int> lambdavec = {lambda, 1-lambda};

  //typedef boost::mt19937 RNGType;
  //boost::random::uniform_int_distribution<> sample_splitvardist(1,num_split_vars);
  //boost::variate_generator< RNGType, boost::uniform_int<> >  sample_splitvars(rng, sample_splitvardist);

  //boost::random::uniform_real_distribution<double> b_unifdist(0,1);
  //boost::variate_generator< RNGType, boost::uniform_real<> >  b_unif_point(rng, b_unifdist);



  std::random_device device;
  std::mt19937 gen(device());

  //possibly use seed?
  //// std::mt19937 gen(seed);


  std::bernoulli_distribution coin_flip(lambda);

  std::uniform_int_distribution<> distsampvar(1, num_split_vars);
  std::uniform_real_distribution<> dis_cont_unif(0, 1);




  //arma::mat arma_test_data(testdat_trans.begin(), testdat_trans.nrow(), testdat_trans.ncol(), false);


  arma::mat pred_mat_overall=arma::zeros<arma::mat>(arma_test_data.n_rows,num_cats);


  //List overall_treetables(num_trees);
  arma::vec overall_liks(num_trees);


  //overall_treetables[i]= wrap(tree_table1);
  //double templik = as<double>(treepred_output[1]);
  //overall_liks[i]= pow(lik_prod,beta_pow);



  for(int j=0; j<num_trees;j++){

    //If parallelizing, define the distributinos before this loop
    //and use lrng and the following two lines
    //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
    //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


    //NumericVector treenodes_bin(0);
    //arma::uvec treenodes_bin(0);

    std::vector<int> treenodes_bin;


    int count_terminals = 0;
    int count_internals = 0;

    //int count_treebuild = 0;

    while(count_internals > (count_terminals -1)){

      //Also consider standard library and random header
      // std::random_device device;
      // std::mt19937 gen(device());
      // std::bernoulli_distribution coin_flip(lambda);
      // bool outcome = coin_flip(gen);


      int tempdraw = coin_flip(gen);

      //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


      //int tempdraw = Rcpp::rbinom(1,lambda,1);
      //int tempdraw = R::rbinom(1,lambda);
      treenodes_bin.push_back(tempdraw);

      //Rcout << "tempdraw = " << tempdraw << ".\n" ;

      //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;
      //need to update rng if use boost?
      //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));
      if(tempdraw==1){
        count_internals=count_internals+1;
      }else{
        count_terminals=count_terminals+1;
      }

    }//end of while loop creating parent vector treenodes_bin

    //Consider making this an armadillo vector
    //IntegerVector split_var_vec(treenodes_bin.size());
    //arma::uvec split_var_vec(treenodes_bin.size());
    std::vector<int> split_var_vec(treenodes_bin.size());

    //loop drawing splitting variables
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_var_vec[i] = -1;
      }else{
        // also consider the standard library function uniform_int_distribution
        // might need random header
        // This uses the Mersenne twister

        //Three lines below should probably be outside all the loops
        // std::random_device rd;
        // std::mt19937 engine(rd());
        // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
        //
        // split_var_vec[i] = distsampvar(engine);

        split_var_vec[i] = distsampvar(gen);


        //consider using boost
        //might need to update rng
        //split_var_vec[i] <- sample_splitvars(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

        //not sure if this returns an integer or a vector?
        //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
        //could try
        //split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
        //could also try RcppArmadillo::rmultinom

      }

    }// end of for-loop drawing split variables


    //Consider making this an armadillo vector
    //NumericVector split_point_vec(treenodes_bin.size());
    //arma::vec split_point_vec(treenodes_bin.size());
    std::vector<double> split_point_vec(treenodes_bin.size());


    //loop drawing splitting points
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_point_vec[i] = -1;
      }else{


        //////////////////////////////////////////////////////////
        //following function not reccommended
        //split_point_vec[i] = std::rand();
        //////////////////////////////////////////////////////////
        ////Standard library:
        ////This should probably be outside all the loops
        ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
        ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
        ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

        split_point_vec[i] = dis_cont_unif(gen);

        //////////////////////////////////////////////////////////
        //from armadillo
        //split_point_vec[i] = arma::randu();

        //////////////////////////////////////////////////////////
        //probably not adviseable for paralelization
        //From Rcpp
        //split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

        //////////////////////////////////////////////////////////
        //consider using boost
        //might need to update rng
        //split_point_vec[i] <- b_unif_point(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

        //not sure if this returns an integer or a vector?





      }

    }// end of for-loop drawing split points


    //Create tree table matrix

    //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

    //Rcout << "Line 1037. \n";
    //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

    //initialize with zeros. Not sure if this is necessary
    arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),5+num_cats);
    //Rcout << "Line 1040. \n";


    //tree_table1(_,2) = wrap(split_var_vec);
    //tree_table1(_,3) = wrap(split_point_vec);
    //tree_table1(_,4) = wrap(treenodes_bin);

    //It might be more efficient to make everything an armadillo object initially
    // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
    arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
    arma::colvec split_point_vec_arma(split_point_vec);
    arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


    //Rcout << "Line 1054. \n";

    tree_table1.col(2) = split_var_vec_arma;
    tree_table1.col(3) = split_point_vec_arma;
    tree_table1.col(4) = treenodes_bin_arma;


    //Rcout << "Line 1061. j = " << j << ". \n";



    // Now start filling in left daughter and right daughter columns
    std::vector<int> rd_spaces;
    int prev_node = -1;

    for(unsigned int i=0; i<treenodes_bin.size();i++){
      //Rcout << "Line 1061. i = " << i << ". \n";
      if(prev_node==0){
        //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
        //Rcout << "Line 1073. j = " << j << ". \n";

        tree_table1(rd_spaces.back(), 1)=i+1;
        //Rcout << "Line 1076. j = " << j << ". \n";

        rd_spaces.pop_back();
      }
      if(treenodes_bin[i]==1){
        //Rcout << "Line 1081. j = " << j << ". \n";

        tree_table1(i,0) = i+2;
        rd_spaces.push_back(i);
        prev_node = 1;
        //Rcout << "Line 185. j = " << j << ". \n";

      }else{                  // These 2 lines unnecessary if begin with matrix of zeros
        //Rcout << "Line 1089. j = " << j << ". \n";
        tree_table1(i,0)=0 ;
        tree_table1(i,1) = 0 ;
        prev_node = 0;
        //Rcout << "Line 1093. j = " << j << ". \n";

      }
    }//
    //Rcout << "Line 1097. j = " << j << ". \n";





    //List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
    //                                     originaldata,
    //                                     treetable_list[i]  );


    //use armadillo object tree_table1

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////



    double lik_prod=1;
    double alph_prod=1;
    for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
      alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
    }
    double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
    double alph_term=gam_alph_sum/alph_prod;

    //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
    //arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


    //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
    //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

    //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

    //NumericVector terminal_nodes=find_term_nodes(treetable);

    //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

    //arma::vec colmat=arma_tree.col(4);
    //arma::uvec term_nodes=arma::find(colmat==-1);

    //arma::vec colmat=arma_tree.col(2);
    //arma::uvec term_nodes=arma::find(colmat==0);

    arma::vec colmat=tree_table1.col(4);
    arma::uvec term_nodes=arma::find(colmat==0);

    term_nodes=term_nodes+1;

    //NumericVector terminal_nodes= wrap(term_nodes);




    //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
    //NumericVector tree_predictions;

    //now for each internal node find the observations that belong to the terminal nodes

    //NumericVector predictions(test_data.nrow());
    //List term_obs(term_nodes.n_elem);
    if(term_nodes.n_elem==1){
      //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
      //predictions=rep(nodemean,test_data.nrow());
      //Rcout << "Line 67 .\n";

      //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
      //term_obs[0]= temp_obsvec;
      double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

      double num_prod=1;
      double num_sum=0;
      //Rcout << "Line 129.\n";

      for(int k=0; k<num_cats; k++){
        //assuming categories of y are from 1 to num_cats
        arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
        double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        tree_table1(0,5+k)= m_plus_alph/denom_temp ;

        //for likelihood calculation
        num_prod=num_prod*tgamma(m_plus_alph);
        num_sum=num_sum +m_plus_alph ;
      }

      lik_prod= alph_term*num_prod/tgamma(num_sum);

    }
    else{
      for(unsigned int i=0;i<term_nodes.n_elem;i++){
        //arma::mat subdata=testd;
        int curr_term=term_nodes(i);

        int row_index;
        int term_node=term_nodes(i);
        //Rcout << "Line 152.\n";


        //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
        //Why should the ro index be different for a right daughter?
        //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
        row_index=0;

        // if(curr_term % 2==0){
        //   //term node is left daughter
        //   row_index=terminal_nodes[i];
        // }else{
        //   //term node is right daughter
        //   row_index=terminal_nodes[i]-1;
        // }




        //save the left and right node data into arma uvec

        //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
        //arma::vec left_nodes=arma_tree.col(0);
        //arma::vec right_nodes=arma_tree.col(1);

        arma::vec left_nodes=tree_table1.col(0);
        arma::vec right_nodes=tree_table1.col(1);



        arma::mat node_split_mat;
        node_split_mat.set_size(0,3);
        //Rcout << "Line 182. i = " << i << " .\n";

        while(row_index!=1){
          //for each terminal node work backwards and see if the parent node was a left or right node
          //append split info to a matrix
          int rd=0;
          arma::uvec parent_node=arma::find(left_nodes == term_node);

          if(parent_node.size()==0){
            parent_node=arma::find(right_nodes == term_node);
            rd=1;
          }

          //want to cout parent node and append to node_split_mat

          node_split_mat.insert_rows(0,1);

          //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
          //node_split_mat(0,0)=treetable(parent_node[0],2);
          //node_split_mat(0,1)=treetable(parent_node[0],3);

          //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
          //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

          node_split_mat(0,0)=tree_table1(parent_node(0),2);
          node_split_mat(0,1)=tree_table1(parent_node(0),3);

          node_split_mat(0,2)=rd;
          row_index=parent_node(0)+1;
          term_node=parent_node(0)+1;
        }

        //once we have the split info, loop through rows and find the subset indexes for that terminal node!
        //then fill in the predicted value for that tree
        //double prediction = tree_data(term_node,5);
        arma::uvec pred_indices;
        int split= node_split_mat(0,0)-1;

        //Rcout << "Line 224.\n";
        //Rcout << "split = " << split << ".\n";
        //arma::vec tempvec = testd.col(split);
        arma::vec tempvec = arma_orig_data.col(split);
        //Rcout << "Line 227.\n";


        double temp_split = node_split_mat(0,1);

        if(node_split_mat(0,2)==0){
          pred_indices = arma::find(tempvec <= temp_split);
        }else{
          pred_indices = arma::find(tempvec > temp_split);
        }
        //Rcout << "Line 236.\n";

        arma::uvec temp_pred_indices;

        //arma::vec data_subset = testd.col(split);
        arma::vec data_subset = arma_orig_data.col(split);

        data_subset=data_subset.elem(pred_indices);

        //now loop through each row of node_split_mat
        int n=node_split_mat.n_rows;
        //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
        //Rcout << "Line 248.\n";

        for(int j=1;j<n;j++){
          int curr_sv=node_split_mat(j,0);
          double split_p = node_split_mat(j,1);

          //data_subset = testd.col(curr_sv-1);
          //Rcout << "Line 255.\n";
          //Rcout << "curr_sv = " << curr_sv << ".\n";
          data_subset = arma_orig_data.col(curr_sv-1);
          //Rcout << "Line 258.\n";

          data_subset=data_subset.elem(pred_indices);

          if(node_split_mat(j,2)==0){
            //split is to the left
            temp_pred_indices=arma::find(data_subset <= split_p);
          }else{
            //split is to the right
            temp_pred_indices=arma::find(data_subset > split_p);
          }
          pred_indices=pred_indices.elem(temp_pred_indices);

          if(pred_indices.size()==0){
            continue;
          }

        }
        //Rcout << "Line 199. i = " << i <<  ".\n";

        //double nodemean=tree_data(terminal_nodes[i]-1,5);
        //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
        //predictions[predind]= nodemean;
        //term_obs[i]=predind;

        double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
        //Rcout << "Line 207. predind = " << predind <<  ".\n";
        //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
        // << "Line 207. term_node = " << term_node <<  ".\n";

        double num_prod=1;
        double num_sum=0;

        for(int k=0; k<num_cats; k++){
          //assuming categories of y are from 1 to num_cats
          arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
          double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);

          tree_table1(curr_term-1,5+k)= m_plus_alph/denom_temp ;

          num_prod=num_prod*tgamma(m_plus_alph);
          num_sum=num_sum +m_plus_alph ;
        }


        lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);
        //Rcout << "Line 297.\n";


      }
      //Rcout << "Line 301.\n";

    }
    //List ret(1);
    //ret[0] = term_obs;

    //ret[0] = terminal_nodes;
    //ret[1] = term_obs;
    //ret[2] = predictions;
    //return(term_obs);
    //Rcout << "Line 309";

    //return(wrap(arma_tree_table));

    //List ret(2);
    //ret[0]=wrap(arma_tree_table);
    //ret[1]=lik_prod;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////






    //overall_treetables[j]= wrap(tree_table1);


    //double templik = as<double>(treepred_output[1]);

    double templik = pow(lik_prod,beta_par);
    overall_liks(j)= templik;






    //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
    //arma::mat arma_test_data(testdata.begin(), testdata.nrow(), testdata.ncol(), false);


    //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
    //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

    //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

    //NumericVector terminal_nodes=find_term_nodes(treetable);
    //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
    //NumericVector tree_predictions;

    //now for each internal node find the observations that belong to the terminal nodes

    //NumericVector predictions(test_data.nrow());

    arma::mat pred_mat(testdata_arma.n_rows,num_cats);
    //arma::vec filled_in(testdata.nrow());


    //List term_obs(terminal_nodes.size());
    if(term_nodes.size()==1){

      //Rcout << "Line 422. \n";


      pred_mat=repmat(tree_table1(0,arma::span(5,5+num_cats-1)),testdata_arma.n_rows,1);


      //Rcout << "Line 424. \n";


      // for(int k=0; k<num_cats; k++){
      // pred_mat(_,k)=rep(treetable(0,5+k),testdata.nrow());
      // }
      //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
      //predictions=rep(nodemean,test_data.nrow());
      //Rcout << "Line 67 .\n";

      //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
      //term_obs[0]= temp_obsvec;
      // double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);
      //
      // double num_prod=1;
      // double num_sum=0;

      // for(int k=0; k<num_cats; k++){
      //   //assuming categories of y are from 1 to num_cats
      //   arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
      //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
      //   arma_tree_table(0,5+k)= m_plus_alph/denom_temp ;
      //
      //   //for likelihood calculation
      //   num_prod=num_prod*tgamma(m_plus_alph);
      //   num_sum=num_sum +m_plus_alph ;
      // }
      //
      // lik_prod= alph_term*num_prod/tgamma(num_sum);
      //
    }
    else{
      for(unsigned int i=0;i<term_nodes.size();i++){
        //arma::mat subdata=testd;
        int curr_term=term_nodes(i);

        int row_index;
        int term_node=term_nodes(i);


        //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
        //Why should the ro index be different for a right daughter?
        //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
        row_index=0;

        // if(curr_term % 2==0){
        //   //term node is left daughter
        //   row_index=terminal_nodes[i];
        // }else{
        //   //term node is right daughter
        //   row_index=terminal_nodes[i]-1;
        // }








        //save the left and right node data into arma uvec

        //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
        //arma::vec left_nodes=arma_tree.col(0);
        //arma::vec right_nodes=arma_tree.col(1);

        arma::vec left_nodes=tree_table1.col(0);
        arma::vec right_nodes=tree_table1.col(1);



        arma::mat node_split_mat;
        node_split_mat.set_size(0,3);
        //Rcout << "Line 124. i = " << i << " .\n";

        while(row_index!=1){
          //for each terminal node work backwards and see if the parent node was a left or right node
          //append split info to a matrix
          int rd=0;
          arma::uvec parent_node=arma::find(left_nodes == term_node);

          if(parent_node.size()==0){
            parent_node=arma::find(right_nodes == term_node);
            rd=1;
          }

          //want to cout parent node and append to node_split_mat

          node_split_mat.insert_rows(0,1);

          //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
          //node_split_mat(0,0)=treetable(parent_node[0],2);
          //node_split_mat(0,1)=treetable(parent_node[0],3);

          //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
          //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

          node_split_mat(0,0)=tree_table1(parent_node(0),2);
          node_split_mat(0,1)=tree_table1(parent_node(0),3);

          node_split_mat(0,2)=rd;
          row_index=parent_node(0)+1;
          term_node=parent_node(0)+1;
        }

        //once we have the split info, loop through rows and find the subset indexes for that terminal node!
        //then fill in the predicted value for that tree
        //double prediction = tree_data(term_node,5);
        arma::uvec pred_indices;
        int split= node_split_mat(0,0)-1;

        //arma::vec tempvec = testd.col(split);
        arma::vec tempvec = arma_test_data.col(split);


        double temp_split = node_split_mat(0,1);

        if(node_split_mat(0,2)==0){
          pred_indices = arma::find(tempvec <= temp_split);
        }else{
          pred_indices = arma::find(tempvec > temp_split);
        }

        arma::uvec temp_pred_indices;

        //arma::vec data_subset = testd.col(split);
        arma::vec data_subset = arma_test_data.col(split);

        data_subset=data_subset.elem(pred_indices);

        //now loop through each row of node_split_mat
        int n=node_split_mat.n_rows;
        //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
        //Rcout << "Line 174. node_split_mat= " << node_split_mat << ". n = " << n << ".\n";


        for(int j=1;j<n;j++){
          int curr_sv=node_split_mat(j,0);
          double split_p = node_split_mat(j,1);

          //data_subset = testd.col(curr_sv-1);
          data_subset = arma_test_data.col(curr_sv-1);

          data_subset=data_subset.elem(pred_indices);

          if(node_split_mat(j,2)==0){
            //split is to the left
            temp_pred_indices=arma::find(data_subset <= split_p);
          }else{
            //split is to the right
            temp_pred_indices=arma::find(data_subset > split_p);
          }
          pred_indices=pred_indices.elem(temp_pred_indices);

          if(pred_indices.size()==0){
            continue;
          }

        }
        //Rcout << "Line 199. i = " << i <<  ".\n";

        //double nodemean=tree_data(terminal_nodes[i]-1,5);
        //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
        //predictions[predind]= nodemean;
        //term_obs[i]=predind;

        //Rcout << "Line 635. \n";
        //Rcout << "pred_indices = " << pred_indices << ".\n";

        //pred_mat.rows(pred_indices)=arma::repmat(arma_tree_table(curr_term-1,arma::span(5,5+num_cats-1)),pred_indices.n_elem,1);
        pred_mat.each_row(pred_indices)=tree_table1(curr_term-1,arma::span(5,4+num_cats));



        //Rcout << "Line 588. \n";

        // for(int k=0; k<num_cats; k++){
        //   pred_mat(predind,k)=rep(treetable(curr_term-1,5+k),predind.size());
        // }


        // double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
        // //Rcout << "Line 207. predind = " << predind <<  ".\n";
        // //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
        // // << "Line 207. term_node = " << term_node <<  ".\n";
        //
        // double num_prod=1;
        // double num_sum=0;
        //
        // for(int k=0; k<num_cats; k++){
        //   //assuming categories of y are from 1 to num_cats
        //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
        //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        //
        //   arma_tree_table(curr_term-1,5+k)= m_plus_alph/denom_temp ;
        //
        //   num_prod=num_prod*tgamma(m_plus_alph);
        //   num_sum=num_sum +m_plus_alph ;
        // }


        //lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);


      }
    }






    //THIS SHOULD BE DIFFERENT IF THE CODE IS TO BE PARALLELIZED
    //EACH THREAD SHOULD OUTPUT ITS OWN MATRIX AND SUM OF LIKELIHOODS
    //THEN ADD THE MATRICES TOGETHER AND DIVIDE BY THE TOTAL SUM OF LIKELIHOODS
    //OR JUST SAVE ALL MATRICES TO ONE LIST

    pred_mat_overall = pred_mat_overall + templik*pred_mat;




    //arma::mat treeprob_output = get_test_probs(weights, num_cats,
    //                                           testdata,
    //                                           treetable_list[i]  );

    //Rcout << "Line 688. i== " << i << ". \n";

    //double weighttemp = weights[i];
    //Rcout << "Line 691. i== " << i << ". \n";

    //pred_mat_overall = pred_mat_overall + weighttemp*treeprob_output;



  }//end of loop over all trees




  ///////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////



  double sumlik_total= arma::sum(overall_liks);
  pred_mat_overall=pred_mat_overall*(1/sumlik_total);
  //Rcout << "Line 1141 . \n";
  //Rcout << "Line 1146 . \n";

  return(wrap(pred_mat_overall));

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]
#include <xoshiro.h>
#include <dqrng_distribution.h>
//#include <dqrng.h>

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

//' @title Parallel Safe-Bayesian Random Forest
//'
//' @description A parallelized implementation of the Safe-Bayesian Random Forest described by Quadrianto and Ghahramani (2015)
//' @param lambda A real number between 0 and 1 that determines the splitting probability in the prior (which is used as the importance sampler of tree models). Quadrianto and Ghahramani (2015) recommend a value less than 0.5 .
//' @param num_trees The number of trees to be sampled.
//' @param seed The seed for random number generation.
//' @param num_cats The number of possible values for the outcome variable.
//' @param y The training data vector of outcomes. This must be a vector of integers between 1 and num_cats.
//' @param original_datamat The original training data. Currently all variables must be continuous. The training data does not need to be transformed before being entered to this function.
//' @param alpha_parameters Vector of prior parameters.
//' @param beta_par The power to which the likelihood is to be raised. For BMA, set beta_par=1.
//' @param original_datamat The original test data. This matrix must have the same number of columns (variables) as the training data. Currently all variables must be continuous. The test data does not need to be transformed before being entered to this function.
//' @param ncores The number of cores to be used in parallelization.
//' @return A matrix of probabilities with the number of rows equl to the number of test observations and the number of columns equal to the number of possible outcome categories.
//' @export
// [[Rcpp::export]]
NumericVector sBART_onefunc_parallel(double lambda,
                                     int num_models,
                                     int num_trees,
                                        int seed,
                                        NumericVector ytrain,
                                        NumericMatrix original_datamat,
                                        double beta_par,
                                        NumericMatrix test_datamat,
                                        int ncores,
                                        int outsamppreds,
                                        double nu,
                                        double a,
                                        double lambdaBART,
                                        int valid_trees,
                                        int tree_prior,
                                        int imp_sampler,
                                        double alpha_BART,
                                        double beta_BART){


  //Rcout << "imp_sampler = " << imp_sampler << ".\n";

  NumericVector y_scaled=scale_response(min(ytrain),max(ytrain),-0.5,0.5,ytrain);

  int num_split_vars= original_datamat.ncol();
  arma::mat data_arma= as<arma::mat>(original_datamat);
  arma::mat testdata_arma= as<arma::mat>(test_datamat);
  arma::vec orig_y_arma= as<arma::vec>(y_scaled);
  //arma::vec alpha_pars_arma= as<arma::vec>(alpha_parameters);
  int num_obs = data_arma.n_rows;
  int num_test_obs = testdata_arma.n_rows;


  //calculations for likelihood
  arma::mat y(num_obs,1);
  y.col(0)=orig_y_arma;
  //get exponent
  double expon=(num_obs+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;


  ///////////////////////
  //NumericMatrix Data_transformed = cpptrans_cdf(original_datamat);
  // NumericMatrix Data_transformed(original_datamat.nrow(), original_datamat.ncol());
  // for(int i=0; i<original_datamat.ncol();i++){
  //   NumericVector samp= original_datamat(_,i);
  //   NumericVector sv(clone(samp));
  //   std::sort(sv.begin(), sv.end());
  //   double nobs = samp.size();
  //   NumericVector ans(nobs);
  //   for (int k = 0; k < samp.size(); ++k)
  //     ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
  //   //NumericVector ansnum = ans;
  //   Data_transformed(_,i) = (ans+1)/nobs;
  // }



  //arma::mat arma_orig_data(Data_transformed.begin(), Data_transformed.nrow(), Data_transformed.ncol(), false);



  //NumericMatrix transformedData(originaldata.nrow(), originaldata.ncol());

  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  arma::mat arma_orig_data(data_arma.n_rows,data_arma.n_cols);
  for(unsigned int k=0; k<data_arma.n_cols;k++){
    arma::vec samp= data_arma.col(k);
    arma::vec sv=arma::sort(samp);
    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      while (sv(j) < ssampi && j < sv.size()) ++j;
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }
    arma_orig_data.col(k)=(ans+1)/nobs;
  }





  /////////////////////////////////////
  // NumericMatrix testdat_trans = cpptrans_cdf_test(original_datamat,test_datamat);
  // //NumericMatrix testdat_trans(test_datamat.nrow(), test_datamat.ncol());
  // for(int i=0; i<test_datamat.ncol();i++){
  //   NumericVector samp= test_datamat(_,i);
  //   NumericVector svtest = original_datamat(_,i);
  //   NumericVector sv(clone(svtest));
  //   std::sort(sv.begin(), sv.end());
  //   double nobs = samp.size();
  //   NumericVector ans(nobs);
  //   double nobsref = svtest.size();
  //   for (int k = 0; k < samp.size(); ++k){
  //     ans[k] = std::lower_bound(sv.begin(), sv.end(), samp[k]) - sv.begin();
  //   }
  //   //NumericVector ansnum = ans;
  //   testdat_trans(_,i) = (ans)/nobsref;
  // }





  //NumericMatrix transformedData(originaldata.nrow(), originaldata.ncol());
  //arma::mat data_arma= as<arma::mat>(originaldata);

  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  arma::mat arma_test_data(testdata_arma.n_rows,testdata_arma.n_cols);
  for(unsigned int k=0; k<data_arma.n_cols;k++){
    arma::vec ref= data_arma.col(k);
    arma::vec samp= testdata_arma.col(k);

    arma::vec sv=arma::sort(samp);
    arma::vec sref=arma::sort(ref);

    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    double nobsref = ref.n_elem;

    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      if(j+1>sref.size()){
      }else{
        while (sref(j) < ssampi && j < sref.size()){
          ++j;
          if(j==sref.size()) break;
        }
      }
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }

    arma_test_data.col(k)=(ans)/nobsref;

  }







  /////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////////////
  //List table_list = draw_trees(lambda, num_trees, seed, num_split_vars, num_cats );



  //dqrng::dqRNGkind("Xoroshiro128+");
  //dqrng::dqset_seed(IntegerVector::create(seed));

  //use following with binomial?
  //dqrng::xoshiro256plus rng(seed);

  std::vector<double> lambdavec = {lambda, 1-lambda};

  //typedef boost::mt19937 RNGType;
  //boost::random::uniform_int_distribution<> sample_splitvardist(1,num_split_vars);
  //boost::variate_generator< RNGType, boost::uniform_int<> >  sample_splitvars(rng, sample_splitvardist);

  //boost::random::uniform_real_distribution<double> b_unifdist(0,1);
  //boost::variate_generator< RNGType, boost::uniform_real<> >  b_unif_point(rng, b_unifdist);



  std::random_device device;
  //std::mt19937 gen(device());

  //possibly use seed?
  //// std::mt19937 gen(seed);

  dqrng::xoshiro256plus gen(device());              // properly seeded rng

  //dqrng::xoshiro256plus gen(seed);              // properly seeded rng




  std::bernoulli_distribution coin_flip(lambda);

  std::uniform_int_distribution<> distsampvar(1, num_split_vars);
  std::uniform_real_distribution<> dis_cont_unif(0, 1);


  //dqrng::uniform_distribution dis_cont_unif(0.0, 1.0); // Uniform distribution [0,1)

  //Following three functions can't be used in parallel
  //dqrng::dqsample_int coin_flip2(2, 1, true,lambdavec );
  //dqrng::dqsample_int distsampvar(num_split_vars, 1, true);
  //dqrng::dqrunif dis_cont_unif(1, 0, 1);



  //arma::mat arma_test_data(testdat_trans.begin(), testdat_trans.nrow(), testdat_trans.ncol(), false);


  arma::vec pred_vec_overall=arma::zeros<arma::vec>(arma_test_data.n_rows);


  //arma::field<arma::mat> overall_treetables(num_models);

  arma::field<arma::vec> overall_preds(num_models);

  arma::vec overall_liks(num_models);


  //overall_treetables[i]= wrap(tree_table1);
  //double templik = as<double>(treepred_output[1]);
  //overall_liks[i]= pow(lik_prod,beta_pow);

  //Rcout << "Line 3338. \n";


#pragma omp parallel num_threads(ncores)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps

#pragma omp for
  for(int j=0; j<num_models;j++){

    arma::mat Wmat(num_obs,0);

    //maybe use line below, depends how Jmat joined to Wmat
    //int upsilon=0;

    arma::mat W_tilde(num_test_obs,0);

    //maybe use line below, depends how Jmat joined to Wmat
    //int upsilon2=0;

    //double sum_tree_samp_prob=1;
    //double sum_tree_prior_prob=1;

    double sum_prior_over_samp_prob=1;

    for(int q=0; q<num_trees;q++){  //start of loop over trees in sum


    //If parallelizing, define the distributinos before this loop
    //and use lrng and the following two lines
    //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
    //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


    //NumericVector treenodes_bin(0);
    //arma::uvec treenodes_bin(0);

    std::vector<int> treenodes_bin;


    int count_terminals = 0;
    int count_internals = 0;

    //int count_treebuild = 0;


    if(imp_sampler==1){ //If sampling from BART prior

      //std::bernoulli_distribution coin_flip2(lambda);
      double depth1=0;
      int prev_node=0; //1 if previous node splits, zero otherwise

      double samp_prob;

      while(count_internals > (count_terminals -1)){
        samp_prob=alpha_BART*pow(double(depth1+1),-beta_BART);
        std::bernoulli_distribution coin_flip2(samp_prob);

        int tempdraw = coin_flip2(lgen);
        treenodes_bin.push_back(tempdraw);

        if(tempdraw==1){

          depth1=depth1+1; //after a split, the depth will increase by 1
          prev_node=1;
          count_internals=count_internals+1;

        }else{

          if(prev_node==1){//zero following a 1, therefore at same depth.
            //Don't change depth. Do nothing
          }else{ //zero following a zero, therefore the depth will decrease by 1
            depth1=depth1-1;
          }
          prev_node=0;
          count_terminals=count_terminals+1;

        }

      }

    }else{  //If not sampling from BAT prior
      if(imp_sampler==2){//If sampling from spike and tree prior
        throw std::range_error("code not yet written for spike and tree sampling");

      }else{//If sampling from default Q+G prior. i.e. not sampling from BART nor spike and tree prior

        while(count_internals > (count_terminals -1)){

          //Also consider standard library and random header
          // std::random_device device;
          // std::mt19937 gen(device());
          // std::bernoulli_distribution coin_flip(lambda);
          // bool outcome = coin_flip(gen);


          int tempdraw = coin_flip(lgen);

          //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


          //int tempdraw = Rcpp::rbinom(1,lambda,1);
          //int tempdraw = R::rbinom(1,lambda);

          //Rcout << "tempdraw = " << tempdraw << ".\n" ;

          //int tempdraw = coin_flip2(lgen)-1;

          //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;


          //need to update rng if use boost?
          //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));

          treenodes_bin.push_back(tempdraw);


          if(tempdraw==1){
            count_internals=count_internals+1;
          }else{
            count_terminals=count_terminals+1;
          }

        }//end of while loop creating parent vector treenodes_bin
      }

    }



    //Consider making this an armadillo vector
    //IntegerVector split_var_vec(treenodes_bin.size());
    //arma::uvec split_var_vec(treenodes_bin.size());
    std::vector<int> split_var_vec(treenodes_bin.size());

    //loop drawing splitting variables
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_var_vec[i] = -1;
      }else{
        // also consider the standard library function uniform_int_distribution
        // might need random header
        // This uses the Mersenne twister

        //Three lines below should probably be outside all the loops
        // std::random_device rd;
        // std::mt19937 engine(rd());
        // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
        //
        // split_var_vec[i] = distsampvar(engine);

        split_var_vec[i] = distsampvar(lgen);


        //consider using boost
        //might need to update rng
        //split_var_vec[i] <- sample_splitvars(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

        //not sure if this returns an integer or a vector?
        //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
        //could try
        //split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
        //could also try RcppArmadillo::rmultinom

      }

    }// end of for-loop drawing split variables


    //Consider making this an armadillo vector
    //NumericVector split_point_vec(treenodes_bin.size());
    //arma::vec split_point_vec(treenodes_bin.size());
    std::vector<double> split_point_vec(treenodes_bin.size());


    //loop drawing splitting points
    //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

    //if using armadillo, it might be faster to subset to split nodes
    //then use a vector of draws
    for(unsigned int i=0; i<treenodes_bin.size();i++){
      if(treenodes_bin[i]==0){
        split_point_vec[i] = -1;
      }else{


        //////////////////////////////////////////////////////////
        //following function not reccommended
        //split_point_vec[i] = std::rand();
        //////////////////////////////////////////////////////////
        ////Standard library:
        ////This should probably be outside all the loops
        ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
        ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
        ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

        split_point_vec[i] = dis_cont_unif(lgen);

        //////////////////////////////////////////////////////////
        //from armadillo
        //split_point_vec[i] = arma::randu();

        //////////////////////////////////////////////////////////
        //probably not adviseable for paralelization
        //From Rcpp
        //split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

        //////////////////////////////////////////////////////////
        //consider using boost
        //might need to update rng
        //split_point_vec[i] <- b_unif_point(rng);

        //or use dqrng
        //not sure if have to update the random number
        //check if the following line is written properly
        //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

        //not sure if this returns an integer or a vector?





      }

    }// end of for-loop drawing split points





    //CODE FOR ADJUSTING SPLITTING POINTS SO THAT THE TREES ARE VALID
    if(valid_trees==1){
      for(unsigned int i=0; i<treenodes_bin.size();i++){ //loop over all nodes
        if(treenodes_bin[i]==1){ // if it is an internal node, then check for further splits on the same variable and update
          double first_split_var=split_var_vec[i];      //splitting variable to check for
          double first_split_point=split_point_vec[i];  //splitting point to use in updates

          double sub_int_nodes=0;       //this internal node count will be used to determine if in subtree relevant to sub_int_nodes
          double sub_term_nodes=0;      //this terminal node count will be used to determine if in subtree relevant to sub_int_nodes
          double preventing_updates=0; //indicates if still within subtree that is not to be updated
          double prevent_int_count=0;   //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
          double prevent_term_count=0;  //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
          for(unsigned int k=i+1; k<treenodes_bin.size();k++){
            if(treenodes_bin[k]==1){
              sub_int_nodes=sub_int_nodes+1;
              if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
                prevent_int_count=prevent_int_count+1;
              }
            }else{
              sub_term_nodes=sub_term_nodes+1;
              if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
                prevent_term_count=prevent_term_count+1;
              }
            }
            if(sub_int_nodes<=sub_term_nodes-2){
              break;
            }


            if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
              if(prevent_int_count>prevent_term_count-1){ //if this rule is satisfied then in subtree that is not to be updated
                continue; //still in subtree, therefore continue instead of checking for splits to be updates
              }else{
                preventing_updates=0; // no longer in subtree, therefore reset preventing_updates to zero
              }
            }


            if(sub_int_nodes>sub_term_nodes-1){
              if(treenodes_bin[k]==1){
                if(split_var_vec[k]==first_split_var){
                  split_point_vec[k]=split_point_vec[k]*first_split_point;
                  //beginning count of subtree that should not have
                  //further splits on first_split_var updated
                  preventing_updates=1; //indicates if still within subtree that is not to be updated
                  prevent_int_count=1;
                  prevent_term_count=0;
                }
              }
            }else{
              if(treenodes_bin[k]==1){
                if(split_var_vec[k]==first_split_var){
                  split_point_vec[k]=split_point_vec[k]+first_split_point-first_split_point*split_point_vec[k];
                  //beginning count of subtree that should not have
                  //further splits on first_split_var updated
                  preventing_updates=1; //indicates if still within subtree that is not to be updated
                  prevent_int_count=1;
                  prevent_term_count=0;
                }
              }
            }



          }//end of inner loop over k
        }//end of if statement treenodes_bin[i]==1)
      }//end of loop over i
    }//end of if statement valid_trees==1










    //Create tree table matrix

    //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

    //Rcout << "Line 1037. \n";
    //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

    //initialize with zeros. Not sure if this is necessary
    arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),6);
    //Rcout << "Line 1040. \n";


    //tree_table1(_,2) = wrap(split_var_vec);
    //tree_table1(_,3) = wrap(split_point_vec);
    //tree_table1(_,4) = wrap(treenodes_bin);

    //It might be more efficient to make everything an armadillo object initially
    // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
    arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
    arma::colvec split_point_vec_arma(split_point_vec);
    arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


    //Rcout << "Line 1054. \n";

    //Fill in splitting variable column
    tree_table1.col(2) = split_var_vec_arma;
    //Fill in splitting point column
    tree_table1.col(3) = split_point_vec_arma;
    //Fill in split/parent column
    tree_table1.col(4) = treenodes_bin_arma;


    //Rcout << "Line 1061. j = " << j << ". \n";



    // Now start filling in left daughter and right daughter columns
    std::vector<int> rd_spaces;
    int prev_node = -1;

    for(unsigned int i=0; i<treenodes_bin.size();i++){
      //Rcout << "Line 1061. i = " << i << ". \n";
      if(prev_node==0){
        //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
        //Rcout << "Line 1073. j = " << j << ". \n";

        tree_table1(rd_spaces.back(), 1)=i+1;
        //Rcout << "Line 1076. j = " << j << ". \n";

        rd_spaces.pop_back();
      }
      if(treenodes_bin[i]==1){
        //Rcout << "Line 1081. j = " << j << ". \n";

        tree_table1(i,0) = i+2;
        rd_spaces.push_back(i);
        prev_node = 1;
        //Rcout << "Line 185. j = " << j << ". \n";

      }else{                  // These 2 lines unnecessary if begin with matrix of zeros
        //Rcout << "Line 1089. j = " << j << ". \n";
        tree_table1(i,0)=0 ;
        tree_table1(i,1) = 0 ;
        prev_node = 0;
        //Rcout << "Line 1093. j = " << j << ". \n";

      }
    }//
    //Rcout << "Line 1097. j = " << j << ". \n";





    //List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
    //                                     originaldata,
    //                                     treetable_list[i]  );


    //use armadillo object tree_table1

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //create variables for likelihood calcuations
    // double lik_prod=1;
    // double alph_prod=1;
    // for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
    //   alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
    // }
    // double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
    // double alph_term=gam_alph_sum/alph_prod;

    //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
    //arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


    //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
    //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

    //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

    //NumericVector terminal_nodes=find_term_nodes(treetable);

    //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

    //arma::vec colmat=arma_tree.col(4);
    //arma::uvec term_nodes=arma::find(colmat==-1);

    //arma::vec colmat=arma_tree.col(2);
    //arma::uvec term_nodes=arma::find(colmat==0);

    //arma::vec colmat=tree_table1.col(4);
    //arma::uvec term_nodes=arma::find(colmat==0);

    //4th column is treenodes_bin_arma
    arma::uvec term_nodes=arma::find(treenodes_bin_arma==0);

    term_nodes=term_nodes+1;

    //NumericVector terminal_nodes= wrap(term_nodes);



    //GET J MATRIX

    arma::mat Jmat(num_obs,term_nodes.n_elem);
    arma::mat Jtilde(num_test_obs,term_nodes.n_elem);

    //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
    //NumericVector tree_predictions;

    //now for each internal node find the observations that belong to the terminal nodes

    //NumericVector predictions(test_data.nrow());
    //List term_obs(term_nodes.n_elem);

    //GET J MATRIX

    if(term_nodes.n_elem==1){
      //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
      //predictions=rep(nodemean,test_data.nrow());
      //Rcout << "Line 67 .\n";

      //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
      //term_obs[0]= temp_obsvec;
      //double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

      //double num_prod=1;
      //double num_sum=0;
      //Rcout << "Line 129.\n";
      Jmat.col(0) = arma::ones<arma::vec>(num_obs);
      Jtilde.col(0) = arma::ones<arma::vec>(num_test_obs);

      //for(int k=0; k<num_cats; k++){
        //assuming categories of y are from 1 to num_cats
        //arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
        //double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        //tree_table1(0,5+k)= m_plus_alph/denom_temp ;

        //for likelihood calculation
        //num_prod=num_prod*tgamma(m_plus_alph);
        //num_sum=num_sum +m_plus_alph ;
      //}

      //lik_prod= alph_term*num_prod/tgamma(num_sum);

    }
    else{
      for(unsigned int i=0;i<term_nodes.n_elem;i++){
        //arma::mat subdata=testd;
        //int curr_term=term_nodes(i);

        int row_index;
        int term_node=term_nodes(i);
        //Rcout << "Line 152.\n";


        //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
        //Why should the ro index be different for a right daughter?
        //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
        row_index=0;

        // if(curr_term % 2==0){
        //   //term node is left daughter
        //   row_index=terminal_nodes[i];
        // }else{
        //   //term node is right daughter
        //   row_index=terminal_nodes[i]-1;
        // }




        //save the left and right node data into arma uvec

        //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
        //arma::vec left_nodes=arma_tree.col(0);
        //arma::vec right_nodes=arma_tree.col(1);

        arma::vec left_nodes=tree_table1.col(0);
        arma::vec right_nodes=tree_table1.col(1);



        arma::mat node_split_mat;
        node_split_mat.set_size(0,3);
        //Rcout << "Line 182. i = " << i << " .\n";

        while(row_index!=1){
          //for each terminal node work backwards and see if the parent node was a left or right node
          //append split info to a matrix
          int rd=0;
          arma::uvec parent_node=arma::find(left_nodes == term_node);

          if(parent_node.size()==0){
            parent_node=arma::find(right_nodes == term_node);
            rd=1;
          }

          //want to cout parent node and append to node_split_mat

          node_split_mat.insert_rows(0,1);

          //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
          //node_split_mat(0,0)=treetable(parent_node[0],2);
          //node_split_mat(0,1)=treetable(parent_node[0],3);

          //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
          //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

          node_split_mat(0,0)=tree_table1(parent_node(0),2);
          node_split_mat(0,1)=tree_table1(parent_node(0),3);

          node_split_mat(0,2)=rd;
          row_index=parent_node(0)+1;
          term_node=parent_node(0)+1;
        }

        //once we have the split info, loop through rows and find the subset indexes for that terminal node!
        //then fill in the predicted value for that tree
        //double prediction = tree_data(term_node,5);
        arma::uvec pred_indices;
        arma::uvec pred_test_indices;
        int split= node_split_mat(0,0)-1;

        //Rcout << "Line 224.\n";
        //Rcout << "split = " << split << ".\n";
        //arma::vec tempvec = testd.col(split);
        arma::vec tempvec = arma_orig_data.col(split);
        arma::vec temptest_vec = arma_test_data.col(split);
        //Rcout << "Line 227.\n";


        double temp_split = node_split_mat(0,1);

        if(node_split_mat(0,2)==0){
          pred_indices = arma::find(tempvec <= temp_split);
          pred_test_indices = arma::find(temptest_vec <= temp_split);
        }else{
          pred_indices = arma::find(tempvec > temp_split);
          pred_test_indices = arma::find(temptest_vec > temp_split);
        }
        //Rcout << "Line 236.\n";

        arma::uvec temp_pred_indices;
        arma::uvec temp_test_pred_indices;

        //arma::vec data_subset = testd.col(split);
        arma::vec data_subset = arma_orig_data.col(split);
        arma::vec data_test_subset = arma_test_data.col(split);

        data_subset=data_subset.elem(pred_indices);
        data_test_subset=data_test_subset.elem(pred_test_indices);

        //now loop through each row of node_split_mat
        int n=node_split_mat.n_rows;
        //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
        //Rcout << "Line 248.\n";

        for(int j=1;j<n;j++){
          int curr_sv=node_split_mat(j,0);
          double split_p = node_split_mat(j,1);

          //data_subset = testd.col(curr_sv-1);
          //Rcout << "Line 255.\n";
          //Rcout << "curr_sv = " << curr_sv << ".\n";
          data_subset = arma_orig_data.col(curr_sv-1);
          data_test_subset = arma_test_data.col(curr_sv-1);
          //Rcout << "Line 258.\n";

          data_subset=data_subset.elem(pred_indices);
          data_test_subset=data_test_subset.elem(pred_test_indices);

          if(node_split_mat(j,2)==0){
            //split is to the left
            temp_pred_indices=arma::find(data_subset <= split_p);
            temp_test_pred_indices=arma::find(data_test_subset <= split_p);
          }else{
            //split is to the right
            temp_pred_indices=arma::find(data_subset > split_p);
            temp_test_pred_indices=arma::find(data_test_subset > split_p);
          }
          pred_indices=pred_indices.elem(temp_pred_indices);
          pred_test_indices=pred_test_indices.elem(temp_test_pred_indices);

          //if(pred_indices.size()==0){
          //  continue;
          //}

        }
        //Rcout << "Line 199. i = " << i <<  ".\n";

        //There is probably a more efficient way of doing this
        //e.g. initialize J matrix so that all elements are equal to zero
        arma::vec tempcol_J=arma::zeros<arma::vec>(num_obs);
        tempcol_J(pred_indices) = arma::ones<arma::vec>(pred_indices.size());
        Jmat.col(i) = tempcol_J;

        arma::vec tempcol_Jtilde=arma::zeros<arma::vec>(num_test_obs);
        tempcol_Jtilde(pred_test_indices) = arma::ones<arma::vec>(pred_test_indices.size());
        Jtilde.col(i) = tempcol_Jtilde;

        //double nodemean=tree_data(terminal_nodes[i]-1,5);
        //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
        //predictions[predind]= nodemean;
        //term_obs[i]=predind;

        //double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
        //Rcout << "Line 207. predind = " << predind <<  ".\n";
        //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
        // << "Line 207. term_node = " << term_node <<  ".\n";

        //double num_prod=1;
        //double num_sum=0;

        // for(int k=0; k<num_cats; k++){
        //   //assuming categories of y are from 1 to num_cats
        //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
        //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        //
        //   tree_table1(curr_term-1,5+k)= m_plus_alph/denom_temp ;
        //
        //   num_prod=num_prod*tgamma(m_plus_alph);
        //   num_sum=num_sum +m_plus_alph ;
        // }
        //
        //
        // lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);
        //Rcout << "Line 297.\n";


      }//End of loop over terminal nodes.
    }// end of else statement (for when more than one terminal node)
    // Now have J matrix

    Wmat=join_rows(Wmat,Jmat);
    //or
    //Wmat.insert_cols(Wmat.n_cols,Jmat);
    //or
    //int b_j=term_nodes.n_elem;
    //Wmat.insert_cols(upsilon,Jmat);
    //upsilon+=b_j;


    //Obtain test W_tilde, i.e. W matrix for test data

    W_tilde=join_rows(W_tilde,Jtilde);
    //or
    //W_tilde.insert_cols(W_tilde.n_cols,Jtilde);
    //or
    //int b_jtest=term_nodes.n_elem;
    //W_tilde.insert_cols(upsilon2,Jtilde);
    //upsilon2+=b_jtest;


    if(imp_sampler!=tree_prior){//check if importance sampler is not equal to the prior
      // //get impportance sampler probability and tree prior
      // long double temp_samp_prob;
      // long double temp_prior_prob;
      // //get sampler tree probability
      // if(imp_sampler==1){//If sample from BART prior
      //
      //
      //
      //   temp_samp_prob=1;
      //
      //   double depth1=0;
      //   int prev_node=0; //1 if previous node splits, zero otherwise
      //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
      //     if(treenodes_bin[i_2]==1){
      //       temp_samp_prob=temp_samp_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
      //       depth1=depth1+1; //after a split, the depth will increase by 1
      //       prev_node=1;
      //     }else{
      //       temp_samp_prob=temp_samp_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
      //       if(prev_node==1){//zero following a 1, therefore at same depth.
      //         //Don't change depth. Do nothing
      //       }else{ //zero following a zero, therefore the depth will decrease by 1
      //         depth1=depth1-1;
      //       }
      //       prev_node=0;
      //
      //     }
      //   }
      //
      //   //end of calculating BART tree probability
      // }else{
      //   if(imp_sampler==2){//If sample from spike and tree prior
      //     throw std::range_error("code not yet written for spike and tree prior");
      //
      //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
      //     double tempexp1=treenodes_bin.size()-arma::sum(treenodes_bin_arma);
      //     double tempexp2=arma::sum(treenodes_bin_arma);
      //     temp_samp_prob=pow(lambda,tempexp2)*
      //       pow(1-lambda,tempexp1);
      //       //(1/pow(double(num_split_vars),tempexp2));
      //
      //       temp_samp_prob=exp(log(lambda)*tempexp2+
      //         log(1-lambda)*tempexp1);
      //
      //     //temp_samp_prob=pow(lambda,arma::sum(treenodes_bin_arma))*
      //     //  pow(1-lambda,treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
      //     //  pow((1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
      //   }
      // }
      //
      // sum_tree_samp_prob=sum_tree_samp_prob*temp_samp_prob;
      // //end of getting importance sampler probability
      //
      // //get prior tree probability
      // if(tree_prior==1){//If sample from BART prior
      //
      //
      //
      //   temp_prior_prob=1;
      //
      //   double depth1=0;
      //   int prev_node=0; //1 if previous node splits, zero otherwise
      //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
      //
      //     if(treenodes_bin[i_2]==1){
      //       temp_prior_prob=temp_prior_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
      //       depth1=depth1+1; //after a split, the depth will increase by 1
      //       prev_node=1;
      //     }else{
      //       temp_prior_prob=temp_prior_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
      //       if(prev_node==1){//zero following a 1, therefore at same depth.
      //         //Don't change depth. Do nothing
      //       }else{ //zero following a zero, therefore the depth will decrease by 1
      //         depth1=depth1-1;
      //       }
      //       prev_node=0;
      //
      //     }
      //     //if(alpha_BART==0){
      //     //  Rcout << "alpha_BART equals zero!!!!.\n";
      //     //}
      //   }
      //
      //   //end of calculating BART tree probability
      // }else{
      //   if(tree_prior==2){//If sample from spike and tree prior
      //     throw std::range_error("code not yet written for spike and tree prior");
      //
      //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
      //     temp_prior_prob=pow((long double)(lambda),arma::sum(treenodes_bin_arma))*
      //       pow((long double)(1-lambda),treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
      //       pow((long double)(1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
      //   }
      // }
      //
      // sum_tree_prior_prob=sum_tree_prior_prob*temp_prior_prob;
      // if(temp_prior_prob==0){
      //   Rcout << "Line 4097, j= " << j << ". \n";
      //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
      //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
      // }
      // if(temp_samp_prob==0){
      //   Rcout << "Line 4102, j= " << j << ". \n";
      //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
      //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
      // }
      //
      // if(sum_tree_samp_prob==0){
      //   Rcout << "Line 4102, j= " << j << ". \n";
      //   Rcout << "sum_tree_samp_prob= " << sum_tree_samp_prob << ". \n";
      //   //Rcout << "treenodes_bin_arma= " << treenodes_bin_arma << ". \n";
      //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
      //
      // }




      //get tree prior over impportance sampler probability
      double tree_prior_over_samp_prob=1;
      if(imp_sampler==1){   //If sample from BART prior


        if(tree_prior==1){  //If tree prior is BART prior
          throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

        }else{
          if(tree_prior==2){  //If tree prior is spike-and-tree prior
            throw std::range_error("code not yet written for spike and tree prior");

          }else{//otherwise the tree prior is the Quadrianto and Ghahramani prior
            double depth1=0;
            int prev_node=0; //1 if previous node splits, zero otherwise
            for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
              if(treenodes_bin[i_2]==1){
                tree_prior_over_samp_prob=tree_prior_over_samp_prob*(lambda/(alpha_BART*pow(double(depth1+1),-beta_BART)));
                depth1=depth1+1; //after a split, the depth will increase by 1
                prev_node=1;
              }else{
                tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-lambda)/(1-alpha_BART*pow(double(depth1+1),-beta_BART)));
                if(prev_node==1){//zero following a 1, therefore at same depth.
                  //Don't change depth. Do nothing
                }else{ //zero following a zero, therefore the depth will decrease by 1
                  depth1=depth1-1;
                }
                prev_node=0;

              }
            }

          }
        }


      }else{// if not sampling from BART prior
        if(imp_sampler==2){//If sample from spike and tree prior
          throw std::range_error("code not yet written for sampling from spike and tree prior");

        }else{//otherwise sampling from Quadrianto and Ghahramani prior
          if(tree_prior==1){  //If tree prior is BART prior

            double depth1=0;
            int prev_node=0; //1 if previous node splits, zero otherwise
            for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
              if(treenodes_bin[i_2]==1){
                tree_prior_over_samp_prob=tree_prior_over_samp_prob*((alpha_BART*pow(double(depth1+1),-beta_BART))/lambda);
                depth1=depth1+1; //after a split, the depth will increase by 1
                prev_node=1;
              }else{
                tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-alpha_BART*pow(double(depth1+1),-beta_BART))/(1-lambda));
                if(prev_node==1){//zero following a 1, therefore at same depth.
                  //Don't change depth. Do nothing
                }else{ //zero following a zero, therefore the depth will decrease by 1
                  depth1=depth1-1;
                }
                prev_node=0;

              }//close (zero node) else stattement

            }//end for loop over i_2

          }else{
            if(tree_prior==2){  //If tree prior is spike-and-tree prior
              throw std::range_error("code not yet written for spike and tree prior");

            }else{
              throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

            }//close (not BART nor spike and tree prior) else statement
          }// close (not BART prior) else statememt




        }//close (not sampling from BART or spike and tree)  else statement

      }//close (not sampling from BART) else statement

      sum_prior_over_samp_prob=sum_prior_over_samp_prob*tree_prior_over_samp_prob;
      //end of getting tree prior over impportance sampler probability

      // if(sum_prior_over_samp_prob==0){
      //   Rcout << "Line 4266, j= " << j << ". \n";
      //   Rcout << "Line 4267, q= " << q << ". \n";
      //   Rcout << "sum_prior_over_samp_prob= " << sum_prior_over_samp_prob << ". \n";
      //
      // }else{
      //   Rcout << "Line 4266, j= " << j << ". \n";
      //   Rcout << "Line 4267, q= " << q << ". \n";
      //   Rcout << "sum_prior_over_samp_prob= " << sum_prior_over_samp_prob << ". \n";
      // }

    }//end of tree prior and importance sampler calculations



    } //end of loop over trees in sum


    //Obtain W matrix. If more than one tree in sum, need to join J matrices, possibly in loop over model trees above
    // i.e. add a loop from just within the start of the outer loop to here of length equal to the number of trees within the model
    // Create a Wmat with zero columns at start of loop, and join the Jmat at the end of each loop

    //for now, testing a one-tree model
    //replace Jmat with Wmat later


    //Obtain likelihood

    double b=Wmat.n_cols;

    //get t(y)inv(psi)J
    arma::mat ytW=y.t()*Wmat;
    //get t(J)inv(psi)J
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y;
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    //arma::mat rel=(b/2)*log(a)-(1/2)*log(det(sec_term))-expon*log(nu*lambdaBART - mvm +yty);


    //double val1;
    //double sign1;

    //log_det(val1, sign1, sec_term);
    //double templik0=exp(arma::as_scalar((b*0.5)*log(a)-0.5*val1-expon*log(nu*lambdaBART - mvm +yty)));

    double templik0=exp(arma::as_scalar((b*0.5)*log(a)-0.5*real(arma::log_det(sec_term))-expon*log(nu*lambdaBART - mvm +yty)));

    //double templik0=exp(arma::as_scalar((b*0.5)*log(a)-0.5*log(det(sec_term))-expon*log(nu*lambdaBART - mvm +yty)));



    // Rcout << "templik0= " << templik0 << ". \n";
    // Rcout << "b= " << b << ". \n";
    // Rcout << "(b*0.5)*log(a)= " << (b*0.5)*log(a) << ". \n";
    //
    // Rcout << "-0.5*log(det(sec_term))= " << -0.5*log(det(sec_term)) << ". \n";
    // Rcout << "det(sec_term)= " << det(sec_term) << ". \n";
    // Rcout << "arma::det(sec_term)= " << arma::det(sec_term) << ". \n";
    // Rcout << "arma::log_det(sec_term)= " << arma::log_det(sec_term) << ". \n";
    // Rcout << "real(arma::log_det(sec_term))= " << real(arma::log_det(sec_term)) << ". \n";
    // Rcout << "log(det(sec_term))= " << log(det(sec_term)) << ". \n";
    // Rcout << "log(arma::det(sec_term))= " << log(arma::det(sec_term)) << ". \n";
    //
    // Rcout << "-expon*log(nu*lambdaBART - mvm +yty)= " << -expon*log(nu*lambdaBART - mvm +yty) << ". \n";
    //
    // Rcout << "val= " << val << ". \n";


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //overall_treetables[j]= wrap(tree_table1);


    //double templik = as<double>(treepred_output[1]);

    double templik = pow(templik0,beta_par);
    if(imp_sampler!=tree_prior){//check if importance sampler is not equal to the prior
      //templik=templik*(sum_tree_prior_prob/sum_tree_samp_prob);
      templik=templik*sum_prior_over_samp_prob;

    }
    overall_liks(j)= templik;

    // if(std::isnan(templik)){
    // Rcout << "Line 3943, j= " << j << ". \n";
    // Rcout << "templik= " << templik << ". \n";
    // Rcout << "sum_tree_prior_prob= " << sum_tree_prior_prob << ". \n";
    // Rcout << "sum_tree_samp_prob= " << sum_tree_samp_prob << ". \n";
    // }


    //now fill in the predictions

    //If want tree tables with predictions filled in, use
    // arma::vec term_node_par_means = sec_term_inv*third_term;
    // //and would need to save a field of tree tables,
    // //add add a column, or begin with one more column
    // //then the first treetableF[0].n_rows elements of term_node_par_means
    // //give the first
    // int row_count1=0;
    // for(int tree_i=0; tree_i < treetableF.n_elem; tree_i++){
    //   tabletemp= treetableF(i);
    //   tabletemp.col(5) = term_node_par_means(arma::span(row_count1,tabletemp.n_rows));
    //   treetableF(i)=tabletemp;
    //   row_count1+=tabletemp.n_rows;
    // }
    //This would give an alternative method for obtaining test data predictions
    //Look up the terminal nodes and add the relevant terminal node parameters




    //arma::vec pred_vec(testdata_arma.n_rows);
    arma::vec preds_temp_arma= W_tilde*sec_term_inv*third_term;





    //THIS SHOULD BE DIFFERENT IF THE CODE IS TO BE PARALLELIZED
    //EACH THREAD SHOULD OUTPUT ITS OWN MATRIX AND SUM OF LIKELIHOODS
    //THEN ADD THE MATRICES TOGETHER AND DIVIDE BY THE TOTAL SUM OF LIKELIHOODS
    //OR JUST SAVE ALL MATRICES TO ONE LIST


    //pred_mat_overall = pred_mat_overall + templik*pred_mat;
    //overall_treetables(j)= pred_mat*templik;
    overall_preds(j)=preds_temp_arma*templik;

    //Rcout << "Line 3985, j= " << j << ". \n";
    //Rcout << "preds_temp_arma= " << preds_temp_arma << ". \n";

    //overall_treetables(j)= pred_mat;
    //overall_liks(j) =templik;

    //arma::mat treeprob_output = get_test_probs(weights, num_cats,
    //                                           testdata,
    //                                           treetable_list[i]  );

    //Rcout << "Line 688. i== " << i << ". \n";

    //double weighttemp = weights[i];
    //Rcout << "Line 691. i== " << i << ". \n";

    //pred_mat_overall = pred_mat_overall + weighttemp*treeprob_output;



  }//end of loop over all trees

}//end of pragma omp code


///////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////


//for(unsigned int i=0; i<overall_treetables.n_elem;i++){
//  pred_mat_overall = pred_mat_overall + overall_liks(i)*overall_treetables(i);
//}


#pragma omp parallel
{
  arma::vec result_private=arma::zeros<arma::vec>(arma_test_data.n_rows);
#pragma omp for nowait //fill result_private in parallel
  for(unsigned int i=0; i<overall_preds.size(); i++) result_private += overall_preds(i);
#pragma omp critical
  pred_vec_overall += result_private;
}



//Rcout << "Line 4030. \n";






double sumlik_total= arma::sum(overall_liks);
//Rcout << "sumlik_total = " << sumlik_total << ". \n";

pred_vec_overall=pred_vec_overall*(1/sumlik_total);
//Rcout << "Line 1141 . \n";
//Rcout << "Line 1146 . \n";

//Rcout << "Line 4042. \n";
NumericVector orig_preds=get_original(min(ytrain),max(ytrain),-0.5,0.5,wrap(pred_vec_overall)) ;


return(orig_preds);

}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]
#include <xoshiro.h>
#include <dqrng_distribution.h>
//#include <dqrng.h>

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

//' @title Parallel Safe-Bayesian Causal Forest
//'
//' @description A parallelized implementation of the Safe-Bayesian Random Forest described by Quadrianto and Ghahramani (2015)
//' @param lambda A real number between 0 and 1 that determines the splitting probability in the prior (which is used as the importance sampler of tree models). Quadrianto and Ghahramani (2015) recommend a value less than 0.5 .
//' @param num_trees The number of trees to be sampled.
//' @param seed The seed for random number generation.
//' @param num_cats The number of possible values for the outcome variable.
//' @param y The training data vector of outcomes. This must be a vector of integers between 1 and num_cats.
//' @param original_datamat The original training data. Currently all variables must be continuous. The training data does not need to be transformed before being entered to this function.
//' @param alpha_parameters Vector of prior parameters.
//' @param beta_par The power to which the likelihood is to be raised. For BMA, set beta_par=1.
//' @param original_datamat The original test data. This matrix must have the same number of columns (variables) as the training data. Currently all variables must be continuous. The test data does not need to be transformed before being entered to this function.
//' @param ncores The number of cores to be used in parallelization.
//' @return A matrix of probabilities with the number of rows equl to the number of test observations and the number of columns equal to the number of possible outcome categories.
//' @export
// [[Rcpp::export]]
NumericVector sBCF_onefunc_parallel(double lambda_mu,
                                    double lambda_tau,
                                     int num_models,
                                     int num_trees_mu,
                                     int num_trees_tau,
                                     int seed,
                                     NumericVector ytrain,
                                     NumericMatrix original_datamat,
                                     NumericVector ztrain,
                                     NumericMatrix pihat_train,
                                     double beta_par,
                                     NumericMatrix test_datamat,
                                     NumericVector test_z,
                                     NumericMatrix test_pihat,
                                     int ncores,
                                     int outsamppreds,
                                     double nu,
                                     double a_mu,
                                     double a_tau,
                                     double lambdaBCF,
                                     int valid_trees,
                                     int tree_prior,
                                     int imp_sampler,
                                     double alpha_BCF_mu,
                                     double beta_BCF_mu,
                                     double alpha_BCF_tau,
                                     double beta_BCF_tau,
                                     int include_pi2){


  //Check that various input vectors and matrices have consistent dimensions

  //Rcout << "Line 4528.\n";

  bool is_test_data=0;					// create bool is_test_data. Initialize equal to 0.
  if(test_datamat.nrow()>0){					// If test data has non-zero number of rows.
    is_test_data=1;						// set is_test_data equal to 1.
  }
  if(ytrain.size() !=original_datamat.nrow()){				// If the length of input vector y is not equal to the nunber of rows in the input data (covariates)
    if(ytrain.size()<original_datamat.nrow()){			// If the length of y is less than the number of rows in data
      throw std::range_error("Response length is smaller than the number of observations in the data");
    }else{								// If the length of y is greater than the number of rows in data
      throw std::range_error("Response length is greater than the number of observations in the data");
    }
  }
  if(ztrain.size() !=original_datamat.nrow()){				// If the length of input vector z is not equal to the nunber of rows in the input data (covariates)
    if(ztrain.size()<original_datamat.nrow()){			// If the length of z is less than the number of rows in data
      throw std::range_error("Treatment indicator vector length is smaller than the number of observations in the data");
    }else{								// If the length of z is greater than the number of rows in data
      throw std::range_error("Treatment indicator vector length is greater than the number of observations in the data");
    }
  }
  if(pihat_train.nrow() !=original_datamat.nrow()){				// If the nunber of rows in the input matrix pihat is not equal to the nunber of rows in the input data (covariates)
    if(pihat_train.nrow()<original_datamat.nrow()){			// If the nunber of rows in the input matrix pihat is less than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat_train is smaller than the number of observations in the data");
    }else{								// If the nunber of rows in the input matrix pihat is greater than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat_train is greater than the number of observations in the data");
    }
  }
  //check test data has the same number of variables as training data
  if(test_datamat.nrow()>0 && (original_datamat.ncol() != test_datamat.ncol())){	// If the number of rows in the test data is >0 AND the number of columns (variables) is not equal to that of data (the training data)
    throw std::range_error("Test data and training data must have the same number of variables. BART BMA assumes variables are in the same order.");
  }
  if(test_z.size() != test_datamat.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data treatment indicator variable
    throw std::range_error("Test data covariates and test data treatment indicator variable must have the same number of observations.");
  }
  if(test_datamat.nrow() != test_pihat.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data propensity score estimates matrix
    throw std::range_error("Test data covariates and test data propensity score estimates must have the same number of observations.");
  }
  if(test_pihat.nrow()>0 && (pihat_train.ncol() != test_pihat.ncol())){	// If the number of rows in the test data propensity score estimates is >0 AND the number of columns (variables) is not equal to that of the training data propensity score estimates
    throw std::range_error("Test data propensity score estimates and training data propensity score estimates must have the same number of columns. BART BMA assumes variables are in the same order.");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////



  // Now add propensity score estimates matrix as new leftmost column of data matrix. Call the resulting matrix x_control (to be consistent with terminology used by bcf package).
  arma::mat D1(original_datamat.begin(), original_datamat.nrow(), original_datamat.ncol(), false);				// copy the covariate data matrix into an arma mat
  arma::mat pihat_1(pihat_train.begin(), pihat_train.nrow(), pihat_train.ncol(), false);				// copy the pihat matrix into an arma mat
  //arma::mat x_control_a=D1;				// create a copy of data arma mat called x_control_a


  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  arma::mat x_control_a_temp(D1.n_rows,D1.n_cols);
  for(unsigned int k=0; k<D1.n_cols;k++){
    arma::vec samp= D1.col(k);
    arma::vec sv=arma::sort(samp);
    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      while (sv(j) < ssampi && j < sv.size()) ++j;
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }
    x_control_a_temp.col(k)=(ans+1)/nobs;
  }

  arma::mat x_control_a=x_control_a_temp;			// create arma mat copy of x_control_a_temp.

  arma::mat x_moderate_a=x_control_a_temp;			// create arma mat copy of x_control_a_temp.

  arma::mat pihat_a(pihat_1.n_rows,pihat_1.n_cols);

  //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
  for(unsigned int k=0; k<pihat_1.n_cols;k++){
    arma::vec samp= pihat_1.col(k);
    arma::vec sv=arma::sort(samp);
    //std::sort(sv.begin(), sv.end());
    arma::uvec ord = arma::sort_index(samp);
    double nobs = samp.n_elem;
    arma::vec ans(nobs);
    for (unsigned int i = 0, j = 0; i < nobs; ++i) {
      int ind=ord(i);
      double ssampi(samp[ind]);
      while (sv(j) < ssampi && j < sv.size()) ++j;
      ans(ind) = j;     // j is the 1-based index of the lower bound
    }
    pihat_a.col(k)=(ans+1)/nobs;
  }



  if((include_pi2==0) | (include_pi2==2) ){
    if(pihat_train.nrow()>0 ){
      x_control_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
  // Rcout << "Number of columns of matrix" << x_control_a.n_cols << ".\n";


  //NumericMatrix x_control=wrap(x_control_a);	// convert x_control_a to a NumericMatrix called x_control

  // Name the matrix without the estimated propensity scores x_moderate.[CAN REMOVE THE DUPLICATION AND ADD x_control, x_moderate, and include_pi as input parameters later]
  //NumericMatrix x_moderate = data;	// x_moderate matrix is the covariate data without the propensity scores
  //arma::mat x_moderate_a=D1;			// create arma mat copy of x_moderate.
  if((include_pi2==1)| (include_pi2==2) ){
    if(pihat_train.nrow()>0 ){
      x_moderate_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
    }
  }


  //NumericMatrix x_moderate=wrap(x_moderate_a);	// convert x_control_a to a NumericMatrix called x_control


  // Rcout << "Get to Line 7139  "  << ".\n";
  // Add test propensity scores to test data matrix
  arma::mat T1(test_datamat.begin(), test_datamat.nrow(), test_datamat.ncol(), false);				// copy the covariate test_data matrix into an arma mat
  arma::mat pihat_1_test(test_pihat.begin(), test_pihat.nrow(), test_pihat.ncol(), false);				// copy the test_pihat matrix into an arma mat
  //arma::mat x_control_test_a=T1;				// create a copy of test_data arma mat called x_control_test_a

  arma::mat x_control_test_a(T1.n_rows,T1.n_cols);
  arma::mat x_moderate_test_a(T1.n_rows,T1.n_cols);
  arma::mat pihat_a_test(pihat_1_test.n_rows,pihat_1_test.n_cols);

  if(is_test_data==1){
    //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
    arma::mat x_control_a_test_temp(T1.n_rows,T1.n_cols);

    for(unsigned int k=0; k<T1.n_cols;k++){
      arma::vec samp= T1.col(k);
      arma::vec sv=arma::sort(samp);
      //std::sort(sv.begin(), sv.end());
      arma::uvec ord = arma::sort_index(samp);
      double nobs = samp.n_elem;
      arma::vec ans(nobs);
      for (unsigned int i = 0, j = 0; i < nobs; ++i) {
        int ind=ord(i);
        double ssampi(samp[ind]);
        while (sv(j) < ssampi && j < sv.size()) ++j;
        ans(ind) = j;     // j is the 1-based index of the lower bound
      }
      x_control_a_test_temp.col(k)=(ans+1)/nobs;
    }

    arma::mat x_control_test_a=x_control_a_test_temp;			// create arma mat copy of x_control_a_temp.

    arma::mat x_moderate_test_a=x_control_a_test_temp;			// create arma mat copy of x_control_a_temp.


    //THIS CAN BE PARALLELIZED IF THERE ARE MANY VARIABLES
    for(unsigned int k=0; k<pihat_1_test.n_cols;k++){
      arma::vec samp= pihat_1_test.col(k);
      arma::vec sv=arma::sort(samp);
      //std::sort(sv.begin(), sv.end());
      arma::uvec ord = arma::sort_index(samp);
      double nobs = samp.n_elem;
      arma::vec ans(nobs);
      for (unsigned int i = 0, j = 0; i < nobs; ++i) {
        int ind=ord(i);
        double ssampi(samp[ind]);
        while (sv(j) < ssampi && j < sv.size()) ++j;
        ans(ind) = j;     // j is the 1-based index of the lower bound
      }
      pihat_a_test.col(k)=(ans+1)/nobs;
    }

  }



  if((include_pi2==0)| (include_pi2==2) ){
    if(test_pihat.nrow()>0 ){
      x_control_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_test_a
    }
  }


  //NumericMatrix x_control_test=wrap(x_control_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test


  // Name the matrix without the estimated propensity scores x_moderate_test.[CAN REMOVE THE DUPLICATION AND ADD x_control_test, x_moderate_test, and include_pi as input parameters later]
  //NumericMatrix x_moderate_test = test_data;	// x_moderate_test matrix is the covariate test_data without the propensity scores
  //arma::mat x_moderate_test_a=T1;			// create arma mat copy of x_moderate_test.
  if((include_pi2==1)| (include_pi2==2) ){
    if(test_pihat.nrow()>0 ){
      x_moderate_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_a
    }
  }

  //NumericMatrix x_moderate_test=wrap(x_moderate_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test



  //////////////////////////////////////////////////////////////////////////////////////////
  //Rcout << "Line 4715.\n";

  //////////////////////////////////////////////////////////////////////////////////////////

  //End of checks and adding propensity scores to matrices

  NumericVector y_scaled=scale_response(min(ytrain),max(ytrain),-0.5,0.5,ytrain);

  arma::vec z_ar=Rcpp::as<arma::vec>(ztrain);		// converts to arma vec


  int num_split_vars_mu= x_control_a.n_cols;

  int num_split_vars_tau= x_moderate_a.n_cols;


  //Rcout << "num_split_vars_mu = " << num_split_vars_mu << ".\n" ;
  //Rcout << "num_split_vars_tau = " << num_split_vars_tau << ".\n" ;

  //arma::mat data_arma= as<arma::mat>(original_datamat);
  //arma::mat testdata_arma= as<arma::mat>(test_datamat);


  arma::vec orig_y_arma= as<arma::vec>(y_scaled);
  //arma::vec alpha_pars_arma= as<arma::vec>(alpha_parameters);
  int num_obs = x_control_a.n_rows;
  int num_test_obs = x_control_test_a.n_rows;


  //calculations for likelihood
  arma::mat y(num_obs,1);
  y.col(0)=orig_y_arma;
  //get exponent
  double expon=(num_obs+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;







  /////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////////////
  //List table_list = draw_trees(lambda, num_trees, seed, num_split_vars, num_cats );



  //dqrng::dqRNGkind("Xoroshiro128+");
  //dqrng::dqset_seed(IntegerVector::create(seed));

  //use following with binomial?
  //dqrng::xoshiro256plus rng(seed);

  std::vector<double> lambdavec_mu = {lambda_mu, 1-lambda_mu};
  std::vector<double> lambdavec_tau = {lambda_tau, 1-lambda_tau};

  //typedef boost::mt19937 RNGType;
  //boost::random::uniform_int_distribution<> sample_splitvardist(1,num_split_vars);
  //boost::variate_generator< RNGType, boost::uniform_int<> >  sample_splitvars(rng, sample_splitvardist);

  //boost::random::uniform_real_distribution<double> b_unifdist(0,1);
  //boost::variate_generator< RNGType, boost::uniform_real<> >  b_unif_point(rng, b_unifdist);



  std::random_device device;
  //std::mt19937 gen(device());

  //possibly use seed?
  //// std::mt19937 gen(seed);

  dqrng::xoshiro256plus gen(device());              // properly seeded rng

  //dqrng::xoshiro256plus gen(seed);              // properly seeded rng




  std::bernoulli_distribution coin_flip_mu(lambda_mu);
  std::bernoulli_distribution coin_flip_tau(lambda_tau);

  std::uniform_int_distribution<> distsampvar_mu(1, num_split_vars_mu);
  std::uniform_int_distribution<> distsampvar_tau(1, num_split_vars_tau);

  std::uniform_real_distribution<> dis_cont_unif(0, 1);


  //dqrng::uniform_distribution dis_cont_unif(0.0, 1.0); // Uniform distribution [0,1)

  //Following three functions can't be used in parallel
  //dqrng::dqsample_int coin_flip2(2, 1, true,lambdavec );
  //dqrng::dqsample_int distsampvar(num_split_vars, 1, true);
  //dqrng::dqrunif dis_cont_unif(1, 0, 1);



  //arma::mat arma_test_data(testdat_trans.begin(), testdat_trans.nrow(), testdat_trans.ncol(), false);


  arma::vec pred_vec_overall;
  if(is_test_data==1){
    pred_vec_overall=arma::zeros<arma::vec>(x_moderate_test_a.n_rows);
  }else{
    pred_vec_overall=arma::zeros<arma::vec>(x_moderate_a.n_rows);
  }

  //arma::field<arma::mat> overall_treetables(num_models);

  arma::field<arma::vec> overall_preds(num_models);

  arma::vec overall_liks(num_models);


  //overall_treetables[i]= wrap(tree_table1);
  //double templik = as<double>(treepred_output[1]);
  //overall_liks[i]= pow(lik_prod,beta_pow);

  //Rcout << "Line 3338. \n";

  //Rcout << "Line 4836.\n";

#pragma omp parallel num_threads(ncores)
{//start of pragma omp code
  dqrng::xoshiro256plus lgen(gen);      // make thread local copy of rng
  lgen.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... ncores jumps

#pragma omp for
  for(int j=0; j<num_models;j++){

    arma::mat Wmat_mu(num_obs,0);
    arma::mat Wmat_tau(num_obs,0);

    //maybe use line below, depends how Jmat joined to Wmat
    //int upsilon=0;

    arma::mat W_tilde_mu(num_test_obs,0);
    arma::mat W_tilde_tau(num_test_obs,0);

    //maybe use line below, depends how Jmat joined to Wmat
    //int upsilon2=0;

    //double sum_tree_samp_prob=1;
    //double sum_tree_prior_prob=1;

    double sum_prior_over_samp_prob=1;

    for(int q=0; q<num_trees_mu;q++){  //start of loop over trees in sum


      //If parallelizing, define the distributinos before this loop
      //and use lrng and the following two lines
      //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
      //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


      //NumericVector treenodes_bin(0);
      //arma::uvec treenodes_bin(0);

      std::vector<int> treenodes_bin;


      int count_terminals = 0;
      int count_internals = 0;

      //int count_treebuild = 0;


      if(imp_sampler==1){ //If sampling from BART prior

        double depth1=0;
        int prev_node=0; //1 if previous node splits, zero otherwise

        double samp_prob;

        while(count_internals > (count_terminals -1)){
          samp_prob=alpha_BCF_mu*pow(double(depth1+1),-beta_BCF_mu);
          std::bernoulli_distribution coin_flip2(samp_prob);

          int tempdraw = coin_flip2(lgen);
          treenodes_bin.push_back(tempdraw);

          if(tempdraw==1){

            depth1=depth1+1; //after a split, the depth will increase by 1
            prev_node=1;
            count_internals=count_internals+1;

          }else{

            if(prev_node==1){//zero following a 1, therefore at same depth.
              //Don't change depth. Do nothing
            }else{ //zero following a zero, therefore the depth will decrease by 1
              depth1=depth1-1;
            }
            prev_node=0;
            count_terminals=count_terminals+1;

          }

        }

      }else{  //If not sampling from BAT prior
        if(imp_sampler==2){//If sampling from spike and tree prior
          throw std::range_error("code not yet written for spike and tree sampling");

        }else{//If sampling from default Q+G prior. i.e. not sampling from BART nor spike and tree prior

          while(count_internals > (count_terminals -1)){

            //Also consider standard library and random header
            // std::random_device device;
            // std::mt19937 gen(device());
            // std::bernoulli_distribution coin_flip(lambda);
            // bool outcome = coin_flip(gen);


            int tempdraw = coin_flip_mu(lgen);

            //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


            //int tempdraw = Rcpp::rbinom(1,lambda,1);
            //int tempdraw = R::rbinom(1,lambda);

            ////Rcout << "tempdraw = " << tempdraw << ".\n" ;

            //int tempdraw = coin_flip2(lgen)-1;

            //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;


            //need to update rng if use boost?
            //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));

            treenodes_bin.push_back(tempdraw);


            if(tempdraw==1){
              count_internals=count_internals+1;
            }else{
              count_terminals=count_terminals+1;
            }

          }//end of while loop creating parent vector treenodes_bin
        }

      }



      //Consider making this an armadillo vector
      //IntegerVector split_var_vec(treenodes_bin.size());
      //arma::uvec split_var_vec(treenodes_bin.size());
      std::vector<int> split_var_vec(treenodes_bin.size());

      //loop drawing splitting variables
      //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

      //if using armadillo, it might be faster to subset to split nodes
      //then use a vector of draws
      for(unsigned int i=0; i<treenodes_bin.size();i++){
        if(treenodes_bin[i]==0){
          split_var_vec[i] = -1;
        }else{
          // also consider the standard library function uniform_int_distribution
          // might need random header
          // This uses the Mersenne twister

          //Three lines below should probably be outside all the loops
          // std::random_device rd;
          // std::mt19937 engine(rd());
          // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
          //
          // split_var_vec[i] = distsampvar(engine);

          split_var_vec[i] = distsampvar_mu(lgen);


          //consider using boost
          //might need to update rng
          //split_var_vec[i] <- sample_splitvars(rng);

          //or use dqrng
          //not sure if have to update the random number
          //check if the following line is written properly
          //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

          //not sure if this returns an integer or a vector?
          //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
          //could try
          //split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
          //could also try RcppArmadillo::rmultinom

        }

      }// end of for-loop drawing split variables


      //Consider making this an armadillo vector
      //NumericVector split_point_vec(treenodes_bin.size());
      //arma::vec split_point_vec(treenodes_bin.size());
      std::vector<double> split_point_vec(treenodes_bin.size());


      //loop drawing splitting points
      //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

      //if using armadillo, it might be faster to subset to split nodes
      //then use a vector of draws
      for(unsigned int i=0; i<treenodes_bin.size();i++){
        if(treenodes_bin[i]==0){
          split_point_vec[i] = -1;
        }else{


          //////////////////////////////////////////////////////////
          //following function not reccommended
          //split_point_vec[i] = std::rand();
          //////////////////////////////////////////////////////////
          ////Standard library:
          ////This should probably be outside all the loops
          ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
          ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
          ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

          split_point_vec[i] = dis_cont_unif(lgen);

          //////////////////////////////////////////////////////////
          //from armadillo
          //split_point_vec[i] = arma::randu();

          //////////////////////////////////////////////////////////
          //probably not adviseable for paralelization
          //From Rcpp
          //split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

          //////////////////////////////////////////////////////////
          //consider using boost
          //might need to update rng
          //split_point_vec[i] <- b_unif_point(rng);

          //or use dqrng
          //not sure if have to update the random number
          //check if the following line is written properly
          //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

          //not sure if this returns an integer or a vector?





        }

      }// end of for-loop drawing split points





      //CODE FOR ADJUSTING SPLITTING POINTS SO THAT THE TREES ARE VALID
      if(valid_trees==1){
        for(unsigned int i=0; i<treenodes_bin.size();i++){ //loop over all nodes
          if(treenodes_bin[i]==1){ // if it is an internal node, then check for further splits on the same variable and update
            double first_split_var=split_var_vec[i];      //splitting variable to check for
            double first_split_point=split_point_vec[i];  //splitting point to use in updates

            double sub_int_nodes=0;       //this internal node count will be used to determine if in subtree relevant to sub_int_nodes
            double sub_term_nodes=0;      //this terminal node count will be used to determine if in subtree relevant to sub_int_nodes
            double preventing_updates=0; //indicates if still within subtree that is not to be updated
            double prevent_int_count=0;   //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
            double prevent_term_count=0;  //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
            for(unsigned int k=i+1; k<treenodes_bin.size();k++){
              if(treenodes_bin[k]==1){
                sub_int_nodes=sub_int_nodes+1;
                if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
                  prevent_int_count=prevent_int_count+1;
                }
              }else{
                sub_term_nodes=sub_term_nodes+1;
                if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
                  prevent_term_count=prevent_term_count+1;
                }
              }
              if(sub_int_nodes<=sub_term_nodes-2){
                break;
              }


              if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
                if(prevent_int_count>prevent_term_count-1){ //if this rule is satisfied then in subtree that is not to be updated
                  continue; //still in subtree, therefore continue instead of checking for splits to be updates
                }else{
                  preventing_updates=0; // no longer in subtree, therefore reset preventing_updates to zero
                }
              }


              if(sub_int_nodes>sub_term_nodes-1){
                if(treenodes_bin[k]==1){
                  if(split_var_vec[k]==first_split_var){
                    split_point_vec[k]=split_point_vec[k]*first_split_point;
                    //beginning count of subtree that should not have
                    //further splits on first_split_var updated
                    preventing_updates=1; //indicates if still within subtree that is not to be updated
                    prevent_int_count=1;
                    prevent_term_count=0;
                  }
                }
              }else{
                if(treenodes_bin[k]==1){
                  if(split_var_vec[k]==first_split_var){
                    split_point_vec[k]=split_point_vec[k]+first_split_point-first_split_point*split_point_vec[k];
                    //beginning count of subtree that should not have
                    //further splits on first_split_var updated
                    preventing_updates=1; //indicates if still within subtree that is not to be updated
                    prevent_int_count=1;
                    prevent_term_count=0;
                  }
                }
              }



            }//end of inner loop over k
          }//end of if statement treenodes_bin[i]==1)
        }//end of loop over i
      }//end of if statement valid_trees==1





      //Rcout << "Line 5150.\n";





      //Create tree table matrix

      //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

      //Rcout << "Line 1037. \n";
      //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

      //initialize with zeros. Not sure if this is necessary
      arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),6);
      //Rcout << "Line 1040. \n";


      //tree_table1(_,2) = wrap(split_var_vec);
      //tree_table1(_,3) = wrap(split_point_vec);
      //tree_table1(_,4) = wrap(treenodes_bin);

      //It might be more efficient to make everything an armadillo object initially
      // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
      arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
      arma::colvec split_point_vec_arma(split_point_vec);
      arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


      //Rcout << "Line 1054. \n";

      //Fill in splitting variable column
      tree_table1.col(2) = split_var_vec_arma;
      //Fill in splitting point column
      tree_table1.col(3) = split_point_vec_arma;
      //Fill in split/parent column
      tree_table1.col(4) = treenodes_bin_arma;


      //Rcout << "Line 5189. j = " << j << ". \n";
      //Rcout << "Line 5190. tree_table1 mu = " << tree_table1 << ". \n";



      // Now start filling in left daughter and right daughter columns
      std::vector<int> rd_spaces;
      int prev_node = -1;

      for(unsigned int i=0; i<treenodes_bin.size();i++){
        //Rcout << "Line 1061. i = " << i << ". \n";
        if(prev_node==0){
          //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
          //Rcout << "Line 1073. j = " << j << ". \n";

          tree_table1(rd_spaces.back(), 1)=i+1;
          //Rcout << "Line 1076. j = " << j << ". \n";

          rd_spaces.pop_back();
        }
        if(treenodes_bin[i]==1){
          //Rcout << "Line 1081. j = " << j << ". \n";

          tree_table1(i,0) = i+2;
          rd_spaces.push_back(i);
          prev_node = 1;
          //Rcout << "Line 185. j = " << j << ". \n";

        }else{                  // These 2 lines unnecessary if begin with matrix of zeros
          //Rcout << "Line 1089. j = " << j << ". \n";
          tree_table1(i,0)=0 ;
          tree_table1(i,1) = 0 ;
          prev_node = 0;
          //Rcout << "Line 1093. j = " << j << ". \n";

        }
      }//
      //Rcout << "Line 1097. j = " << j << ". \n";





      //List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
      //                                     originaldata,
      //                                     treetable_list[i]  );


      //use armadillo object tree_table1

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////


      //create variables for likelihood calcuations
      // double lik_prod=1;
      // double alph_prod=1;
      // for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
      //   alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
      // }
      // double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
      // double alph_term=gam_alph_sum/alph_prod;

      //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
      //arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


      //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
      //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

      //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

      //NumericVector terminal_nodes=find_term_nodes(treetable);

      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

      //arma::vec colmat=arma_tree.col(4);
      //arma::uvec term_nodes=arma::find(colmat==-1);

      //arma::vec colmat=arma_tree.col(2);
      //arma::uvec term_nodes=arma::find(colmat==0);

      //arma::vec colmat=tree_table1.col(4);
      //arma::uvec term_nodes=arma::find(colmat==0);

      //4th column is treenodes_bin_arma
      arma::uvec term_nodes=arma::find(treenodes_bin_arma==0);

      term_nodes=term_nodes+1;

      //NumericVector terminal_nodes= wrap(term_nodes);


      //Rcout << "Line 5282.\n";

      //GET J MATRIX

      arma::mat Jmat(num_obs,term_nodes.n_elem);
      arma::mat Jtilde(num_test_obs,term_nodes.n_elem);

      //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
      //NumericVector tree_predictions;

      //now for each internal node find the observations that belong to the terminal nodes

      //NumericVector predictions(test_data.nrow());
      //List term_obs(term_nodes.n_elem);

      //GET J MATRIX

      if(term_nodes.n_elem==1){
        //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
        //predictions=rep(nodemean,test_data.nrow());
        //Rcout << "Line 67 .\n";

        //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
        //term_obs[0]= temp_obsvec;
        //double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

        //double num_prod=1;
        //double num_sum=0;
        //Rcout << "Line 129.\n";
        Jmat.col(0) = arma::ones<arma::vec>(num_obs);

        if(is_test_data==1){
          Jtilde.col(0) = arma::ones<arma::vec>(num_test_obs);
        }

        //for(int k=0; k<num_cats; k++){
        //assuming categories of y are from 1 to num_cats
        //arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
        //double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
        //tree_table1(0,5+k)= m_plus_alph/denom_temp ;

        //for likelihood calculation
        //num_prod=num_prod*tgamma(m_plus_alph);
        //num_sum=num_sum +m_plus_alph ;
        //}

        //lik_prod= alph_term*num_prod/tgamma(num_sum);

      }
      else{
        for(unsigned int i=0;i<term_nodes.n_elem;i++){
          //arma::mat subdata=testd;
          //int curr_term=term_nodes(i);

          int row_index;
          int term_node=term_nodes(i);
          //Rcout << "Line 152.\n";


          //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
          //Why should the ro index be different for a right daughter?
          //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
          row_index=0;

          // if(curr_term % 2==0){
          //   //term node is left daughter
          //   row_index=terminal_nodes[i];
          // }else{
          //   //term node is right daughter
          //   row_index=terminal_nodes[i]-1;
          // }




          //save the left and right node data into arma uvec

          //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
          //arma::vec left_nodes=arma_tree.col(0);
          //arma::vec right_nodes=arma_tree.col(1);

          arma::vec left_nodes=tree_table1.col(0);
          arma::vec right_nodes=tree_table1.col(1);



          arma::mat node_split_mat;
          node_split_mat.set_size(0,3);
          //Rcout << "Line 182. i = " << i << " .\n";

          while(row_index!=1){
            //for each terminal node work backwards and see if the parent node was a left or right node
            //append split info to a matrix
            int rd=0;
            arma::uvec parent_node=arma::find(left_nodes == term_node);

            if(parent_node.size()==0){
              parent_node=arma::find(right_nodes == term_node);
              rd=1;
            }

            //want to cout parent node and append to node_split_mat

            node_split_mat.insert_rows(0,1);

            //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
            //node_split_mat(0,0)=treetable(parent_node[0],2);
            //node_split_mat(0,1)=treetable(parent_node[0],3);

            //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
            //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

            node_split_mat(0,0)=tree_table1(parent_node(0),2);
            node_split_mat(0,1)=tree_table1(parent_node(0),3);

            node_split_mat(0,2)=rd;
            row_index=parent_node(0)+1;
            term_node=parent_node(0)+1;
          }

          //once we have the split info, loop through rows and find the subset indexes for that terminal node!
          //then fill in the predicted value for that tree
          //double prediction = tree_data(term_node,5);
          arma::uvec pred_indices;
          arma::uvec pred_test_indices;
          int split= node_split_mat(0,0)-1;

          //Rcout << "Line 224.\n";
          //Rcout << "split = " << split << ".\n";
          //arma::vec tempvec = testd.col(split);
          arma::vec tempvec = x_control_a.col(split);
          //Rcout << "Line 227.\n";


          double temp_split = node_split_mat(0,1);

          if(node_split_mat(0,2)==0){
            pred_indices = arma::find(tempvec <= temp_split);
          }else{
            pred_indices = arma::find(tempvec > temp_split);
          }

          if(is_test_data==1){
            arma::vec temptest_vec = x_control_test_a.col(split);

            if(node_split_mat(0,2)==0){
              pred_test_indices = arma::find(temptest_vec <= temp_split);
            }else{
              pred_test_indices = arma::find(temptest_vec > temp_split);
            }
          }


          //Rcout << "Line 236.\n";

          arma::uvec temp_pred_indices;
          arma::uvec temp_test_pred_indices;

          //arma::vec data_subset = testd.col(split);
          arma::vec data_subset = x_control_a.col(split);
          data_subset=data_subset.elem(pred_indices);

          arma::vec data_test_subset;
          if(is_test_data==1){
            data_test_subset =x_control_test_a.col(split);
            data_test_subset=data_test_subset.elem(pred_test_indices);
          }

          //now loop through each row of node_split_mat
          int n=node_split_mat.n_rows;
          //Rcout << "Line 174. i = " << i << ". n = " << n << ".\n";
          //Rcout << "Line 248.\n";

          for(int j=1;j<n;j++){
            int curr_sv=node_split_mat(j,0);
            double split_p = node_split_mat(j,1);

            //data_subset = testd.col(curr_sv-1);
            //Rcout << "Line 255.\n";
            //Rcout << "curr_sv = " << curr_sv << ".\n";
            data_subset = x_control_a.col(curr_sv-1);
            //Rcout << "Line 258.\n";

            data_subset=data_subset.elem(pred_indices);


            if(node_split_mat(j,2)==0){
              //split is to the left
              temp_pred_indices=arma::find(data_subset <= split_p);
            }else{
              //split is to the right
              temp_pred_indices=arma::find(data_subset > split_p);
            }
            pred_indices=pred_indices.elem(temp_pred_indices);


            if(is_test_data==1){
              data_test_subset = x_control_test_a.col(curr_sv-1);
              data_test_subset=data_test_subset.elem(pred_test_indices);

              if(node_split_mat(j,2)==0){
                //split is to the left
                temp_test_pred_indices=arma::find(data_test_subset <= split_p);
              }else{
                //split is to the right
                temp_test_pred_indices=arma::find(data_test_subset > split_p);
              }
              pred_test_indices=pred_test_indices.elem(temp_test_pred_indices);

            }


            //if(pred_indices.size()==0){
            //  continue;
            //}

          }
          //Rcout << "Line 199. i = " << i <<  ".\n";

          //There is probably a more efficient way of doing this
          //e.g. initialize J matrix so that all elements are equal to zero
          arma::vec tempcol_J=arma::zeros<arma::vec>(num_obs);
          tempcol_J(pred_indices) = arma::ones<arma::vec>(pred_indices.size());
          Jmat.col(i) = tempcol_J;

          if(is_test_data==1){
            arma::vec tempcol_Jtilde=arma::zeros<arma::vec>(num_test_obs);
            tempcol_Jtilde(pred_test_indices) = arma::ones<arma::vec>(pred_test_indices.size());
            Jtilde.col(i) = tempcol_Jtilde;
          }

          //double nodemean=tree_data(terminal_nodes[i]-1,5);
          //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
          //predictions[predind]= nodemean;
          //term_obs[i]=predind;

          //double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
          //Rcout << "Line 207. predind = " << predind <<  ".\n";
          //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
          // << "Line 207. term_node = " << term_node <<  ".\n";

          //double num_prod=1;
          //double num_sum=0;

          // for(int k=0; k<num_cats; k++){
          //   //assuming categories of y are from 1 to num_cats
          //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
          //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
          //
          //   tree_table1(curr_term-1,5+k)= m_plus_alph/denom_temp ;
          //
          //   num_prod=num_prod*tgamma(m_plus_alph);
          //   num_sum=num_sum +m_plus_alph ;
          // }
          //
          //
          // lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);
          //Rcout << "Line 297.\n";


        }//End of loop over terminal nodes.
      }// end of else statement (for when more than one terminal node)
      // Now have J matrix

      Wmat_mu=join_rows(Wmat_mu,Jmat);
      //or
      //Wmat.insert_cols(Wmat.n_cols,Jmat);
      //or
      //int b_j=term_nodes.n_elem;
      //Wmat.insert_cols(upsilon,Jmat);
      //upsilon+=b_j;


      //Obtain test W_tilde, i.e. W matrix for test data
      if(is_test_data==1){
        W_tilde_mu=join_rows(W_tilde_mu,Jtilde);
      }

      //or
      //W_tilde.insert_cols(W_tilde.n_cols,Jtilde);
      //or
      //int b_jtest=term_nodes.n_elem;
      //W_tilde.insert_cols(upsilon2,Jtilde);
      //upsilon2+=b_jtest;
      //Rcout << "Line 5566.\n";


      if(imp_sampler!=tree_prior){//check if importance sampler is not equal to the prior
        // //get impportance sampler probability and tree prior
        // long double temp_samp_prob;
        // long double temp_prior_prob;
        // //get sampler tree probability
        // if(imp_sampler==1){//If sample from BART prior
        //
        //
        //
        //   temp_samp_prob=1;
        //
        //   double depth1=0;
        //   int prev_node=0; //1 if previous node splits, zero otherwise
        //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
        //     if(treenodes_bin[i_2]==1){
        //       temp_samp_prob=temp_samp_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
        //       depth1=depth1+1; //after a split, the depth will increase by 1
        //       prev_node=1;
        //     }else{
        //       temp_samp_prob=temp_samp_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
        //       if(prev_node==1){//zero following a 1, therefore at same depth.
        //         //Don't change depth. Do nothing
        //       }else{ //zero following a zero, therefore the depth will decrease by 1
        //         depth1=depth1-1;
        //       }
        //       prev_node=0;
        //
        //     }
        //   }
        //
        //   //end of calculating BART tree probability
        // }else{
        //   if(imp_sampler==2){//If sample from spike and tree prior
        //     throw std::range_error("code not yet written for spike and tree prior");
        //
        //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
        //     double tempexp1=treenodes_bin.size()-arma::sum(treenodes_bin_arma);
        //     double tempexp2=arma::sum(treenodes_bin_arma);
        //     temp_samp_prob=pow(lambda,tempexp2)*
        //       pow(1-lambda,tempexp1);
        //       //(1/pow(double(num_split_vars),tempexp2));
        //
        //       temp_samp_prob=exp(log(lambda)*tempexp2+
        //         log(1-lambda)*tempexp1);
        //
        //     //temp_samp_prob=pow(lambda,arma::sum(treenodes_bin_arma))*
        //     //  pow(1-lambda,treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
        //     //  pow((1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
        //   }
        // }
        //
        // sum_tree_samp_prob=sum_tree_samp_prob*temp_samp_prob;
        // //end of getting importance sampler probability
        //
        // //get prior tree probability
        // if(tree_prior==1){//If sample from BART prior
        //
        //
        //
        //   temp_prior_prob=1;
        //
        //   double depth1=0;
        //   int prev_node=0; //1 if previous node splits, zero otherwise
        //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
        //
        //     if(treenodes_bin[i_2]==1){
        //       temp_prior_prob=temp_prior_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
        //       depth1=depth1+1; //after a split, the depth will increase by 1
        //       prev_node=1;
        //     }else{
        //       temp_prior_prob=temp_prior_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
        //       if(prev_node==1){//zero following a 1, therefore at same depth.
        //         //Don't change depth. Do nothing
        //       }else{ //zero following a zero, therefore the depth will decrease by 1
        //         depth1=depth1-1;
        //       }
        //       prev_node=0;
        //
        //     }
        //     //if(alpha_BART==0){
        //     //  Rcout << "alpha_BART equals zero!!!!.\n";
        //     //}
        //   }
        //
        //   //end of calculating BART tree probability
        // }else{
        //   if(tree_prior==2){//If sample from spike and tree prior
        //     throw std::range_error("code not yet written for spike and tree prior");
        //
        //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
        //     temp_prior_prob=pow((long double)(lambda),arma::sum(treenodes_bin_arma))*
        //       pow((long double)(1-lambda),treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
        //       pow((long double)(1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
        //   }
        // }
        //
        // sum_tree_prior_prob=sum_tree_prior_prob*temp_prior_prob;
        // if(temp_prior_prob==0){
        //   Rcout << "Line 4097, j= " << j << ". \n";
        //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
        //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
        // }
        // if(temp_samp_prob==0){
        //   Rcout << "Line 4102, j= " << j << ". \n";
        //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
        //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
        // }
        //
        // if(sum_tree_samp_prob==0){
        //   Rcout << "Line 4102, j= " << j << ". \n";
        //   Rcout << "sum_tree_samp_prob= " << sum_tree_samp_prob << ". \n";
        //   //Rcout << "treenodes_bin_arma= " << treenodes_bin_arma << ". \n";
        //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
        //
        // }




        //get tree prior over impportance sampler probability
        double tree_prior_over_samp_prob=1;
        if(imp_sampler==1){   //If sample from BART prior


          if(tree_prior==1){  //If tree prior is BART prior
            throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

          }else{
            if(tree_prior==2){  //If tree prior is spike-and-tree prior
              throw std::range_error("code not yet written for spike and tree prior");

            }else{//otherwise the tree prior is the Quadrianto and Ghahramani prior
              double depth1=0;
              int prev_node=0; //1 if previous node splits, zero otherwise
              for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
                if(treenodes_bin[i_2]==1){
                  tree_prior_over_samp_prob=tree_prior_over_samp_prob*(lambda_mu/(alpha_BCF_mu*pow(double(depth1+1),-beta_BCF_mu)));
                  depth1=depth1+1; //after a split, the depth will increase by 1
                  prev_node=1;
                }else{
                  tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-lambda_mu)/(1-alpha_BCF_mu*pow(double(depth1+1),-beta_BCF_mu)));
                  if(prev_node==1){//zero following a 1, therefore at same depth.
                    //Don't change depth. Do nothing
                  }else{ //zero following a zero, therefore the depth will decrease by 1
                    depth1=depth1-1;
                  }
                  prev_node=0;

                }
              }

            }
          }


        }else{// if not sampling from BART prior
          if(imp_sampler==2){//If sample from spike and tree prior
            throw std::range_error("code not yet written for sampling from spike and tree prior");

          }else{//otherwise sampling from Quadrianto and Ghahramani prior
            if(tree_prior==1){  //If tree prior is BART prior

              double depth1=0;
              int prev_node=0; //1 if previous node splits, zero otherwise
              for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
                if(treenodes_bin[i_2]==1){
                  tree_prior_over_samp_prob=tree_prior_over_samp_prob*((alpha_BCF_mu*pow(double(depth1+1),-beta_BCF_mu))/lambda_mu);
                  depth1=depth1+1; //after a split, the depth will increase by 1
                  prev_node=1;
                }else{
                  tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-alpha_BCF_mu*pow(double(depth1+1),-beta_BCF_mu))/(1-lambda_mu));
                  if(prev_node==1){//zero following a 1, therefore at same depth.
                    //Don't change depth. Do nothing
                  }else{ //zero following a zero, therefore the depth will decrease by 1
                    depth1=depth1-1;
                  }
                  prev_node=0;

                }//close (zero node) else stattement

              }//end for loop over i_2

            }else{
              if(tree_prior==2){  //If tree prior is spike-and-tree prior
                throw std::range_error("code not yet written for spike and tree prior");

              }else{
                throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

              }//close (not BART nor spike and tree prior) else statement
            }// close (not BART prior) else statememt




          }//close (not sampling from BART or spike and tree)  else statement

        }//close (not sampling from BART) else statement

        sum_prior_over_samp_prob=sum_prior_over_samp_prob*tree_prior_over_samp_prob;
        //end of getting tree prior over impportance sampler probability



      }//end of tree prior and importance sampler calculations




    } //end of loop over mu trees in sum


    /////////////////////////////////////////////////////////////////////////////////////////
    //Rcout << "Line 5782. TAU TREES. \n";

    /////////////////////////////////////////////////////////////////////////////////////////
    // NOW LOOP OVER TAU TREES


 for(int q=0; q<num_trees_tau;q++){  //start of loop over trees in sum


   //If parallelizing, define the distributinos before this loop
   //and use lrng and the following two lines
   //dqrng::xoshiro256plus lrng(rng);      // make thread local copy of rng
   //lrng.jump(omp_get_thread_num() + 1);  // advance rng by 1 ... nthreads jumps


   //NumericVector treenodes_bin(0);
   //arma::uvec treenodes_bin(0);

   std::vector<int> treenodes_bin;


   int count_terminals = 0;
   int count_internals = 0;

   //int count_treebuild = 0;


   if(imp_sampler==1){ //If sampling from BART prior

     double depth1=0;
     int prev_node=0; //1 if previous node splits, zero otherwise

     double samp_prob;

     while(count_internals > (count_terminals -1)){
       samp_prob=alpha_BCF_tau*pow(double(depth1+1),-beta_BCF_tau);
       std::bernoulli_distribution coin_flip2(samp_prob);

       int tempdraw = coin_flip2(lgen);
       treenodes_bin.push_back(tempdraw);

       if(tempdraw==1){

         depth1=depth1+1; //after a split, the depth will increase by 1
         prev_node=1;
         count_internals=count_internals+1;

       }else{

         if(prev_node==1){//zero following a 1, therefore at same depth.
           //Don't change depth. Do nothing
         }else{ //zero following a zero, therefore the depth will decrease by 1
           depth1=depth1-1;
         }
         prev_node=0;
         count_terminals=count_terminals+1;

       }

     }

   }else{  //If not sampling from BART prior
     if(imp_sampler==2){//If sampling from spike and tree prior
       throw std::range_error("code not yet written for spike and tree sampling");

     }else{//If sampling from default Q+G prior. i.e. not sampling from BART nor spike and tree prior

       while(count_internals > (count_terminals -1)){

         //Also consider standard library and random header
         // std::random_device device;
         // std::mt19937 gen(device());
         // std::bernoulli_distribution coin_flip(lambda);
         // bool outcome = coin_flip(gen);


         int tempdraw = coin_flip_tau(lgen);

         //int tempdraw = rbinom(n = 1, prob = lambda,size=1);


         //int tempdraw = Rcpp::rbinom(1,lambda,1);
         //int tempdraw = R::rbinom(1,lambda);

         ////Rcout << "tempdraw = " << tempdraw << ".\n" ;

         //int tempdraw = coin_flip2(lgen)-1;

         //int tempdraw = dqrng::dqsample_int(2, 1, true,lambdavec )-1;


         //need to update rng if use boost?
         //int tempdraw = bernoulli(rng, binomial::param_type(1, lambda));

         treenodes_bin.push_back(tempdraw);


         if(tempdraw==1){
           count_internals=count_internals+1;
         }else{
           count_terminals=count_terminals+1;
         }

       }//end of while loop creating parent vector treenodes_bin
     }

   }



   //Consider making this an armadillo vector
   //IntegerVector split_var_vec(treenodes_bin.size());
   //arma::uvec split_var_vec(treenodes_bin.size());
   std::vector<int> split_var_vec(treenodes_bin.size());

   //loop drawing splitting variables
   //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

   //if using armadillo, it might be faster to subset to split nodes
   //then use a vector of draws
   for(unsigned int i=0; i<treenodes_bin.size();i++){
     if(treenodes_bin[i]==0){
       split_var_vec[i] = -1;
     }else{
       // also consider the standard library function uniform_int_distribution
       // might need random header
       // This uses the Mersenne twister

       //Three lines below should probably be outside all the loops
       // std::random_device rd;
       // std::mt19937 engine(rd());
       // std::uniform_int_distribution<> distsampvar(1, num_split_vars);
       //
       // split_var_vec[i] = distsampvar(engine);

       split_var_vec[i] = distsampvar_tau(lgen);


       //consider using boost
       //might need to update rng
       //split_var_vec[i] <- sample_splitvars(rng);

       //or use dqrng
       //not sure if have to update the random number
       //check if the following line is written properly
       //split_var_vec[i] = dqrng::dqsample_int(num_split_vars, 1, true);

       //not sure if this returns an integer or a vector?
       //split_var_vec[i] = RcppArmadillo::sample(num_split_vars, 1,true);
       //could try
       //split_var_vec[i] = as<int>(Rcpp::sample(num_split_vars, 1,true));
       //could also try RcppArmadillo::rmultinom

     }

   }// end of for-loop drawing split variables


   //Consider making this an armadillo vector
   //NumericVector split_point_vec(treenodes_bin.size());
   //arma::vec split_point_vec(treenodes_bin.size());
   std::vector<double> split_point_vec(treenodes_bin.size());


   //loop drawing splitting points
   //REPLACE SQUARE BRACKETS WITH "( )" if using ARMADILLO vector for split_var_vec or treenodes_bin

   //if using armadillo, it might be faster to subset to split nodes
   //then use a vector of draws
   for(unsigned int i=0; i<treenodes_bin.size();i++){
     if(treenodes_bin[i]==0){
       split_point_vec[i] = -1;
     }else{


       //////////////////////////////////////////////////////////
       //following function not reccommended
       //split_point_vec[i] = std::rand();
       //////////////////////////////////////////////////////////
       ////Standard library:
       ////This should probably be outside all the loops
       ////std::random_device rd;  //Will be used to obtain a seed for the random number engine
       ////std::mt19937 gen2(rd()); //Standard mersenne_twister_engine seeded with rd()
       ////std::uniform_real_distribution<> dis_cont_unif(0, 1);

       split_point_vec[i] = dis_cont_unif(lgen);

       //////////////////////////////////////////////////////////
       //from armadillo
       //split_point_vec[i] = arma::randu();

       //////////////////////////////////////////////////////////
       //probably not adviseable for paralelization
       //From Rcpp
       //split_point_vec[i] = as<double>(Rcpp::runif(1,0,1));

       //////////////////////////////////////////////////////////
       //consider using boost
       //might need to update rng
       //split_point_vec[i] <- b_unif_point(rng);

       //or use dqrng
       //not sure if have to update the random number
       //check if the following line is written properly
       //split_point_vec[i] = dqrng::dqrunif(1, 0, 1);

       //not sure if this returns an integer or a vector?





     }

   }// end of for-loop drawing split points



   //Rcout << "Line 6000.\n";


   //CODE FOR ADJUSTING SPLITTING POINTS SO THAT THE TREES ARE VALID
   if(valid_trees==1){
     for(unsigned int i=0; i<treenodes_bin.size();i++){ //loop over all nodes
       if(treenodes_bin[i]==1){ // if it is an internal node, then check for further splits on the same variable and update
         double first_split_var=split_var_vec[i];      //splitting variable to check for
         double first_split_point=split_point_vec[i];  //splitting point to use in updates

         double sub_int_nodes=0;       //this internal node count will be used to determine if in subtree relevant to sub_int_nodes
         double sub_term_nodes=0;      //this terminal node count will be used to determine if in subtree relevant to sub_int_nodes
         double preventing_updates=0; //indicates if still within subtree that is not to be updated
         double prevent_int_count=0;   //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
         double prevent_term_count=0;  //this internal node count will be used to determine if in sub-sub-tree that is not to be updated within the k loop
         for(unsigned int k=i+1; k<treenodes_bin.size();k++){
           if(treenodes_bin[k]==1){
             sub_int_nodes=sub_int_nodes+1;
             if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
               prevent_int_count=prevent_int_count+1;
             }
           }else{
             sub_term_nodes=sub_term_nodes+1;
             if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
               prevent_term_count=prevent_term_count+1;
             }
           }
           if(sub_int_nodes<=sub_term_nodes-2){
             break;
           }


           if(preventing_updates==1){ //indicates if still within subtree that is not to be updated
             if(prevent_int_count>prevent_term_count-1){ //if this rule is satisfied then in subtree that is not to be updated
               continue; //still in subtree, therefore continue instead of checking for splits to be updates
             }else{
               preventing_updates=0; // no longer in subtree, therefore reset preventing_updates to zero
             }
           }


           if(sub_int_nodes>sub_term_nodes-1){
             if(treenodes_bin[k]==1){
               if(split_var_vec[k]==first_split_var){
                 split_point_vec[k]=split_point_vec[k]*first_split_point;
                 //beginning count of subtree that should not have
                 //further splits on first_split_var updated
                 preventing_updates=1; //indicates if still within subtree that is not to be updated
                 prevent_int_count=1;
                 prevent_term_count=0;
               }
             }
           }else{
             if(treenodes_bin[k]==1){
               if(split_var_vec[k]==first_split_var){
                 split_point_vec[k]=split_point_vec[k]+first_split_point-first_split_point*split_point_vec[k];
                 //beginning count of subtree that should not have
                 //further splits on first_split_var updated
                 preventing_updates=1; //indicates if still within subtree that is not to be updated
                 prevent_int_count=1;
                 prevent_term_count=0;
               }
             }
           }



         }//end of inner loop over k
       }//end of if statement treenodes_bin[i]==1)
     }//end of loop over i
   }//end of if statement valid_trees==1







   //Rcout << "Line 6078.\n";



   //Create tree table matrix

   //NumericMatrix tree_table1(treenodes_bin.size(),5+num_cats);

   //Rcout << "Line 1037. \n";
   //arma::mat tree_table1(treenodes_bin.size(),5+num_cats);

   //initialize with zeros. Not sure if this is necessary
   arma::mat tree_table1=arma::zeros<arma::mat>(treenodes_bin.size(),6);
   //Rcout << "Line 1040. \n";


   //tree_table1(_,2) = wrap(split_var_vec);
   //tree_table1(_,3) = wrap(split_point_vec);
   //tree_table1(_,4) = wrap(treenodes_bin);

   //It might be more efficient to make everything an armadillo object initially
   // but then would need to replace push_back etc with a different approach (but this might be more efficient anyway)
   arma::colvec split_var_vec_arma=arma::conv_to<arma::colvec>::from(split_var_vec);
   arma::colvec split_point_vec_arma(split_point_vec);
   arma::colvec treenodes_bin_arma=arma::conv_to<arma::colvec>::from(treenodes_bin);


   //Rcout << "Line 1054. \n";

   //Fill in splitting variable column
   tree_table1.col(2) = split_var_vec_arma;
   //Fill in splitting point column
   tree_table1.col(3) = split_point_vec_arma;
   //Fill in split/parent column
   tree_table1.col(4) = treenodes_bin_arma;


   //Rcout << "Line 1061. j = " << j << ". \n";
   //Rcout << "Line 6117. j = " << j << ". \n";
   //Rcout << "Line 6118. tree_table1 tau = " << tree_table1 << ". \n";



   // Now start filling in left daughter and right daughter columns
   std::vector<int> rd_spaces;
   int prev_node = -1;

   for(unsigned int i=0; i<treenodes_bin.size();i++){
     //Rcout << "Line 1061. i = " << i << ". \n";
     if(prev_node==0){
       //tree_table1(rd_spaces[rd_spaces.size()-1], 1)=i;
       //Rcout << "Line 1073. j = " << j << ". \n";

       tree_table1(rd_spaces.back(), 1)=i+1;
       //Rcout << "Line 1076. j = " << j << ". \n";

       rd_spaces.pop_back();
     }
     if(treenodes_bin[i]==1){
       //Rcout << "Line 1081. j = " << j << ". \n";

       tree_table1(i,0) = i+2;
       rd_spaces.push_back(i);
       prev_node = 1;
       //Rcout << "Line 185. j = " << j << ". \n";

     }else{                  // These 2 lines unnecessary if begin with matrix of zeros
       //Rcout << "Line 1089. j = " << j << ". \n";
       tree_table1(i,0)=0 ;
       tree_table1(i,1) = 0 ;
       prev_node = 0;
       //Rcout << "Line 1093. j = " << j << ". \n";

     }
   }//
   //Rcout << "Line 1097. j = " << j << ". \n";





   //List treepred_output = get_treepreds(original_y, num_cats, alpha_pars,
   //                                     originaldata,
   //                                     treetable_list[i]  );


   //use armadillo object tree_table1

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //create variables for likelihood calcuations
   // double lik_prod=1;
   // double alph_prod=1;
   // for(unsigned int i=0; i<alpha_pars_arma.n_elem;i++){
   //   alph_prod=alph_prod*tgamma(alpha_pars_arma(i));
   // }
   // double gam_alph_sum= tgamma(arma::sum(alpha_pars_arma));
   // double alph_term=gam_alph_sum/alph_prod;

   //arma::mat arma_tree_table(treetable.begin(), treetable.nrow(), treetable.ncol(), false);
   //arma::mat arma_orig_data(originaldata.begin(), originaldata.nrow(), originaldata.ncol(), false);


   //arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
   //arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);

   //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);

   //NumericVector terminal_nodes=find_term_nodes(treetable);

   //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);

   //arma::vec colmat=arma_tree.col(4);
   //arma::uvec term_nodes=arma::find(colmat==-1);

   //arma::vec colmat=arma_tree.col(2);
   //arma::uvec term_nodes=arma::find(colmat==0);

   //arma::vec colmat=tree_table1.col(4);
   //arma::uvec term_nodes=arma::find(colmat==0);

   //4th column is treenodes_bin_arma
   arma::uvec term_nodes=arma::find(treenodes_bin_arma==0);

   term_nodes=term_nodes+1;

   //NumericVector terminal_nodes= wrap(term_nodes);

   //Rcout << "Line 6207.\n";


   //GET J MATRIX

   arma::mat Jmat(num_obs,term_nodes.n_elem);
   arma::mat Jtilde(num_test_obs,term_nodes.n_elem);

   //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
   //NumericVector tree_predictions;

   //now for each internal node find the observations that belong to the terminal nodes

   //NumericVector predictions(test_data.nrow());
   //List term_obs(term_nodes.n_elem);

   //GET J MATRIX

   if(term_nodes.n_elem==1){
     //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
     //predictions=rep(nodemean,test_data.nrow());
     //Rcout << "Line 67 .\n";

     //IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
     //term_obs[0]= temp_obsvec;
     //double denom_temp= orig_y_arma.n_elem+arma::sum(alpha_pars_arma);

     //double num_prod=1;
     //double num_sum=0;
     //Rcout << "Line 129.\n";
     Jmat.col(0) = arma::ones<arma::vec>(num_obs);

     if(is_test_data==1){
       Jtilde.col(0) = arma::ones<arma::vec>(num_test_obs);
     }

     //for(int k=0; k<num_cats; k++){
     //assuming categories of y are from 1 to num_cats
     //arma::uvec cat_inds= arma::find(orig_y_arma==k+1);
     //double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
     //tree_table1(0,5+k)= m_plus_alph/denom_temp ;

     //for likelihood calculation
     //num_prod=num_prod*tgamma(m_plus_alph);
     //num_sum=num_sum +m_plus_alph ;
     //}

     //lik_prod= alph_term*num_prod/tgamma(num_sum);

   }
   else{
     for(unsigned int i=0;i<term_nodes.n_elem;i++){
       //arma::mat subdata=testd;
       //int curr_term=term_nodes(i);

       int row_index;
       int term_node=term_nodes(i);
       //Rcout << "Line 152.\n";


       //WHAT IS THE PURPOSE OF THIS IF-STATEMENT?
       //Why should the ro index be different for a right daughter?
       //Why not just initialize row_index to any number not equal to 1 (e.g. 0)?
       row_index=0;

       // if(curr_term % 2==0){
       //   //term node is left daughter
       //   row_index=terminal_nodes[i];
       // }else{
       //   //term node is right daughter
       //   row_index=terminal_nodes[i]-1;
       // }




       //save the left and right node data into arma uvec

       //CHECK THAT THIS REFERS TO THE CORRECT COLUMNS
       //arma::vec left_nodes=arma_tree.col(0);
       //arma::vec right_nodes=arma_tree.col(1);

       arma::vec left_nodes=tree_table1.col(0);
       arma::vec right_nodes=tree_table1.col(1);



       arma::mat node_split_mat;
       node_split_mat.set_size(0,3);
       //Rcout << "Line 6296. i = " << i << " .\n";

       while(row_index!=1){
         //for each terminal node work backwards and see if the parent node was a left or right node
         //append split info to a matrix
         int rd=0;
         arma::uvec parent_node=arma::find(left_nodes == term_node);

         if(parent_node.size()==0){
           parent_node=arma::find(right_nodes == term_node);
           rd=1;
         }

         //want to cout parent node and append to node_split_mat

         node_split_mat.insert_rows(0,1);

         //CHECK THAT COLUMNS OF TREETABLE ARE CORRECT
         //node_split_mat(0,0)=treetable(parent_node[0],2);
         //node_split_mat(0,1)=treetable(parent_node[0],3);

         //node_split_mat(0,0)=arma_tree_table(parent_node[0],3);
         //node_split_mat(0,1)=arma_tree_table(parent_node[0],4);

         node_split_mat(0,0)=tree_table1(parent_node(0),2);
         node_split_mat(0,1)=tree_table1(parent_node(0),3);

         node_split_mat(0,2)=rd;
         row_index=parent_node(0)+1;
         term_node=parent_node(0)+1;
       }

       //once we have the split info, loop through rows and find the subset indexes for that terminal node!
       //then fill in the predicted value for that tree
       //double prediction = tree_data(term_node,5);
       arma::uvec pred_indices;
       arma::uvec pred_test_indices;
       int split= node_split_mat(0,0)-1;

       //Rcout << "Line 6335.\n";
       //Rcout << "split = " << split << ".\n";
       //Rcout << "x_moderate_a.n_cols = " << x_moderate_a.n_cols << ".\n";


       //arma::vec tempvec = testd.col(split);
       arma::vec tempvec = x_moderate_a.col(split);
       ////Rcout << "Line 227.\n";

       //Rcout << "Line 6341.\n";

       double temp_split = node_split_mat(0,1);

       if(node_split_mat(0,2)==0){
         pred_indices = arma::find(tempvec <= temp_split);
       }else{
         pred_indices = arma::find(tempvec > temp_split);
       }

       //Rcout << "Line 6351.\n";

       if(is_test_data==1){
         arma::vec temptest_vec = x_moderate_test_a.col(split);
         //Rcout << "Line 6355.\n";

         if(node_split_mat(0,2)==0){
           pred_test_indices = arma::find(temptest_vec <= temp_split);
         }else{
           pred_test_indices = arma::find(temptest_vec > temp_split);
         }
       }


       //Rcout << "Line 6361.\n";

       arma::uvec temp_pred_indices;
       arma::uvec temp_test_pred_indices;

       //arma::vec data_subset = testd.col(split);
       arma::vec data_subset = x_moderate_a.col(split);
       data_subset=data_subset.elem(pred_indices);

       arma::vec data_test_subset;
       if(is_test_data==1){
         data_test_subset =x_moderate_test_a.col(split);
         data_test_subset=data_test_subset.elem(pred_test_indices);
       }

       //now loop through each row of node_split_mat
       int n=node_split_mat.n_rows;

       //Rcout << "Line 6378. i = " << i << ". n = " << n << ".\n";

       for(int j=1;j<n;j++){
         int curr_sv=node_split_mat(j,0);
         double split_p = node_split_mat(j,1);

         //data_subset = testd.col(curr_sv-1);
         //Rcout << "Line 255.\n";
         //Rcout << "curr_sv = " << curr_sv << ".\n";
         data_subset = x_moderate_a.col(curr_sv-1);
         //Rcout << "Line 258.\n";

         data_subset=data_subset.elem(pred_indices);


         if(node_split_mat(j,2)==0){
           //split is to the left
           temp_pred_indices=arma::find(data_subset <= split_p);
         }else{
           //split is to the right
           temp_pred_indices=arma::find(data_subset > split_p);
         }
         pred_indices=pred_indices.elem(temp_pred_indices);


         if(is_test_data==1){
           data_test_subset = x_moderate_test_a.col(curr_sv-1);
           data_test_subset=data_test_subset.elem(pred_test_indices);

           if(node_split_mat(j,2)==0){
             //split is to the left
             temp_test_pred_indices=arma::find(data_test_subset <= split_p);
           }else{
             //split is to the right
             temp_test_pred_indices=arma::find(data_test_subset > split_p);
           }
           pred_test_indices=pred_test_indices.elem(temp_test_pred_indices);

         }


         //if(pred_indices.size()==0){
         //  continue;
         //}

       }
       //Rcout << "Line 6425. i = " << i <<  ".\n";

       //There is probably a more efficient way of doing this
       //e.g. initialize J matrix so that all elements are equal to zero
       arma::vec tempcol_J=arma::zeros<arma::vec>(num_obs);
       tempcol_J(pred_indices) = arma::ones<arma::vec>(pred_indices.size());
       Jmat.col(i) = tempcol_J;

       if(is_test_data==1){
         arma::vec tempcol_Jtilde=arma::zeros<arma::vec>(num_test_obs);
         tempcol_Jtilde(pred_test_indices) = arma::ones<arma::vec>(pred_test_indices.size());
         Jtilde.col(i) = tempcol_Jtilde;
       }

       //double nodemean=tree_data(terminal_nodes[i]-1,5);
       //IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
       //predictions[predind]= nodemean;
       //term_obs[i]=predind;

       //double denom_temp= pred_indices.n_elem+arma::sum(alpha_pars_arma);
       //Rcout << "Line 207. predind = " << predind <<  ".\n";
       //Rcout << "Line 207. denom_temp = " << denom_temp <<  ".\n";
       // << "Line 207. term_node = " << term_node <<  ".\n";

       //double num_prod=1;
       //double num_sum=0;

       // for(int k=0; k<num_cats; k++){
       //   //assuming categories of y are from 1 to num_cats
       //   arma::uvec cat_inds= arma::find(orig_y_arma(pred_indices)==k+1);
       //   double m_plus_alph=cat_inds.n_elem +alpha_pars_arma(k);
       //
       //   tree_table1(curr_term-1,5+k)= m_plus_alph/denom_temp ;
       //
       //   num_prod=num_prod*tgamma(m_plus_alph);
       //   num_sum=num_sum +m_plus_alph ;
       // }
       //
       //
       // lik_prod= lik_prod*alph_term*num_prod/tgamma(num_sum);

       //Rcout << "Line 6466.\n";


     }//End of loop over terminal nodes.
   }// end of else statement (for when more than one terminal node)
   // Now have J matrix

   //Rcout << "Line 6472.\n";

   Wmat_tau=join_rows(Wmat_tau,Jmat);
   //or
   //Wmat.insert_cols(Wmat.n_cols,Jmat);
   //or
   //int b_j=term_nodes.n_elem;
   //Wmat.insert_cols(upsilon,Jmat);
   //upsilon+=b_j;


   //Obtain test W_tilde, i.e. W matrix for test data
   if(is_test_data==1){
     W_tilde_tau=join_rows(W_tilde_tau,Jtilde);
   }

   //or
   //W_tilde.insert_cols(W_tilde.n_cols,Jtilde);
   //or
   //int b_jtest=term_nodes.n_elem;
   //W_tilde.insert_cols(upsilon2,Jtilde);
   //upsilon2+=b_jtest;


   if(imp_sampler!=tree_prior){//check if importance sampler is not equal to the prior
     // //get impportance sampler probability and tree prior
     // long double temp_samp_prob;
     // long double temp_prior_prob;
     // //get sampler tree probability
     // if(imp_sampler==1){//If sample from BART prior
     //
     //
     //
     //   temp_samp_prob=1;
     //
     //   double depth1=0;
     //   int prev_node=0; //1 if previous node splits, zero otherwise
     //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
     //     if(treenodes_bin[i_2]==1){
     //       temp_samp_prob=temp_samp_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
     //       depth1=depth1+1; //after a split, the depth will increase by 1
     //       prev_node=1;
     //     }else{
     //       temp_samp_prob=temp_samp_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
     //       if(prev_node==1){//zero following a 1, therefore at same depth.
     //         //Don't change depth. Do nothing
     //       }else{ //zero following a zero, therefore the depth will decrease by 1
     //         depth1=depth1-1;
     //       }
     //       prev_node=0;
     //
     //     }
     //   }
     //
     //   //end of calculating BART tree probability
     // }else{
     //   if(imp_sampler==2){//If sample from spike and tree prior
     //     throw std::range_error("code not yet written for spike and tree prior");
     //
     //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
     //     double tempexp1=treenodes_bin.size()-arma::sum(treenodes_bin_arma);
     //     double tempexp2=arma::sum(treenodes_bin_arma);
     //     temp_samp_prob=pow(lambda,tempexp2)*
     //       pow(1-lambda,tempexp1);
     //       //(1/pow(double(num_split_vars),tempexp2));
     //
     //       temp_samp_prob=exp(log(lambda)*tempexp2+
     //         log(1-lambda)*tempexp1);
     //
     //     //temp_samp_prob=pow(lambda,arma::sum(treenodes_bin_arma))*
     //     //  pow(1-lambda,treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
     //     //  pow((1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
     //   }
     // }
     //
     // sum_tree_samp_prob=sum_tree_samp_prob*temp_samp_prob;
     // //end of getting importance sampler probability
     //
     // //get prior tree probability
     // if(tree_prior==1){//If sample from BART prior
     //
     //
     //
     //   temp_prior_prob=1;
     //
     //   double depth1=0;
     //   int prev_node=0; //1 if previous node splits, zero otherwise
     //   for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
     //
     //     if(treenodes_bin[i_2]==1){
     //       temp_prior_prob=temp_prior_prob*alpha_BART*pow(double(depth1+1),-beta_BART);
     //       depth1=depth1+1; //after a split, the depth will increase by 1
     //       prev_node=1;
     //     }else{
     //       temp_prior_prob=temp_prior_prob*(1-alpha_BART*pow(double(depth1+1),-beta_BART));
     //       if(prev_node==1){//zero following a 1, therefore at same depth.
     //         //Don't change depth. Do nothing
     //       }else{ //zero following a zero, therefore the depth will decrease by 1
     //         depth1=depth1-1;
     //       }
     //       prev_node=0;
     //
     //     }
     //     //if(alpha_BART==0){
     //     //  Rcout << "alpha_BART equals zero!!!!.\n";
     //     //}
     //   }
     //
     //   //end of calculating BART tree probability
     // }else{
     //   if(tree_prior==2){//If sample from spike and tree prior
     //     throw std::range_error("code not yet written for spike and tree prior");
     //
     //   }else{//otherwise sampling from Quadrianto and Ghahramani prior
     //     temp_prior_prob=pow((long double)(lambda),arma::sum(treenodes_bin_arma))*
     //       pow((long double)(1-lambda),treenodes_bin.size()-arma::sum(treenodes_bin_arma))*
     //       pow((long double)(1/double(num_split_vars)),arma::sum(treenodes_bin_arma));
     //   }
     // }
     //
     // sum_tree_prior_prob=sum_tree_prior_prob*temp_prior_prob;
     // if(temp_prior_prob==0){
     //   Rcout << "Line 4097, j= " << j << ". \n";
     //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
     //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
     // }
     // if(temp_samp_prob==0){
     //   Rcout << "Line 4102, j= " << j << ". \n";
     //   Rcout << "temp_prior_prob= " << temp_prior_prob << ". \n";
     //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
     // }
     //
     // if(sum_tree_samp_prob==0){
     //   Rcout << "Line 4102, j= " << j << ". \n";
     //   Rcout << "sum_tree_samp_prob= " << sum_tree_samp_prob << ". \n";
     //   //Rcout << "treenodes_bin_arma= " << treenodes_bin_arma << ". \n";
     //   Rcout << "temp_samp_prob= " << temp_samp_prob << ". \n";
     //
     // }




     //get tree prior over impportance sampler probability
     double tree_prior_over_samp_prob=1;
     if(imp_sampler==1){   //If sample from BART prior


       if(tree_prior==1){  //If tree prior is BART prior
         throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

       }else{
         if(tree_prior==2){  //If tree prior is spike-and-tree prior
           throw std::range_error("code not yet written for spike and tree prior");

         }else{//otherwise the tree prior is the Quadrianto and Ghahramani prior
           double depth1=0;
           int prev_node=0; //1 if previous node splits, zero otherwise
           for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
             if(treenodes_bin[i_2]==1){
               tree_prior_over_samp_prob=tree_prior_over_samp_prob*(lambda_tau/(alpha_BCF_tau*pow(double(depth1+1),-beta_BCF_tau)));
               depth1=depth1+1; //after a split, the depth will increase by 1
               prev_node=1;
             }else{
               tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-lambda_tau)/(1-alpha_BCF_tau*pow(double(depth1+1),-beta_BCF_tau)));
               if(prev_node==1){//zero following a 1, therefore at same depth.
                 //Don't change depth. Do nothing
               }else{ //zero following a zero, therefore the depth will decrease by 1
                 depth1=depth1-1;
               }
               prev_node=0;

             }
           }

         }
       }


     }else{// if not sampling from BART prior
       if(imp_sampler==2){//If sample from spike and tree prior
         throw std::range_error("code not yet written for sampling from spike and tree prior");

       }else{//otherwise sampling from Quadrianto and Ghahramani prior
         if(tree_prior==1){  //If tree prior is BART prior

           double depth1=0;
           int prev_node=0; //1 if previous node splits, zero otherwise
           for(unsigned int i_2=0; i_2<treenodes_bin.size();i_2++){
             if(treenodes_bin[i_2]==1){
               tree_prior_over_samp_prob=tree_prior_over_samp_prob*((alpha_BCF_tau*pow(double(depth1+1),-beta_BCF_tau))/lambda_tau);
               depth1=depth1+1; //after a split, the depth will increase by 1
               prev_node=1;
             }else{
               tree_prior_over_samp_prob=tree_prior_over_samp_prob*((1-alpha_BCF_tau*pow(double(depth1+1),-beta_BCF_tau))/(1-lambda_tau));
               if(prev_node==1){//zero following a 1, therefore at same depth.
                 //Don't change depth. Do nothing
               }else{ //zero following a zero, therefore the depth will decrease by 1
                 depth1=depth1-1;
               }
               prev_node=0;

             }//close (zero node) else stattement

           }//end for loop over i_2

         }else{
           if(tree_prior==2){  //If tree prior is spike-and-tree prior
             throw std::range_error("code not yet written for spike and tree prior");

           }else{
             throw std::range_error("The code should not calculate the ratio of probabilities if sampler equals prior");

           }//close (not BART nor spike and tree prior) else statement
         }// close (not BART prior) else statememt




       }//close (not sampling from BART or spike and tree)  else statement

     }//close (not sampling from BART) else statement

     sum_prior_over_samp_prob=sum_prior_over_samp_prob*tree_prior_over_samp_prob;
     //end of getting tree prior over impportance sampler probability



   }//end of tree prior and importance sampler calculations



 } //end of loop over trees in sum



 //Rcout << "6708 4528.\n";


      double b_mu=Wmat_mu.n_cols;
      double b_tau=Wmat_tau.n_cols;
      //Wmat_tau.each_col()%=z_ar;


      //Rcout << "Line 14688.\n";
      //Rcout << "b_tau = " << b_tau << ".\n";
      //Rcout << "Wmat_tau.n_rows = " << Wmat_tau.n_rows << ".\n";

      //Rcout << "z_ar.n_elem = " << z_ar.n_elem << ".\n";

      arma::mat DiagZ_Wmat_tau= Wmat_tau.each_col()%z_ar;
      //Rcout << "Line 14693.\n";


      arma::mat Wmat = join_rows(Wmat_mu,DiagZ_Wmat_tau);

      double b=Wmat.n_cols;									// b is number of columns of W_bcf matrix (omega in the paper)

      //get t(orig_y_arma)inv(psi)J_bcf
      arma::mat ytW=orig_y_arma.t()*Wmat;								// orig_y_arma transpose W_bcf
      //get t(J_bcf)inv(psi)J_bcf
      arma::mat WtW=Wmat.t()*Wmat;							// W_bcf transpose W_bcf
      //get jpsij +aI
      arma::mat aI(b,b);									// create b by b matrix called aI. NOT INIIALIZED.
      aI=aI.eye();										// a times b by b identity matrix. The .eye() turns aI into an identity matrix.
      arma::vec a_vec_mu = a_mu*arma::ones<arma::vec>(b_mu);
      arma::vec a_vec_tau = a_tau*arma::ones<arma::vec>(b_tau);
      arma::vec a_vec(b);
      a_vec.head(b_mu) = a_vec_mu;
      a_vec.tail(b_tau) = a_vec_tau;
      aI.diag() = a_vec;

      arma::mat sec_term=WtW+aI;							//
      //arma::mat sec_term_inv=sec_term.i();					// matrix inverse expression in middle of eq 5 in the paper. The .i() obtains the matrix inverse.
      arma::mat sec_term_inv=inv_sympd(sec_term);					// matrix inverse expression in middle of eq 5 in the paper. The .i() obtains the matrix inverse.

      //get t(J_bcf)inv(psi)orig_y_arma
      arma::mat third_term=Wmat.t()*orig_y_arma;						// W_bcf transpose orig_y_arma
      //get m^TV^{-1}m
      arma::mat mvm= ytW*sec_term_inv*third_term;		// matrix expression in middle of equation 5
      //arma::mat rel=(b_mu*0.5)*log(a_mu)+(b_tau*0.5)*log(a_tau)-(1*0.5)*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);		// log of all of equation 5 (i.e. the log of the marginal likelihood of the sum of tree model)

      //Rcout << "Line 14724.\n";



      //arma::vec preds_temp_arma= Vmat*sec_term_inv*Wmat.t()*orig_y_arma;
      //arma::vec preds_temp_arma= Vmat*inv_sympd(sec_term)*Wmat.t()*orig_y_arma;
      //arma::vec preds_temp_arma= Vmat*inv_sympd(sec_term)*third_term;


      //double templik0=exp(arma::as_scalar((b_mu*0.5)*log(a_mu)+(b_tau*0.5)*log(a_tau)-(0.5)*log(det(sec_term))-expon*log(nu*lambdaBCF - mvm +yty)) );

      double templik0=exp(arma::as_scalar((b_mu*0.5)*log(a_mu)+(b_tau*0.5)*log(a_tau)-(0.5)*real(arma::log_det(sec_term))-expon*log(nu*lambdaBCF - mvm +yty)) );


      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      double templik = pow(templik0,beta_par);
      if(imp_sampler!=tree_prior){//check if importance sampler is not equal to the prior
        //templik=templik*(sum_tree_prior_prob/sum_tree_samp_prob);
        templik=templik*sum_prior_over_samp_prob;

      }
      overall_liks(j)= templik;



    if(is_test_data==1){
      //arma::mat zeromat(arma::size(Wmat_mu),arma::fill::zeros);
      //arma::mat zeromat(num_obs ,b_mu ,arma::fill::zeros);
      arma::mat zeromat=arma::zeros<arma::mat>(num_test_obs,b_mu);
      arma::mat Vmat = join_rows(zeromat,W_tilde_tau);

      arma::vec preds_temp_arma= Vmat*sec_term_inv*third_term;

      overall_preds(j)=preds_temp_arma*templik;


      //arma::mat covar_t=as_scalar((1/double(nu+num_obs))*(nu*lambda+yty-mvm))*(Vmat*sec_term_inv*(Vmat.t()));

      //arma::mat catevartemp=averagingvec.t()*covar_t*averagingvec;
      //arma::mat cattvartemp=catt_averagingvec.t()*covar_t*catt_averagingvec;
      //arma::mat catntvartemp=catnt_averagingvec.t()*covar_t*catnt_averagingvec;

      // preds_all_models_arma.col(i)=preds_temp_arma;
      // t_vars_arma.col(i)=covar_t.diag();
      // cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
      // cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i);
      // cate_vars_arma(i)=as_scalar(catevartemp);
      // catt_means_arma(i)=as_scalar(catt_averagingvec.t()*preds_temp_arma);
      // catt_means_weighted_arma(i)=catt_means_arma(i)*post_weights_arma(i);
      // catt_vars_arma(i)=as_scalar(cattvartemp);
      // catnt_means_arma(i)=as_scalar(catnt_averagingvec.t()*preds_temp_arma);
      // catnt_means_weighted_arma(i)=catnt_means_arma(i)*post_weights_arma(i);
      // catnt_vars_arma(i)=as_scalar(catntvartemp);
      //

    }else{

      //arma::mat zeromat(arma::size(Wmat_mu),arma::fill::zeros);
      //arma::mat zeromat(num_obs ,b_mu ,arma::fill::zeros);
      arma::mat zeromat=arma::zeros<arma::mat>(num_obs ,b_mu);
      arma::mat Vmat = join_rows(zeromat,Wmat_tau);

      arma::vec preds_temp_arma= Vmat*sec_term_inv*third_term;

      overall_preds(j)=preds_temp_arma*templik;


      //arma::mat covar_t=as_scalar((1/double(nu+num_obs))*(nu*lambda+yty-mvm))*(Vmat*sec_term_inv*(Vmat.t()));

      //arma::mat catevartemp=averagingvec.t()*covar_t*averagingvec;
      //arma::mat cattvartemp=catt_averagingvec.t()*covar_t*catt_averagingvec;
      //arma::mat catntvartemp=catnt_averagingvec.t()*covar_t*catnt_averagingvec;


      // preds_all_models_arma.col(i)=preds_temp_arma;
      // t_vars_arma.col(i)=covar_t.diag();
      // cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
      // cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i);
      // cate_vars_arma(i)=as_scalar(catevartemp);
      // catt_means_arma(i)=as_scalar(catt_averagingvec.t()*preds_temp_arma);
      // catt_means_weighted_arma(i)=catt_means_arma(i)*post_weights_arma(i);
      // catt_vars_arma(i)=as_scalar(cattvartemp);
      // catnt_means_arma(i)=as_scalar(catnt_averagingvec.t()*preds_temp_arma);
      // catnt_means_weighted_arma(i)=catnt_means_arma(i)*post_weights_arma(i);
      // catnt_vars_arma(i)=as_scalar(catntvartemp);
      //

    }//end of else statement (not test data)




  }//end of loop over all trees

}//end of pragma omp code


///////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
//Rcout << "Line 6852.\n";


//for(unsigned int i=0; i<overall_treetables.n_elem;i++){
//  pred_mat_overall = pred_mat_overall + overall_liks(i)*overall_treetables(i);
//}

if(is_test_data==1){
  #pragma omp parallel
  {
    arma::vec result_private=arma::zeros<arma::vec>(x_control_test_a.n_rows);
  #pragma omp for nowait //fill result_private in parallel
    for(unsigned int i=0; i<overall_preds.size(); i++) result_private += overall_preds(i);
  #pragma omp critical
    pred_vec_overall += result_private;
  }
}else{
  #pragma omp parallel
  {
    arma::vec result_private=arma::zeros<arma::vec>(x_control_a.n_rows);
  #pragma omp for nowait //fill result_private in parallel
    for(unsigned int i=0; i<overall_preds.size(); i++) result_private += overall_preds(i);
  #pragma omp critical
    pred_vec_overall += result_private;
  }
}


//Rcout << "Line 4030. \n";




//Rcout << "overall_liks = " << overall_liks << ". \n";
//Rcout << "max(overall_liks) = " << max(overall_liks) << ". \n";
//Rcout << "overall_liks[14] = " << overall_liks[14] << ". \n";


double sumlik_total= arma::sum(overall_liks);
//Rcout << "sumlik_total = " << sumlik_total << ". \n";

pred_vec_overall=pred_vec_overall*(1/sumlik_total);
//Rcout << "Line 1141 . \n";
//Rcout << "Line 1146 . \n";

//Rcout << "Line 4042. \n";
NumericVector orig_preds=get_original(min(ytrain),max(ytrain),-0.5,0.5,wrap(pred_vec_overall)) ;


return(orig_preds);

}
