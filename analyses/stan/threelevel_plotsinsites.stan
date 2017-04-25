// Aloha!
// By Megan & Lizzie started on 6 January 2017
// Trying to do a simple 3-level model
// Imagine some data: 
    // You have observed densities of something across plots, and plots occur within sites
    // 3 level: 
        // (1) Observation level (density)
        // (2) plot level (multiple observations taken within each plot)
        // (3) Multiple plots for each site
    // Predictor data (aka x) observed only at observation level (imagine you measured moisture around each plant or such)
    // ... (that is, no additional predictors added at plot or site level)
// Model structure does not assume same number of plots across sites
// Code below includes pooled intercepts and slopes at all levels and ...
// assumes that all sites have the same within-site (plot-to-plot) variance.  
// A link (we did not use it, but I am keeping it here for now): http://stackoverflow.com/questions/29379001/nested-model-in-stan

data {
  int N;    //total number of observations (e.g., N=120=3 observations in each plot x 4 plots per site x 10 sites )
  int J;    //total number of plots 
  int S;    //total number of sites
  vector[N] y;    // the data, y, is a vector of length N
  int plotnum[N];  // column of plot number identifiers for each obs
  int sitenum[J];  // column of site number identifiers for each plot
  vector[N] x;        // vector of predictor X for each obs
}

parameters {
  vector[J] a_plot;    //estimated intercept for each plot
  vector[J] b_plot;    //estimated slope for each plot
  vector[S] a_site;    //estimated intercept for each site
  vector[S] b_site;    //estimated slope for each site
  
  real<lower=0> sig_a_site;  //variance in intercept across plots; 
#  real<lower=0> sig_a_site[S];  //variance in intercept across plots; 
      // the intercept for plot j in site s is drawn from a distribution with mean a_site[s]
      // and standard deviation sig_a_site[s]
  real<lower=0> sig_b_site;  //variance in slopes across plots; 
#  real<lower=0> sig_b_site[S];  //variance in slopes across plots; 
      // the slope for plot j in site s is drawn from a distribution with mean b_site[s] 
      /// and standard deviation sig_b_site[s]
  real mu_a;                    //mean intercept across sites; 
      // the site intercept for site s is drawn from distribution with mean mu_a...
  real<lower=0> sig_a;          //...and standard deviation sig_a
  real mu_b;                    //mean slope across sites; 
      // the site slope are drawn from distribution with mean mu_b...
  real<lower=0> sig_b;          //...and standard deviation sig_b
  real<lower=0> sig_y;          // observation error
}

model {
  real yhat[N];                 //vector of predicted y's (`ypred' for Lizzie)
  
  for (i in 1:N) {
    yhat[i] = a_plot[plotnum[i]] + b_plot[plotnum[i]]*x[i];
  }

  // //For estimating separate within-site variacnes
  // for (j in 1:J){
  //   a_plot[j] ~ normal(a_site[sitenum[j]], sig_a_site[sitenum[j]]);
  //   b_plot[j] ~ normal(b_site[sitenum[j]], sig_b_site[sitenum[j]]);
  // }

  //For estimating a single value for all within-site variacnes
  for (j in 1:J){
    a_plot[j] ~ normal(a_site[sitenum[j]], sig_a_site);
    b_plot[j] ~ normal(b_site[sitenum[j]], sig_b_site);
  }
    
  a_site ~ normal(mu_a,sig_a);
  b_site ~ normal(mu_b,sig_b);
  
  y ~ normal(yhat,sig_y);        //data is distributed normally around predicted (yhat) with s.d. sig_y (this is error of data around predicted values)
  

}
