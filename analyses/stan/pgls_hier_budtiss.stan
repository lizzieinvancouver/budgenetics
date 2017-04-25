// Formerly called pgls_lemoine, updated for the budtiss data, now trying to do it right avoid duplicating tips //
// Edits by Lizzie so far //
// Have not tried running this yet ... // 

data{
	int N; // number of spp 
	int K; // 1 
        int J; // total number of species 
	vector[N] y; // doy
	matrix[N, N] V; // will be VCV
	matrix[N, N] Lmat; // a pre-designed smatrix
	vector[N] X1; // warming
	vector[N] X2; // photo
}
transformed data{
	real Ndiv; // a real number (holder)
	matrix[N, N] Ident; // a N by N matrix (holder)

	Ndiv = N; // fill out Ndiv
	Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
}
parameters{
	real<lower=0> sigma;
	real<lower=0, upper=1> lambda;
	vector[J] a; // intercepts
        vector[J] B1; // warming effect
        vector[J] B2; // photo effect
}
transformed parameters{
}
model{
	matrix[N,N] Vlambda; // VCV x lambda I think ... 
	matrix[N,N] Vsigma; // (VCVxlambda) x sigma I think ... 
	vector[N] yhat; // holder for predicted values
	real detV; // determinant of the VCV ... I am not sure if we need this since we do multinormal now
	real cdf;  // cumulative distribution f(x) ... I am not sure if we need this since we do multinormal now
	real logLike_PGLS; // a log likelihood! ... I am not sure if we need this since we do multinormal now

	Vlambda = (lambda*Lmat + Ident) .* V;
	Vsigma = sigma^2*Vlambda;

 for (i in 1:N) {
	yhat[i] = a[i] + B2[i]*X1 + B2[i]*X2; 
   }    
	
        a ~ normal(mu_a, sig_a); // should add priors to these .... 
	y ~ multi_normal(yhat, Vsigma); 
	
	B ~ normal(0, 1);
	sigma ~ cauchy(0, 5); // upping this from last time, sigma may be big.... 
	lambda ~ beta(1, 1);
}
