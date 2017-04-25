// Formerly called pgls_lemoine, updated for the budtiss data //
// Edits by Lizzie so far //

data{
	int N; // number of spp 
	int K; // 1 (WTF?)
	vector[N] y; // doy
	matrix[N, N] V; // will be VCV
	matrix[N, N] Lmat; // some sort of pre-designed matrix
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
	vector[4] B;
}
transformed parameters{
}
model{
	matrix[N,N] Vlambda; // VCV x lambda??
	matrix[N,N] Vsigma; // (VCVxlambda) x sigma??
	vector[N] yhat; // holder for predicted values
	real detV; // determinant of the VCV
	real cdf;  // cumulative distribution f(x)??
	real logLike_PGLS; // a log likelihood!

	Vlambda = (lambda*Lmat + Ident) .* V;
	Vsigma = sigma^2*Vlambda;

	yhat = B[1] + B[2]*X1 + B[3]*X2;


	//detV <- log_determinant(Vsigma);
	//cdf <- ((y-yhat)'*inverse(Vsigma)*(y-yhat));
	//logLike_PGLS <-  -0.5*(detV + cdf);
	//increment_log_prob(logLike_PGLS);
	
	y ~ multi_normal(yhat, Vsigma); // Better than the above!
	
	B ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5);
	lambda ~ beta(1, 1);
}
