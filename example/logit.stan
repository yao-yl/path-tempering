data {
	int n;
	int y[n];
	real x[n];
}

parameters {
	real beta;
}
model {
  for (i in 1:n)
    y[i]~ bernoulli_logit(beta * x[i]);
}
alternative model {
 for (i in 1:n)
    y[i]~  bernoulli(Phi(beta * x[i]));
}
