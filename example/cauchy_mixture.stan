data {
  real gap;
}

parameters {
  real theta;
}

model{
 theta~ cauchy(-gap,0.2);   
 theta~ cauchy(gap,0.2);   
}

alternative model{
  theta~ normal(0,5);   
}

