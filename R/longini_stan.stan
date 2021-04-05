functions {

  real longini_pdf(real B, real Q, int j, int k);
  
  real longini_pdf(real B, real Q, int j, int k) {
    
    if (j == 0) {
      return B^k;
    } else {
      real t1;
      real t2;
      t1 = exp(lchoose(k,j))*(B^(k-j))*(Q^(j*(k-j)));
      t2 = 0;
        for (i in 0:(j-1)) {
          t2 = t2 + longini_pdf(B,Q,i,j);
        }

        return t1*(1-t2);
    }
  }

  real longini_zero_lpmf(int j, int k, real B, real Q) {
    real a;
    real b;

    a = 1.0- longini_pdf(B,Q,0,k);
    b = longini_pdf(B,Q,j,k);
    return log(b)-log(a);
  }


}

data {
  int N; //Number of households
  int k[N]; // Household sizes
  int j[N]; // Number of cases per household

  
}

parameters {
  real<lower=0,upper=1> B;
  real<lower=0,upper=1> Q;
}

model {

  vector[N] pdf_vals;

  for (i in 1:N) {
    pdf_vals[i] = longini_zero_lpmf(j[i] | k[i],B,Q );
  }

  target += (sum(pdf_vals)); 
}