functions {

  real longini_pdf(int j, int k , real B, real Q );
  
  real longini_pdf(int j, int k, real B, real Q) {
    
    if (j == 0) {
      return B^k;
    } else {
      real t1;
      real t2;
      t1 = exp(lchoose(k,j))*(B^(k-j))*(Q^(j*(k-j)));
      t2 = 0;
        for (i in 0:(j-1)) {
          t2 = t2 + longini_pdf(i,j,B,Q);
        }

        return t1*(1-t2);
    }
  }

  real longini_zero_lpmf(int j[ ], int k[ ], real B, real Q) {
    real a;
    real b;

    a = 1.0- longini_pdf(0,k,B,Q);
    b = longini_pdf(j,k,B,Q);
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

    target += longini_zero_lpmf(j | k,B,Q );
  

}