functions{

int cumsum(int i){
    int sum = 0;
    if(i <= 0){
      sum = 0;
    }
    else{
      for(j in 1:i){
        sum += j;
      }
    }
    return sum;
}

real quantile_t(real u, real nu) {
  real y;
  real sign;

  if (u < 0.5) {
    y = 2 * u;
    sign = -1;
  } else {
    y = 2 * (1 - u);
    sign = 1;
  }

  if (u == 0.5) {
    return 0;
  } else {
    return sign * sqrt(nu / inv_inc_beta(nu / 2, 0.5, y) - nu);
  }
}

// https://discourse.mc-stan.org/t/linear-interpolation-and-searchsorted-in-stan/13318/6
real linear_interpolation_v(real x, vector x_pred, vector y_pred, int weight_opt){
    int K = rows(x_pred);
    vector[K] deltas = x - x_pred;
    vector[K] delta_sort;
    real ans;
    real t;
    real w;
    int i;
    real x1;
    real x2;
    real y1;
    real y2;

    if(x<x_pred[1] || x>x_pred[K]) reject("x is outside of the x_pred grid!");
    if(rows(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
    //this is which.min()
    // i = sort_indices_asc(abs(deltas))[3];
    delta_sort = sort_asc(abs(deltas));
    for(l in 1:K){
      if(abs(deltas[l]) == min(abs(deltas))){
        i = l;
      }
    }
    if(deltas[i]<=0) i -= 1;
    ans = y_pred[i];
    x1 = x_pred[i];
    x2 = x_pred[i + 1];
    y1 = y_pred[i];
    y2 = y_pred[i + 1];
    ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
    t = (x-x1)/(x2-x1);

    if(weight_opt == 1){
      w = 1-t;
    }else if(weight_opt == 2){
      w = 1/(1 + exp(1/(1-t) - 1/t));
    }else if(weight_opt == 3){
      w = 1 - 3*pow(t,2.0) + 2*pow(t,3.0);
    }

    ans = w*y1 + (1-w)*y2;

    return ans;
}

real frankF(real par, data vector x_grid, data vector y_grid){
    real franktau;

    franktau = linear_interpolation_v(par, x_grid, y_grid, 1);

    return franktau;
}

real joeF(real par){
    real param1;
    real tem;
    real tau;

    param1 = 2/par + 1;
    tem = digamma(2) - digamma(param1);
    tau = 1 + tem*2/(2-par);

    return tau;
}

vector system_frank(vector y, vector tau, data vector x_grid, data vector y_grid){
    vector[1] ftau = tau;
    vector[1] par;

    par[1] = ftau[1] - frankF(y[1], x_grid, y_grid);

    return par;
}

vector system_joe(vector y, vector tau){
    vector[1] par;

    par[1] = tau[1] - joeF(y[1]);

    return par;
}

vector Tau2Par(vector tau, data array[,] int fam, vector y_guess, data vector x_grid, data vector y_grid, data int d){
  int i;
  int k;
  int itau;
  vector[1] vtau;
  vector[(d*(d-1))%/%2] par;

  for(ifor in 1:(d-1)){
    i = d - ifor;
    for(kfor in (i+1):d){
      k = d + i + 1 - kfor;
      itau = (d - k) * (d - 1) - cumsum(d - k) + d - i;
      if(fam[k][i] == 0){ // independence
        par[itau] = 0;
      }else if((fam[k][i] ==  1) || (fam[k][i] ==  2)){ // Gauss and t
        par[itau] = sin(tau[itau]*pi()/2);
      }else if((fam[k][i] ==  3) || (fam[k][i] == 13)){ // Clayton 0° and 180°
        par[itau] = 2*tau[itau]/(1-tau[itau]);
      }else if((fam[k][i] == 23) || (fam[k][i] == 33)){ // Clayton 90° and 270°
        par[itau] = 2*tau[itau]/(1+tau[itau]);
      }else if((fam[k][i] ==  4) || (fam[k][i] == 14)){ // Gumbel 0° and 180°
        par[itau] = 1/(1-tau[itau]);
      }else if((fam[k][i] == 24) || (fam[k][i] == 34)){ // Gumbel 90° and 270°
        par[itau] = -1/(1+tau[itau]);
      }else if(fam[k][i]  ==  5){ // Frank
        vtau[1] = tau[itau];
        par[itau] = solve_newton(system_frank, y_guess, vtau, x_grid, y_grid)[1];
      }else if((fam[k][i] ==  6) || (fam[k][i] == 16)){ // Joe 0° and 180°
        vtau[1] = tau[itau];
        par[itau] = solve_newton(system_joe, y_guess, vtau)[1];
      }else if((fam[k][i] == 26) || (fam[k][i] == 36)){ // Joe 90° and 270°
        vtau[1] = -tau[itau];
        par[itau] = -solve_newton(system_joe, y_guess, vtau)[1];
      }
    }
  }

  return par;
}

//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// loglik    output
//////////////////////////////////////////////////////////////
real LL(int family, real u_old, real v_old, real par) {
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  real XINFMAX = 1.79769e+308;
  real DBLMIN=1e-4;//?????
  real ll; real t1=0.0; real t2=0.0; real f;vector[2] temp;

  real u = u_old;
  real v = v_old;


  if(u < UMIN) {
    u = UMIN;
  } else if (u > UMAX) {
    u = UMAX;
  }
  if (v < UMIN) {
    v = UMIN;
  } else if (v > UMAX) {
    v = UMAX;
  }

  //Compute log-likelihood:
  if (family == 0) {//independent
    ll = 0;
  } else if (family == 1) {//Gaussian
    t1 = inv_Phi(u);
    t2 = inv_Phi(v);
    f = 1.0 / sqrt(1.0 - pow(par, 2.0)) * exp((pow(t1, 2.0) + pow(t2, 2.0))
        / 2.0 + (2.0 * par * t1 * t2 - pow(t1, 2.0) - pow(t2, 2.0))
        / (2.0 * (1.0 - pow(par, 2.0))));
    if (log(f) > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < DBLMIN) {
      ll = log(DBLMIN);
    } else {
      ll = log(f);
    }
  } else if (family == 3) { //Clayton
    if(par == 0) {
      ll = 0;
    } else if (par < 1e-10) {
      ll = 0;
    } else {
      f = log1p(par) - (1.0 + par) * log(u * v) - (2.0 + 1.0 / par) * log(pow(u, -par) + pow(v, -par) - 1.0);
      if (f > XINFMAX) {
        ll = log(XINFMAX);
      } else if (f < log(DBLMIN)) {
        ll = log(DBLMIN);
      } else {
        ll = f;
      }
    }
  } else if (family == 4) { //Gumbel
    t1 = pow(-log(u), par) + pow(-log(v), par);
    f  = -pow(t1, 1.0 / par) + (2.0 / par - 2.0) * log(t1) + (par - 1.0)
        * log(log(u) * log(v)) - log(u * v) + log1p((par - 1.0)
        * pow(t1, -1.0 / par));
    if (f > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < log(DBLMIN)) {
      ll = log(DBLMIN);
    } else {
      ll = f;
    }
  } else if (family == 5) {// Frank
    if (abs(par) < 1e-10) {
      ll = 0;
    } else {
      f = (par * (exp(par) - 1.0) * exp(par * v + par * u + par))
          / pow(exp(par * v + par * u) - exp(par * v + par)
          - exp(par * u + par) + exp(par), 2.0);
      if( log(f) > XINFMAX) {
        ll = log(XINFMAX);
      } else if (f < DBLMIN) {
        ll = log(DBLMIN);
      } else {
        ll = log(f);
      }
    }
  } else if (family == 6)	{//Joe
    f = pow(pow(1 - u, par) + pow(1 - v, par) - pow(1 - u, par)
        * pow(1 - v, par), 1 / par - 2) * pow(1 - u, par - 1)
        * pow(1 - v, par - 1) * (par - 1 + pow(1 - u, par) + pow(1 - v, par)
        - pow(1 - u, par) * pow(1 - v, par));
    if (log(f) > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < DBLMIN) {
      ll = log(DBLMIN);
    } else {
      ll = log(f);
    }
  } else if (family == 13) {//rotated Clayton (180?)
    if (par == 0) {
      ll = 0;
    } else if (par < XEPS) {
      ll = 0;
    } else {
      f = (1.0 + par) * pow((1 - u) * (1 - v), -1.0 - par)
          * pow(pow((1 - u), -par) + pow((1 - v), -par) - 1.0, -2.0 - 1.0 / par);
      temp[1]=f;
      temp[2]=0;
      f = max(temp);
      if (log(f) > XINFMAX) {
        ll = log(XINFMAX);
      } else if (f < DBLMIN) {
        ll = log(DBLMIN);
      } else {
        ll = log(f);
      }
    }
  } else if (family == 14) {//rotated Gumbel (180?)
    t1 = pow(-log(1 - u), par) + pow(-log(1 - v), par);
    t2 = exp(-pow(t1, 1.0 / par));
    f = t2 / ((1 - u) * (1 - v)) * pow(t1, -2.0 + 2.0 / par)
        * pow(log(1 - u) * log(1 - v), par - 1.0) * (1.0 + (par - 1.0)
        * pow(t1, -1.0 / par));
    if (log(f) > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < DBLMIN) {
      ll = log(DBLMIN);
    } else {
      ll = log(f);
    }
  } else if (family == 16) {//rotated Joe (180?)
    f = pow(pow(u, par) + pow(v, par) - pow(u, par) * pow(v, par), 1 / par - 2)
        * pow(u, par - 1) * pow(v, par - 1) * (par - 1 + pow(u, par)
        + pow(v, par) - pow(u, par) * pow(v, par));
    if (log(f) > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < DBLMIN) {
      ll = log(DBLMIN);
    } else {
      ll = log(f);
    }
  }

  //Write to output vector:
  return ll;
}

real LL(int family, real u_old, real v_old, real par, real par2) {
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  real XINFMAX = 1.79769e+308;
  real DBLMIN=1e-4;//?????
  real ll; real t1=0.0; real t2=0.0; real f;vector[2] temp;

  real u = u_old;
  real v = v_old;


  if(u < UMIN) {
    u = UMIN;
  } else if (u > UMAX) {
    u = UMAX;
  }
  if (v < UMIN) {
    v = UMIN;
  } else if (v > UMAX) {
    v = UMAX;
  }

  //Compute log-likelihood:
  if (family == 2) { // t
    t1 = quantile_t(u, par2);
    t2 = quantile_t(v, par2);
    f = 1 / (2 * pi() * sqrt(1.0 - pow(par, 2.0))
        * exp(student_t_lpdf(t1 | par2, 0, 1))
        * exp(student_t_lpdf(t2 | par2, 0, 1))) * pow(1.0 + (pow(t1, 2.0)
        + pow(t2, 2.0) - 2.0 * par * t1 * t2)
        / (par2 * (1.0 - pow(par, 2.0))), -(par2 + 2.0) / 2.0);
    if (log(f) > XINFMAX) {
      ll = log(XINFMAX);
    } else if (f < DBLMIN) {
      ll = log(DBLMIN);
    } else {
      ll = log(f);
    }
  }

  //Write to output vector:
  return ll;
}



real LL_mod2(int family, real u_old, real v_old, real theta, real par2) {
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  int nfamily;

  real u;
  real v;
  real loglik;

  u = u_old;
  v = v_old;

  if (u < UMIN){
    u = UMIN;
  } else if (u > UMAX) {
    u = UMAX;
  }
  if (v < UMIN) {
    v = UMIN;
  } else if (v > UMAX) {
    v = UMAX;
  }
  if (family == 23 || family == 24 || family == 26) {// 90? rotated copulas
    nfamily = family - 20;
    loglik = LL(nfamily, 1 - u,  v, -theta);
  } else if (family == 33 || family == 34 || family == 36) {// 270? rotated copulas
    nfamily = family - 30;
    loglik = LL(nfamily, u,  1 - v, -theta);
  } else if (family == 2) {
    loglik = LL(family, u,  v, theta, par2);
  }else {
    loglik = LL(family, u,  v, theta);
  }

  return loglik;
}

//////////////////////////////////////////////////////////////
// Function to compute h-function for vine simulation and estimation
// Input:
// family   copula family (0=independent,  1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
// out      output
//////////////////////////////////////////////////////////////
real Hfunc(int family, real u_old, real v_old, real theta) {
  real h;
  vector[2] temp2;
  vector[2] temp1;
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  real x;
  real out;
  int nu = 14;

  real u = u_old;
  real v = v_old;

  if (v == 0 || u == 0) {
    h = 0;
  } else if (v == 1) {
    h = u;
  } else {
    if(family==0) {//independent
      h = u;
    } else if (family == 1) {//gaussian
      x = (inv_Phi(u) - theta * inv_Phi(v)) / sqrt(1.0 - pow(theta, 2.0));
      if (x < positive_infinity()){
        h = Phi(x);
      } else if (inv_Phi(u) - theta * inv_Phi(v) < 0) {
        h = 0;
      } else {
        h = 1;
      }
    } else if(family == 3) {//clayton
      if(theta == 0) {
        h = u;
      }
      if(theta < XEPS) {
        h = u;
      } else {
        x = pow(u, -theta) + pow(v, -theta) - 1.0 ;
        h = pow(v, -theta - 1.0) * pow(x, -1.0 - 1.0 / theta);
        if(theta < 0){
          if(x < 0) {
            h = 0;
          }
        }
      }
    } else if (family == 4) {//gumbel
      if(theta == 1) {
        h = u;
      } else {
        h = -(exp(-pow(pow(-log(v), theta) + pow(-log(u), theta), 1.0/theta))
            * pow(pow(-log(v), theta) + pow(-log(u), theta), 1.0 / theta - 1.0)
            * pow(-log(v), theta)) / (v * log(v));
      }
    } else if (family == 5) {//frank
      if (theta == 0) {
        h = u;
      } else {
        h = -(exp(theta) * (exp(theta * u) - 1.0)) / (exp(theta * v + theta * u)
            - exp(theta * v + theta) - exp(theta * u + theta) + exp(theta));
      }
    } else if (family == 6) {//joe
      if(theta == 1) {
        h = u;
      } else {
        h = pow(pow(1.0 - u, theta) + pow(1.0 - v, theta) - pow(1.0 - u, theta)
        * pow(1.0 - v, theta), 1.0 / theta - 1) * pow(1.0 - v, theta - 1.0)
        * (1 - pow(1 - u, theta));
      }
    } else if (family == 13) {//rotated clayton (180?)
      if(theta == 0) {
        h = u;
      }
      if(theta < XEPS) {
        h = u;
      } else {
        u = 1 - u;
        v = 1 - v;
        x = pow(u, -theta) + pow(v, -theta) - 1.0;
        h = pow(v, -theta - 1.0) * pow(x, -1.0 - 1.0 / theta);
        h = 1 - h;
        u = 1 - u;
        v = 1 - v;
      }
    } else if (family == 14) {//rotated gumbel (180?)
      v = 1 - v;
      u = 1 - u;
      h = -(exp(-pow(pow(-log(v), theta) + pow(-log(u), theta), 1.0 / theta))
          * pow(pow(-log(v), theta) + pow(-log(u), theta), 1.0 / theta - 1.0)
          * pow(-log(v), theta)) / (v *	log(v));
      h = 1 - h;
      u = 1 - u;
      v = 1 - v;
    } else if (family == 16) {
      v = 1 - v;
      u = 1 - u;
      h = pow(pow(1.0 - u, theta) + pow(1.0 - v, theta) - pow(1.0 - u, theta)
          * pow(1.0 - v, theta), 1.0 / theta - 1) * pow(1.0 - v, theta - 1.0)
          * (1 - pow(1 - u, theta));
      h = 1 - h;
      u = 1 - u;
      v = 1 - v;
    }
  }
  temp1[1] = h;
  temp1[2] = UMAX;
  temp2[1] = min(temp1);
  temp2[2] = UMIN;
  out = max(temp2);

  return out;
}

real Hfunc(int family, real u_old, real v_old, real theta, real par2) {
  real h;
  vector[2] temp2;
  vector[2] temp1;
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  real x;
  real out;
  // int nu = 14;

  real u = u_old;
  real v = v_old;

  if (v == 0 || u == 0) {
    h = 0;
  } else if (v == 1) {
    h = u;
  } else {
    if (family == 2) {//student
      real t1=0.0; real t2=0.0; real mu; real sigma2;
      t1 = quantile_t(u, par2);
      t2 = quantile_t(v, par2);
      mu = theta*t2;
      sigma2 = ((par2+t2*t2)*(1.0-theta*(theta)))/(par2+1.0);
      h = student_t_cdf((t1-mu)/sqrt(sigma2)|par2+1.0,0,1);
    }
  }
  temp1[1] = h;
  temp1[2] = UMAX;
  temp2[1] = min(temp1);
  temp2[2] = UMIN;
  out = max(temp2);

  return out;
}


// Since the h function is not symmetric in case of real Gumbel and real Clayton we have two implement both separately,
// i.e. Hfunc1 and Hfunc2
real Hfunc1(int family, real u_old, real v_old, real theta, real par2) {
  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;
  vector[2] temp1;
  vector[2] temp2;
  real u = u_old;
  real v = v_old;
  real out;
  real negv;
  real negu;
  int nfamily;
  int T;

  if(u < UMIN) {
    u = UMIN;
  } else if (u > UMAX) {
    u = UMAX;
  }
  if(v < UMIN) {
    v = UMIN;
  } else if (v > UMAX) {
    v = UMAX;
  }
  if (family == 23 || family == 24 || family == 26) {
    nfamily = family - 20;
    negv = 1 - v;
    out = Hfunc(nfamily, u, negv, -theta);
  } else if (family == 33 || family == 34 || family == 36) {
    nfamily = family - 30;
    negu = 1 - u;
    out = Hfunc(nfamily, negu, v, -theta);
    out = 1 - out;
  } else if (family == 2) {
    out = Hfunc(family, u, v, theta, par2);
  } else {
    out = Hfunc(family, u, v, theta);
  }
  // ensure that results are in [0,1]
  temp1[1] = 0;//UMIN????
  temp1[2] = out;
  temp2[1] = 1;//UMAX????
  temp2[2] = max(temp1);
  out = min(temp2);

  return out;
}

real Hfunc2(int family, real v_old, real u_old, real theta, real par2) {
  real negv;
  real negu;

  real UMAX = 1-1e-10;
  real UMIN = 1e-10;
  real XEPS =1e-4;

  vector[2] temp1;
  vector[2] temp2;
  real u = u_old;
  real v = v_old;
  real out;
  int nfamily;

  if(u < UMIN) {
    u = UMIN;
  } else if (u > UMAX) {
    u = UMAX;
  }
  if (v < UMIN) {
    v = UMIN;
  } else if (v > UMAX) {
    v = UMAX;
  }
  if (family == 23 || family == 24 || family == 26) {
    nfamily = family - 20;
    negv = 1 - v;
    out = Hfunc(nfamily, negv, u, -theta);
    out = 1 - out;
  } else if (family == 33 || family == 34 || family == 36) {
    nfamily = family - 30;
    negu = 1 - u;
    out = Hfunc(nfamily, v, negu, -theta);
  } else if (family == 2) {
    out = Hfunc(family, v, u, theta, par2);
  } else {
    // switch u and v
    out = Hfunc(family, v, u, theta);
  }

  // ensure that results are in [0,1]
  temp1[1] = UMIN;
  temp1[2] = out;
  temp2[1] = UMAX;
  temp2[2] = max(temp1);
  out = min(temp2);

  return out;
}




real VineLogLikRvine2(int T, int d, array[,] int fam, array[] int maxmat,
                      array[] int matri, array[] int condirect,
                      array[] int conindirect, vector par, data array[,] real par2,
                      matrix udata) {
  int l;
  int m;
  vector[T] out;
  array[d,d] real value2;
  array[d,T] real x;
  array[d,d] real vdirect;
  array[d,d] real vindirect;
  int i;
  int k;
  vector[T] sumsplitlog;

  for (t in 1:T) {
    sumsplitlog[t] = 0;
  }

  //Initialize
  // Not needed, just use udata
  for (ifor in 1:d) {
    for (t in 1:T) {
      x[ifor][t] = udata[t,ifor];
    }
  }

  for (t in 1:T) {
    for (ifor in 1:d) {
      vdirect[d][ifor]=x[d-ifor+1][t];
    }
    for(ifor in 1:(d-1)) {
      i = d - ifor;
      for(kfor in (i+1):d) {
        k = d + i + 1 - kfor;
        m = maxmat[k + d * (i - 1)];
        if (m == matri[k + d * (i - 1)]) {
          value2[k][i] = LL_mod2(fam[k][i], vdirect[k][d-m+1], vdirect[k][i],
                           par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                           par2[k][i]);
          if(condirect[k + d * (i - 1) -1] == 1) {
            vdirect[k - 1][i] = Hfunc1(fam[k][i], vdirect[k][i],
                                       vdirect[k][d - m + 1],
                                       par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                                       par2[k][i]);
          }
          if(conindirect[k + d * (i - 1) - 1] == 1) {
            vindirect[k - 1][i] = Hfunc2(fam[k][i], vdirect[k][d - m + 1],
                                         vdirect[k][i],
                                         par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                                         par2[k][i]);
          }
        } else {
          value2[k][i] = LL_mod2(fam[k][i], vindirect[k][d - m + 1],
                           vdirect[k][i],
                           par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                           par2[k][i]);
          if(condirect[k + d * (i - 1) -1] == 1) {
            vdirect[k - 1][i] = Hfunc1(fam[k][i], vdirect[k][i],
                                       vindirect[k][d - m + 1],
                                       par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                                       par2[k][i]);
          }
          if(conindirect[k + d * (i - 1) -1] == 1) {
            vindirect[k - 1][i] = Hfunc2(fam[k][i], vindirect[k][d - m + 1],
                                         vdirect[k][i],
                                         par[(d - k) * (d - 1) - cumsum(d - k) + d - i],
                                         par2[k][i]);
          }
        }
        sumsplitlog[t] += value2[k][i];
      }
    }
    out[t] =  sumsplitlog[t];
  }

  return sum(out);
}

}

data{
  int n; // number of observations
  int d; // dimension
  vector[(d * (d - 1)) %/% 2] L_bound;   // vector of lower bounds for Kendall's tau
  vector[(d * (d - 1)) %/% 2] U_bound;   // vector of upper bounds for Kendall's tau
  array[d*d] real para2; // array for second copula parameters
  array[d*d] int family; // array for copula families
  array[d*d] int maxmat; // array for likelihood calculations
  array[d*d] int matri; // array for likelihood calculations
  array[d*d] int condirect; // array for likelihood calculations
  array[d*d] int conindirect; // array for likelihood calculations
  matrix[n,d] udata; // copula data

  vector[108] x_grid; // grid of Frank copula parameters
  vector[108] y_grid; // grid of Frank copula Kendall's tau values
}

transformed data{
  vector[1] y_guess = [0.5]'; // Guess for algebra solver
  array[d,d] int fam; // adjusted array for copula families
  array[d,d] real par2; // adjusted array for second copula parameters

  for(ifor in 1:d){
    for(j in 1:d){
      fam[ifor][j] = family[ifor+(d)*(j-1)];
      par2[ifor][j] = para2[ifor+(d)*(j-1)];
    }
  }
}

parameters{
    vector<lower=L_bound, upper=U_bound>[(d * (d - 1)) %/% 2] tau; // Kendall's tau
}

transformed parameters{
  vector[(d * (d - 1)) %/% 2] par; // copula parameters

  // transform tau to copula parameter according to family
  par = Tau2Par(tau, fam, y_guess, x_grid, y_grid, d);
}

model{
  // Calculate Loglikelihood
  target += VineLogLikRvine2(n, d, fam, maxmat, matri, condirect, conindirect, par, par2, udata);
}
