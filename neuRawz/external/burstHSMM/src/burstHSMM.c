#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
double maxvec(double *x, int n);
int rdisc(int n,double *prob,int islog);
double norm2(double *vec, int len);
void Rprintvec(char *a, double *x, int len);

int n, L, Lm, len[2], **chpt, tot_nisi[2];
double tot_isi[2], tot_lisi[2], MaxBurstISI;
double *y, ***g, ***el, **gnorm, gtot, gmax[2], acpt_ctr[4];
double isi_var[2], isi_shape[2], isi_mean[2], isi_scale[2], isi_rate[2];
double dur_var[2], dur_shape[2], dur_mean[2], dur_scale[2], dur_rate[2];
double nconst[2], lnconst[2];

double prob_move(double run, double isi, double a, double b, int to_switch, int log_p);
void block_update(int *s, int i0, int m, double *r, double *par);
void par_update2(int *s, double *r, double *par, double *hpar, double *sc, double *logpost, double *acpt, int M);
double lpost(int *s, double *r, double *par, double *hpar, int get_info);
double lpost_incr(int *s, double *r, double *par, double *hpar, int par_set, double *incr);
double disi(double x, int state, int getlog);

#define min(a,b) a < b ? a : b

//++++++ random initialization of states ++++++

void rstate(int *len, double *times, double *par, int *state, double *limit){

	int i, k, n = len[0];
	double lp[2], r;
	
	MaxBurstISI = limit[0];
	
	for(k = 0; k < 2; k++){
		
		dur_mean[k]  = exp(par[5 + 2 * k]);
		dur_shape[k] = exp(par[4 + 2 * k]);
		dur_scale[k] = exp(par[5 + 2 * k] - par[4 + 2 * k]);
		dur_rate[k] = exp(par[4 + 2 * k] - par[5 + 2 * k]);
		
	}
	
	state[0] = 1;
	r = times[0];
	GetRNGstate();	

	for(i = 1; i < n; i++){
		if(times[i] > MaxBurstISI)
			state[i] = 1;
		else {
			for(k = 0; k < 2; k++)
				lp[k] = prob_move(r, times[i - 1], dur_shape[state[i - 1]], dur_scale[state[i - 1]], k ^ state[i - 1], 1);
			
			state[i] = rdisc(2, lp, 1);
		}
		
		r = (1 - (state[i] ^ state[i - 1])) * r + times[i];
	}	

	PutRNGstate();	
}


//++++++++++++++++++ Main MCMC for burstHSMM ++++++++++++++++++++++++

void hsmm(double *isi, double *par, double *hpar, int *state, int *dims, int *mcmc_len, 
		  double *run, int *st_state, double *st_par, double *acpt, double *tune, double *limit){

	int i, j, k, i1, m, M;
	int w_min = dims[1], w_max = dims[2], w_gap = w_max - w_min;
	int nburn = mcmc_len[0], nsamp = mcmc_len[1], nskip = mcmc_len[2], nsweep = nburn + nsamp * nskip;
	
	n = dims[0];
	y = isi;
	L = mcmc_len[3];
	Lm = L - 1;
	MaxBurstISI = limit[0];
//	chpt = imat(2, n);
	chpt = (int **)R_alloc(2, sizeof(int *));
	for(i = 0; i < 2; i++)
		chpt[i] = (int *)R_alloc(n, sizeof(int));
	
	g = (double ***)R_alloc(w_max, sizeof(double **));
	for(i = 0; i < w_max; i++){
		g[i] = (double **)R_alloc(2, sizeof(double *));
		for(k = 0; k < 2; k++)
			g[i][k] = (double *)R_alloc(i + 1, sizeof(double));
	}
	gnorm = (double **)R_alloc(w_max, sizeof(double *));
	for(i = 0; i < w_max; i++)
		gnorm[i] = (double *)R_alloc(2, sizeof(double));
	
	el = (double ***)R_alloc(w_max, sizeof(double **));
	for(i = 0; i < w_max; i++){
		el[i] = (double **)R_alloc(2, sizeof(double *));
		for(k = 0; k < 2; k++)
			el[i][k] = (double *)R_alloc(i + 1, sizeof(double));
	}
	
	run[0] = y[0]; j = 0; 
	for(i = 1; i < n; i++){
		if(state[i] ^ state[i - 1])
			run[i] = y[i];
		else
			run[i] =  run[i - 1] + y[i];
	}
	
	
	double logpost[] = {lpost(state, run, par, hpar, 1)};
	Rprintf("Initial log posterior = %g\n", logpost[0]);
	
	for(k = 0; k < 4; k++){
		acpt[k] = 0.0;
		acpt_ctr[k] = 0.0;
	}
	
	GetRNGstate();
	
	int flag = 1, sweep_count = 1, thik_achey = 1;
	for(i = 0; i < nburn; i++){
		m = w_min + ceil(w_gap * runif(0.0, 1.0));
		i1 = n;
		while(i1 > m){
			block_update(state, i1 - m, m, run, par);
			i1 -= m;
			while(i1 > 0 & state[i1 - 1] == state[i1])
				i1--;
		}
		if(i1 > 0)
			block_update(state, 0, i1, run, par);
		
		logpost[0] = lpost(state, run, par, hpar, 1);
		
		par_update2(state, run, par, hpar, tune, logpost, acpt, L);
		if(flag == 0){
			
			Rprintf("\nSweep = %d, log posterior = %g\n", sweep_count, logpost[0]);
			Rprintvec("Parameters = ", par, 8);
			Rprintf("\n");
		}
		
		flag = (flag + 1)%(nsweep / 10);
		sweep_count++;
	}
	
	
	Rprintf("Done burning\n");
	
	int skip = 0, skip2 = 0, i2;
	
	for(i = 0; i < nsamp; i++){
		for(j = 0; j < nskip; j++){
			
			m = w_min + ceil(w_gap * runif(0.0, 1.0));
			i1 = n;
			while(i1 > m){
				block_update(state, i1 - m, m, run, par);
				i1 -= m;
				while(i1 > 0 & state[i1 - 1] == state[i1])
					i1--;
			}
			if(i1 > 0)
				block_update(state, 0, i1, run, par);
			
			logpost[0] = lpost(state, run, par, hpar, 1);
			
			thik_achey = 1;
			i2 = 0;
			while(thik_achey & (i2 < n)){
				i2++;
				if(!state[i2])
					thik_achey = (y[i2] < MaxBurstISI);
			}
			
			if(i2 < n)
				Rprintf("Thik ney!! state[%d] = %d, isi[%d] = %g, bound = %g\n",
						i2, state[i2], i2, y[i2], MaxBurstISI);
			
			
			par_update2(state, run, par, hpar, tune, logpost, acpt, L);
			if(flag == 0){
				
				Rprintf("\nSweep = %d, log posterior = %g\n", sweep_count, logpost[0]);
				Rprintvec("Parameters = ", par, 8);
				Rprintf("\n");
			}
			
			flag = (flag + 1)%(nsweep / 10);
			sweep_count++;
			
		}
		
		for(j = 0; j < n; j++)
			st_state[skip++] = state[j];
		for(k = 0; k < 8; k++)
			st_par[skip2++] = par[k];
		
	}
	
	
	for(k = 0; k < 4; k++)
		acpt[k] /= acpt_ctr[k];
	Rprintvec("Done looping with acceptance = ", acpt, 4);
	PutRNGstate();
}

// ++++++++ Computes probability of the state run ++++++++++++++++
// ++++++++++++ continuiung beyond current isi +++++++++++++++++++

double prob_move(double run, double isi, double a, double b, int to_switch, int log_p){

	double v, lv;
	lv = pgamma(run, a, b, 0, 1) - pgamma(run - isi, a, b, 0, 1);
	
	if(to_switch & !log_p)
		v = -expm1(lv);
	else if(to_switch & log_p)
		v = log(-expm1(lv));
	else if(!to_switch & !log_p)
		v = exp(lv);
	else 
		v = lv;
	
	return v;
	
}

// +++++++++++++ Block update of states +++++++++++++++++++
// +++++++ Uses forward-recursion-backward-sampling +++++++

void block_update(int *s, int i0, int m, double *r, double *par){
	
	int i, i1, j, k, kc;
	double gmaxg;
	
	isi_mean[0]  = exp(par[1]);    
	isi_shape[0] = exp(par[0]);
	isi_scale[0] = exp(par[1] - par[0]);
	isi_rate[0]  = exp(par[0] - par[1]);
	
	isi_mean[1]  = exp(par[3]);
	isi_shape[1] = exp(par[2]);
	isi_scale[1] = exp(par[3] - par[2]);
	isi_rate[1]  = exp(par[2] - par[3]);
	
	for(k = 0; k < 2; k++){
		
		dur_mean[k]  = exp(par[5 + 2 * k]);
		dur_shape[k] = exp(par[4 + 2 * k]);
		dur_scale[k] = exp(par[5 + 2 * k] - par[4 + 2 * k]);
		dur_rate[k]  = exp(par[4 + 2 * k] - par[5 + 2 * k]);
	}
	
	nconst[0] = pgamma(MaxBurstISI, isi_shape[0], isi_scale[0], 1, 0);
	lnconst[0] = pgamma(MaxBurstISI, isi_shape[0], isi_scale[0], 1, 1);
	
	nconst[1] = 1.0;
	lnconst[1] = 0.0;
	
	
	// Mass functions : forward recursion
	
	i1 = i0;
	if(i0 == 0){
		
		el[0][0][0] = y[i1];
		el[0][1][0] = y[i1];
		g[0][0][0] = log(0.5) + disi(y[i1], 0, 1); // dgamma(y[i1], isi_shape[0], isi_scale[0], 0);
		g[0][1][0] = log(0.5) + disi(y[i1], 1, 1); // dgamma(y[i1], isi_shape[1], isi_scale[1], 0);
		
	} else{
		
		el[0][0][0] = (s[i0 - 1] ^ 1) * r[i0 - 1] + y[i1];
		el[0][1][0] = (s[i0 - 1] ^ 0) * r[i0 - 1] + y[i1];
		g[0][0][0] = prob_move(r[i0 - 1], y[i0 - 1], dur_shape[s[i0 - 1]], dur_scale[s[i0 - 1]], s[i0 - 1] ^ 0, 1) + disi(y[i1], 0, 1);
		g[0][1][0] = prob_move(r[i0 - 1], y[i0 - 1], dur_shape[s[i0 - 1]], dur_scale[s[i0 - 1]], s[i0 - 1] ^ 1, 1) + disi(y[i1], 1, 1); 
		
	}
	
	gmax[0] = g[0][0][0];
	gmax[1] = g[0][1][0];
	
	gmaxg = maxvec(gmax, 2);
	gtot = log(exp(g[0][0][0] - gmaxg) + exp(g[0][1][0] - gmaxg));
	
	gmaxg += gtot;
	g[0][0][0] -= gmaxg; 
	g[0][1][0] -= gmaxg;
	
	gnorm[0][0] = exp(g[0][0][0]);
	gnorm[0][1] = exp(g[0][1][0]);
	
	for(i = 1; i < m; i++){		
		i1++;		

		for(k = 0; k < 2; k++){			
			kc = !k;
			for(j = 1; j <= i; j++)
				g[i][k][j] = g[i - 1][k][j - 1] + prob_move(el[i - 1][k][j - 1], y[i1 - 1], dur_shape[k], dur_scale[k], 0, 1) + disi(y[i1], k, 1);
			
			for(g[i][k][0] = 0.0, j = 0; j < i; j++)
				g[i][k][0] += exp(g[i - 1][kc][j] + prob_move(el[i - 1][kc][j], y[i1 - 1], dur_shape[kc], dur_scale[kc], 1, 1) + disi(y[i1], k, 1));
			g[i][k][0] = log(g[i][k][0]);
			
			
			el[i][k][0] = y[i1];
			for(j = 1; j <= i; j++)
				el[i][k][j] = el[i - 1][k][j - 1] + y[i1]; 
			
			gmax[k] = maxvec(g[i][k], i + 1);
		}
		
		gmaxg = maxvec(gmax, 2);
		for(gtot = 0.0, k = 0; k < 2; k++)
			for(j = 0; j <= i; j++)
				gtot += exp(g[i][k][j] - gmaxg);
		
		gtot = log(gtot);
		
		gmaxg += gtot;
		
		for(k = 0; k < 2; k++)
			for(j = 0; j <= i; j++)
				g[i][k][j] -= gmaxg;
		
		for(k = 0; k < 2; k++)
			for(gnorm[i][k] = 0.0, j = 0; j <= i; j++)
				gnorm[i][k] += exp(g[i][k][j]);
		
	}
	
	
	// Sample : backward recursion
	
	i = m - 1;
	i1 = i0 + i;
	
	if(i1 == n - 1){
		
		k = rdisc(2, gnorm[i], 0);
		j = rdisc(i + 1, g[i][k], 1);
		
		s[i1] = k;
		r[i1] = el[i][k][j];
		
	} else{
		
		k = !s[i1 + 1];
		s[i1] = k;
		for(j = 0; j <= i; j++)
			g[i][k][j] += prob_move(el[i][k][j], y[i1], dur_shape[k], dur_scale[k], 1, 1);
		j = rdisc(i + 1, g[i][k], 1);
		r[i1] = el[i][k][j];
	}
	
	
	while(i > 0){
		i--;
		i1--;
		if(j > 0){
			j--;
			s[i1] = k;
			r[i1] = el[i][k][j];
		} else{
			
			k = !k;
			s[i1] = k;
			for(j = 0; j <= i; j++)
				g[i][k][j] += prob_move(el[i][k][j], y[i1], dur_shape[k], dur_scale[k], 1, 1);
			j = rdisc(i + 1, g[i][k], 1);
			r[i1] = el[i][k][j];
		}
	}
}

// +++++ Update parameters given states via MTM +++++++

void par_update2(int *s, double *r, double *par, double *hpar, double *sc, double *logpost, double *acpt, int M){
	
	int i, k, k2, m, l, lpick, off;
	double u[2], incr[2] , par_old[2], par_new[2], unorm;
	double stepfwd[L], stepbwd[Lm], lpostnew[L], wtfwd, wtbwd, lpost_max;
	
	double a;
	for(m = 0; m < M; m++){
		
		for(k2 = 0; k2 < 4; k2++){
			
			for(k = 0; k < 2; k++)
				par_old[k] = par[2 * k2 + k];
			
			if(k2 > 1){
				u[0] = 0.0;
				u[1] = 1.0;
			}  else{
				u[0] = rnorm(0.0, 1.0); u[1] = rnorm(0.0, 1.0);
				unorm = norm2(u, 2);				
				for(k = 0; k < 2; k++)
					u[k] /= unorm;
			}
			//			Rprintvec("u = ", u, 2);
			
			for(l = 0; l < L; l++)
				stepfwd[l] = rnorm(0.0, 1.0);
			for(l = 0; l < Lm; l++)
				stepbwd[l] = rnorm(0.0, 1.0);
			
			
			for(l = 0; l < L; l++){				
				for(k = 0; k < 2; k++)
					incr[k] = sc[2 * k2 + k] * stepfwd[l] * u[k];				
				lpostnew[l] = lpost_incr(s, r, par_old, hpar, k2, incr);
			}
			
			lpost_max = maxvec(lpostnew, L);
			for(l = 0; l < L; l++)
				lpostnew[l] -= lpost_max;
			
			for(wtfwd = 0.0, l = 0; l < L; l++)
				wtfwd += exp(lpostnew[l]);
			
			lpick = rdisc(L, lpostnew, 1); 
			//			Rprintf("lpost_max = %g, lpostnew[lpick] = %g\n", lpost_max, lpostnew[lpick]);
			
			for(k = 0; k < 2; k++)
				par_new[k] = par_old[k] + sc[2 * k2 + k] * stepfwd[lpick] * u[k];
			
			wtbwd = exp(-lpost_max);
			for(k = 0; k < 2; k++)
				incr[k] = -sc[2 * k2 + k] * stepfwd[lpick] * u[k];
			a = lpost_incr(s, r, par_new, hpar, k2, incr);
			//			Rprintf("Should be = %g, is = %g\n", -lpostnew[lpick] - lpost_max, a);
			for(l = 0; l < Lm; l++){
				
				for(k = 0; k < 2; k++) 
					incr[k] = sc[2 * k2 + k] * stepbwd[l] * u[k];
				a = lpost_incr(s, r, par_new, hpar, k2, incr);
				wtbwd += exp(a + lpostnew[lpick]);
			}
			//      Rprintf("\nwtbwd = %g, wtfwd = %g\n\n", wtbwd, wtfwd);
			
			if(wtbwd * runif(0.0,1.0) < wtfwd){
				
				for(k = 0; k < 2; k++)
					par[2 * k2 + k] = par_new[k];
				logpost[0] += (lpostnew[lpick] + lpost_max);
				
				acpt[k2] += 1.0; 
			}
			acpt_ctr[k2] += 1.0;
		}
	}	
}

// ++++ los posterior density given states and parameters ++++

double lpost(int *s, double *r, double *par, double *hpar, int get_info){
	
	int i, k;
	double ll = 0.0, lp = 0.0, lp1, lp2;
	
	
	isi_shape[0] = exp(par[0]);        
	isi_scale[0] = exp(par[1] - par[0]);
	isi_rate[0]  = exp(par[0] - par[1]);
	
	isi_shape[1] = exp(par[2]);
	isi_scale[1] = exp(par[3] - par[2]);
	isi_rate[1]  = exp(par[2] - par[3]);
	
	nconst[0] = pgamma(MaxBurstISI, isi_shape[0], isi_scale[0], 1, 0);
	lnconst[0] = pgamma(MaxBurstISI, isi_shape[0], isi_scale[0], 1, 1);
	
	nconst[1] = 1.0;
	lnconst[1] = 0.0;
	
	for(k = 0; k < 2; k++){
		dur_shape[k] = exp(par[4 + 2 * k]);
		dur_scale[k] = exp(par[5 + 2 * k] - par[4 + 2 * k]);
		dur_rate[k]  = exp(par[4 + 2 * k] - par[5 + 2 * k]);
	}
	
	if(get_info){
		for(k = 0; k < 2; k++){
			len[k] = 0;
			tot_nisi[k] = 0;
			tot_isi[k] = 0.0;
			tot_lisi[k] = 0.0;
		}
		
		for(i = 1; i < n; i++){
			tot_nisi[s[i - 1]] += 1;
			tot_isi[s[i - 1]] += y[i - 1];
			tot_lisi[s[i - 1]] += log(y[i - 1]);
			
			if(s[i] ^ s[i - 1]){
				chpt[s[i - 1]][len[s[i - 1]]] = i - 1;
				len[s[i - 1]] += 1;
			}
			
		}
		
		tot_nisi[s[n - 1]] += 1;
		tot_isi[s[n - 1]] += y[n - 1];
		tot_lisi[s[n - 1]] += log(y[n - 1]);
		
	}
	
	
	for(k = 0; k < 2; k++){
		ll += (-tot_nisi[k] * (isi_shape[k] * log(isi_scale[k]) 
							   + lnconst[k] 
							   + lgamma(isi_shape[k])) 
			   + (isi_shape[k] - 1) * tot_lisi[k] 
			   - tot_isi[k] * isi_rate[k]);
		
		for(i = 0; i < len[k]; i++){
			ll += (pgamma(r[chpt[k][i]] - y[chpt[k][i]], dur_shape[k], dur_scale[k], 0, 0)
				   - pgamma(r[chpt[k][i]], dur_shape[k], dur_scale[k], 0, 0));
			
		}
		
	}
	
	
	for(k = 0; k < 8; k++)
		lp += dnorm(par[k], hpar[2 * k], hpar[1 + 2 * k], 1);
	
	return (ll + lp);
	
}

// ++++ Compute change in lpost for change in a pair of parameters +++++

double lpost_incr(int *s, double *r, double *par, double *hpar, int par_set, double *incr){
	
	int i;
	double ll_incr = 0.0, lp_incr = 0.0;
	double shape, scale, rate, shape_new, scale_new, rate_new;
	
	shape = exp(par[0]);        
	scale = exp(par[1] - par[0]);
	rate  = exp(par[0] - par[1]);
	shape_new = exp(par[0] + incr[0]);
	scale_new = exp(par[1] - par[0] + incr[1] - incr[0]);
	rate_new  = exp(par[0] - par[1] + incr[0] - incr[1]);
	
	switch(par_set){
			
		case 0 :
			
			ll_incr = (-tot_nisi[0] * (shape_new * log(scale_new)  - shape * log(scale) 
									   + lgamma(shape_new) - lgamma(shape)
									   + (pgamma(MaxBurstISI, shape_new, scale_new, 1, 1)
										  - pgamma(MaxBurstISI, shape, scale, 1, 1))
									   ) 
					   + (shape_new - shape) * tot_lisi[0] 
					   - tot_isi[0] * (rate_new - rate)
					   );
			
			lp_incr = (dnorm(par[0] + incr[0], hpar[0], hpar[1], 1) 
					   - dnorm(par[0], hpar[0], hpar[1], 1)
					   + dnorm(par[1] + incr[1], hpar[2], hpar[3], 1)
					   - dnorm(par[1], hpar[2], hpar[3], 1));
			
			break;
			
		case 1:    
			
			ll_incr = (-tot_nisi[1] * (shape_new * log(scale_new)  - shape * log(scale) 
									   + lgamma(shape_new) - lgamma(shape)) 
					   + (shape_new - shape) * tot_lisi[1] 
					   - tot_isi[1] * (rate_new - rate)
					   );
			
			
			lp_incr = (dnorm(par[0] + incr[0], hpar[4], hpar[5], 1) 
					   - dnorm(par[0], hpar[4], hpar[5], 1)
					   + dnorm(par[1] + incr[1], hpar[6], hpar[7], 1)
					   - dnorm(par[1], hpar[6], hpar[7], 1));
			
			break;
			
		case 2 :
			
			for(ll_incr = 0.0, i = 0; i < len[0]; i++)
				ll_incr += (log(pgamma(r[chpt[0][i]] - y[chpt[0][i]], shape_new, scale_new, 0, 0)
								- pgamma(r[chpt[0][i]], shape_new, scale_new, 0, 0))
							- log(pgamma(r[chpt[0][i]] - y[chpt[0][i]], shape, scale, 0, 0)
								  - pgamma(r[chpt[0][i]], shape, scale, 0, 0))
							);
			
			lp_incr = (dnorm(par[0] + incr[0], hpar[8], hpar[9], 1) 
					   - dnorm(par[0], hpar[8], hpar[9], 1)
					   + dnorm(par[1] + incr[1], hpar[10], hpar[11], 1)
					   - dnorm(par[1], hpar[10], hpar[11], 1));
			
			break;
			
		case 3 :
			
			for(ll_incr = 0.0, i = 0; i < len[1]; i++)
				ll_incr += (log(pgamma(r[chpt[1][i]] - y[chpt[1][i]], shape_new, scale_new, 0, 0)
								- pgamma(r[chpt[1][i]], shape_new, scale_new, 0, 0))
							- log(pgamma(r[chpt[1][i]] - y[chpt[1][i]], shape, scale, 0, 0)
								  - pgamma(r[chpt[1][i]], shape, scale, 0, 0))
							);
			
			lp_incr = (dnorm(par[0] + incr[0], hpar[12], hpar[13], 1) 
					   - dnorm(par[0], hpar[12], hpar[13], 1)
					   + dnorm(par[1] + incr[1], hpar[14], hpar[15], 1)
					   - dnorm(par[1], hpar[14], hpar[15], 1));
			
			break;
			
		default :
			
			Rprintf("Inappropriate par_set value\n");
	}
	
	if(ll_incr == atof("inf"))
		Rprintf("Warning: Inf increment in log posterior\n");
	
	return (ll_incr + lp_incr);
}

// +++++++++ Sample spike train fromm the model +++++++

void rspike(double *pars, double *horizon, int *stt, double *spikes, int *nspikes, double *transition, int *ntrans){
	
	int i, j, k, l, state = 0;
	double t = 0.0, T = horizon[0], t_delta;
	double shape_B_len = pars[0], scale_B_len = pars[1], shape_N_len = pars[2], scale_N_len = pars[3];
	double shape_B_isi = pars[4], scale_B_isi = pars[5], shape_N_isi = pars[6], scale_N_isi = pars[7];
	
	i = 0;
	state = 1;   // start with non bursting
	
	GetRNGstate();
	while(t < T){
		
		transition[i] = rgamma(pars[2 * state], pars[2 * state + 1]);
		state = !state;
		t += transition[i];
		i++;
	}  
	
	ntrans[0] = i;
	transition[ntrans[0] - 1] -= (t - T);
	
	state = 1;
	i = 0;
	j = 0;
	
	t = 0.0;
	double  resid = 0.0;
	while(j < ntrans[0]){
		
		t_delta = resid + rgamma(pars[2 * state + 4], pars[2 * state + 5]);
		
		while(t + t_delta < transition[j]){
			stt[i] = state;
			spikes[i++] = t_delta;
			t += t_delta;
			t_delta = rgamma(pars[2 * state + 4], pars[2 * state + 5]);
		}
		
		resid = transition[j] - t;
		t = -resid;
		state = !state;
		j++;
		
	}
	
	PutRNGstate();
	
	nspikes[0] = i;
	
}

// ++++++ density of isi distribution +++++++++

double disi(double isi_val, int isi_state, int get_log){
	
	double val;
	
	if(isi_state || (isi_val < MaxBurstISI))
		val = dgamma(isi_val, isi_shape[isi_state], isi_scale[isi_state], get_log);
	else
		val = get_log ? atof("-Inf") : 0;
	
	if(get_log)
		val -= lnconst[isi_state];
	else
		val /= nconst[isi_state];
	
	return val;
}



// +++++++ Max of a double vector ++++++++++++++

double maxvec(double *x, int n){
	int i;
	double xmax = x[0];
	for(i = 1; i < n; i++)
		if(x[i] > xmax)
			xmax = x[i];
	
	return xmax;
}

// ++++++ Sampling a discrete probability vector +++++++

int rdisc(int n, double *prob, int islog){
	
	int i, j, i0;
	double u, probmax, *csump;
	
	csump = (double *)R_alloc(n, sizeof(double));	
	if(islog){	
		probmax = maxvec(prob, n);
		csump[0] = exp(prob[0] - probmax);
		for(i = 1; i < n; i++) 
			csump[i] = csump[i - 1] + exp(prob[i] - probmax);
		
	} else{		
		csump[0] = prob[0];		
		for(i = 1; i < n; i++) 
			csump[i] = csump[i - 1] + prob[i];		
	}
	
	if(csump[n - 1] == 0)
		i0 = 0;
	else {		
		u = runif(0.0, 1.0);    
		i = 0;
		while(csump[i] <= u * csump[n - 1]) 
			i++;		
		i0 = i;
	}
	return i0;
}


// ++++++++++++ L-2 norm +++++++++++++++

double norm2(double *vec, int len){
	int i;
	double c_norm = 0.0;         
	for(i = 0; i < len; i++) 
		c_norm += pow(vec[i],2);      
	return sqrt(c_norm);   
}

// +++++++++++++++ print a vector +++++++++++++++

void Rprintvec(char *a, double *x, int len){
	
	int i;
	Rprintf("%s", a);
	for(i = 0; i < len; i++)
		Rprintf("%g ", x[i]);
	
	Rprintf("\n");
}


