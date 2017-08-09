function f = wz_exgauss_mod(X, Mu, Sigma, Tau)
% modified from: https://github.com/bramzandbelt/exgauss


f     = (1./Tau).* ...
        exp(((Mu - X)./Tau) + ((Sigma.^2)./(2.*Tau.^2))).* ...
        .5.*(1+erf((((X-Mu)./Sigma) - (Sigma./Tau))./sqrt(2)));
