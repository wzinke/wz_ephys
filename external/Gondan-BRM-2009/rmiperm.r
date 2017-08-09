#
# Permutation test for the race model inequality
#
# Tested with R-2.9.1, download at http://www.stats.bris.ac.uk/R/
#

# tmax statistic

tmax = function(di)
{
    # t-Tests for each percentile
    tvalues = rowMeans(di) / apply(di, 1, sd) * sqrt(ncol(di))
    max(tvalues)
} # tmax

# The significance test
#
# - d: list including data of all participants
# - quantiles: quantiles of AV distribution at which the RMI is to be tested
# - nperm: number of random permutations used to generate H0 distribution
# - return value: tmax, critical value for one-sided test, P value

rmi_perm = function(d, quantiles=c(.05, .10, .15, .20, .25, .30), nperm=10001)
{
    # Determine quantiles from entire response time distribution 
    tmix     = lapply(d, unlist) 
    q        = lapply(tmix, quantile, probs=quantiles, type=4)

    # Evaluate distribution functions at q
    FA       = mapply(function(ti, qi) ecdf(ti$A)(qi), d, q)
    FV       = mapply(function(ti, qi) ecdf(ti$V)(qi), d, q)
    FAV      = mapply(function(ti, qi) ecdf(ti$AV)(qi), d, q)

    # Determine di = FAV(q) - FA(q) - FV(q)
    di       = matrix(FAV - FA - FV, nrow=length(quantiles), ncol=length(d))

    # Observed test statistic
    stat     = tmax(di)

    # Permutation distribution of test statistic
    stati = numeric(nperm)
    for(i in 1:nperm)
        stati[i] = tmax(di %*% diag(sign(rnorm(d))))

    list(tmax=stat, tcrit=quantile(stati, 0.95), P=mean(stat <= stati))
} # rmi_perm

# Read in data from text file
#
# obs: observer number, cond: condition, rt: response time (Inf: omitted response)
#
# obs cond   rt  # This header should be included in the file
#   1    A  245  # response to A
#   1    V  281  # response to V
#   1   AV  212  # response to AV
#   1    V  Inf  # omitted response

rmi_read = function(fname)
{
    d = read.table(fname, header=TRUE)
    obs = split(d[, c('cond', 'rt')], d$obs)
    for(i in 1:length(obs))
        obs[[i]] = split(obs[[i]]$rt, obs[[i]]$cond)
    obs
} # rmi_read

# Example session with data in C:\RTDATA\go2004.dat

data = rmi_read('C:/RTDATA/go2004.dat')
rmi_perm(data, quantiles=c(.05, .10, .15, .20, .25, .30), nperm=10001)
