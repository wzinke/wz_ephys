random_group_test = function(dt1, dt2, Nrep=1000){
# t-test as permutation test
# from http://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r
    require(boot)

    vals  = c(dt1, dt2)
    group = c(rep('1',length(dt1)), rep('2',length(dt2)))
    grouptest = data.frame(vals, group)

    diff = function(d1,i){
        d = d1;
        d$group <- d$group[i];  # randomly re-assign groups
        Mean= tapply(X=d$vals, INDEX=d$group, mean)
        Diff = Mean[1]-Mean[2]
        Diff
    }

    rsmplboot = boot(data = grouptest, statistic = diff, R = Nrep)
    mean(abs(rsmplboot$t) > abs(rsmplboot$t0))
}
