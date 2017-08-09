function ctr = wz_get_contrast(a,b)
    ctr = (a - b) ./ (a + b);
