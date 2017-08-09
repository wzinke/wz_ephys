function vec = ensure_row(vec)
    if size(vec, 1) ~= 1
        vec = vec';
    end

