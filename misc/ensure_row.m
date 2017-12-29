function vec = ensure_row(vec)
% if column vector transpose to row vector

    if size(vec, 1) ~= 1
        vec = transp(vec);
    end

