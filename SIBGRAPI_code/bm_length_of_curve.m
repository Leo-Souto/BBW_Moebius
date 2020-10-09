function length = bm_length_of_curve(curve, j)
    dx = curve(2:j,1) - curve(1:j-1,1);
    dy = curve(2:j,2) - curve(1:j-1,2);
    ds = sqrt(dx.^2 + dy.^2);
    length = sum(ds);
end