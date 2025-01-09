function temp_Gene = Eg2_ifout(c, ub, lb)
if c < lb || c > ub
    if abs(c - lb) < abs(c - ub)
        temp_Gene = lb;
    else
        temp_Gene = ub;
    end
else
    temp_Gene = c;
end