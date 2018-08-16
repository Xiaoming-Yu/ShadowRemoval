function [r] = sat(d, l, u)
    r = d;
    r(d<l) = l;
    r(d>u) = u;
end