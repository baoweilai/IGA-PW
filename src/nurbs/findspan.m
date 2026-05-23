function i = findspan(U, p, u)
%Locate an index or object used by the computation.

m = length(U);
n = m - p - 1;

if u < U(1) || u > U(m)
    error('Parameter value is outside the knot vector range.');
end

if u == U(n + 1)
    i = n;
    return;
end

low = p + 1;
high = n + 1;
mid = fix((high + low) / 2);

while u < U(mid) || u >= U(mid + 1)
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = fix((high + low) / 2);
end

i = mid;
end
