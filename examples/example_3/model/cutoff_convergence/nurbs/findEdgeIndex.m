function dg_ele_idx = findEdgeIndex(uBreaks,n_edges, knot)
% Find the edge element containing a knot value.

low = 1;
high = n_edges+1;

mid = fix((low+high)/2);

while knot<uBreaks(mid) || knot>uBreaks(mid+1)
    if knot<uBreaks(mid)
        high = mid;
    else
        low = mid;
    end

    mid = fix((low+high)/2);


end

dg_ele_idx = mid;


end
