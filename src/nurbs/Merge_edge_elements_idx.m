function edge_element_idx =  Merge_edge_elements_idx(uBreaks,uBreaks_merge)
%Merge element indices along shared edges.

n_merge_edges = length(uBreaks_merge) - 1;
n_edges = length(uBreaks) - 1;

edge_element_idx = zeros(n_merge_edges,1);

for i = 1:n_merge_edges
    a = uBreaks_merge(i);
    b = uBreaks_merge(i+1);
    mid = (a+b)/2;
    dg_ele_idx = findEdgeIndex(uBreaks,n_edges,mid);
    edge_element_idx(i) = dg_ele_idx;
end


end
