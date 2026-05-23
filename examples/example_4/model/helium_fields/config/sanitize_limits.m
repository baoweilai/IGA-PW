function lim = sanitize_limits(lim)
%Validate plotting limits.
lim = real(lim(:)).';
lim = lim(isfinite(lim));
if numel(lim) < 2
    lim = [0, 1];
end
lim = [min(lim), max(lim)];
if lim(1) == lim(2)
    delta = max(1, abs(lim(1))) * 1e-6;
    lim = [lim(1) - delta, lim(2) + delta];
end
end
