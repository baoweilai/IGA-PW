function cmap = error_colormap(n)
%Return the error colormap.
arguments
    n = 256
end
anchors255 = [
     38,  82,  92;
    103, 170, 158;
    248, 242, 232;
    231, 151, 120;
    178,  56,  48
];
A = anchors255 / 255;
x = linspace(0, 1, size(A, 1));
xi = linspace(0, 1, n);
cmap = zeros(n, 3);
for j = 1:3
    cmap(:, j) = interp1(x, A(:, j), xi, 'pchip');
end
cmap = min(max(cmap, 0), 1);
end
