function [Ubar, wbar] = IGAknotRefineCurve(U, w, p, refinementCount)
% Refine a NURBS curve knot vector.

if refinementCount == 0
    Ubar = U;
    wbar = w;
else
    for refinementStep = 1:refinementCount
        uBreaks = unique(U);
        n = length(w);
        X = (uBreaks(2:end) + uBreaks(1:end-1)) / 2;

        Ubar = sort([U, X]);
        r = length(X);
        wbar = zeros(1, n + r);

        a = findspan(U, p, X(1));
        b = findspan(U, p, X(end)) + 1;
        wbar(1:a-p) = w(1:a-p);
        wbar(b+r-1:n+r) = w(b-1:n);

        knotIndex = b + p - 1;
        outputIndex = b + p + r - 1;
        for insertionIndex = r:-1:1
            while X(insertionIndex) <= U(knotIndex) && knotIndex > a
                wbar(outputIndex-p-1) = w(knotIndex-p-1);
                outputIndex = outputIndex - 1;
                knotIndex = knotIndex - 1;
            end
            wbar(outputIndex-p-1) = wbar(outputIndex-p);

            for degreeIndex = 1:p
                weightIndex = outputIndex - p + degreeIndex;
                alpha = Ubar(outputIndex + degreeIndex) - X(insertionIndex);
                if alpha == 0
                    wbar(weightIndex-1) = wbar(weightIndex);
                else
                    alpha = alpha / ...
                        (Ubar(outputIndex + degreeIndex) ...
                        - U(knotIndex - p + degreeIndex));
                    wbar(weightIndex-1) = alpha * wbar(weightIndex-1) ...
                        + (1 - alpha) * wbar(weightIndex);
                end
            end
            outputIndex = outputIndex - 1;
        end
        U = Ubar;
        w = wbar;
    end
end
