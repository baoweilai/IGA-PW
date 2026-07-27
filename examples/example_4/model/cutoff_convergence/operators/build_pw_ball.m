function [k_list, n_basis] = build_pw_ball(Nc)
% Build the spherical plane-wave index set.

N = floor(Nc);

% Preallocate with a rough upper bound
k_list = zeros((2*N+1)^3, 3);
n_basis = 0;

for k1 = -N:N
    for k2 = -N:N
        rem2 = N^2 - k1^2 - k2^2;
        if rem2 < 0
            continue;
        end

        m = floor(sqrt(rem2));
        for k3 = -m:m
            n_basis = n_basis + 1;
            k_list(n_basis,:) = [k1, k2, k3];
        end
    end
end

k_list = k_list(1:n_basis,:);

end
