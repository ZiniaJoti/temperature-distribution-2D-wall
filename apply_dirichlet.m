function T = apply_dirichlet(K, F, dirNodes, dirVals)
% Dirichlet enforcement using elimination.
% Returns full temperature vector T.

nn = size(K,1);
T = zeros(nn,1);

% Known temperatures
T(dirNodes) = dirVals;

% Modify RHS: F' = F - K(:,dir)*Tdir
Fmod = F - K(:, dirNodes) * T(dirNodes);

free = true(nn,1);
free(dirNodes) = false;

Kff = K(free, free);
Ff  = Fmod(free);

Tf = Kff \ Ff;

T(free) = Tf;
end
