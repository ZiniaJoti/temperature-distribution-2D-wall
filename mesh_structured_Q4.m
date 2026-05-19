function [xn, yn, conn, bnd] = mesh_structured_Q4(L, H, nx, ny)
% Structured Q4 mesh on [0,L]x[0,H]
% Returns:
%   xn, yn : node coordinates (column vectors)
%   conn   : (ne x 4) connectivity [n1 n2 n3 n4]
%   bnd    : structure of boundary node sets

x = linspace(0, L, nx+1);
y = linspace(0, H, ny+1);
[X, Y] = meshgrid(x, y);

xn = X(:);
yn = Y(:);

conn = zeros(nx*ny, 4);
e = 0;
for j = 1:ny
    for i = 1:nx
        e = e + 1;
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n4 = j*(nx+1) + i;
        n3 = n4 + 1;
        conn(e,:) = [n1 n2 n3 n4];
    end
end

tol = 1e-12;
bnd.left   = find(abs(xn - 0) < tol);
bnd.right  = find(abs(xn - L) < tol);
bnd.bottom = find(abs(yn - 0) < tol);
bnd.top    = find(abs(yn - H) < tol);

end
