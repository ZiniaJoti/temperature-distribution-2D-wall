function F = apply_neumann_edges(F, xn, yn, conn, bnd, qTop, qBottom, qLeft, qRight)
% Adds Neumann flux contributions to global F.
% Convention: -k dT/dn = qbar on boundary
% Weak form used in main: RHS gets -∫ w*qbar dΓ
% So we add fe_edge = ∫ N^T * (-qbar) dΓ

% Only applies where q is finite and nonzero.
tol = 1e-12;

% 1D 2-point Gauss
gp = [-1/sqrt(3), +1/sqrt(3)];
w  = [1, 1];

yTop = max(yn);
xRight = max(xn);

for e = 1:size(conn,1)
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);

    % local edges (node pairs):
    % bottom: 1-2 (eta=-1)
    % right : 2-3 (xi=+1)
    % top   : 4-3 (eta=+1)
    % left  : 1-4 (xi=-1)

    % bottom boundary
    if all(abs(ye([1 2]) - 0) < tol) && isfinite(qBottom) && qBottom ~= 0
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'bottom', -qBottom, gp, w);
    end

    % top boundary
    if all(abs(ye([4 3]) - yTop) < tol) && isfinite(qTop) && qTop ~= 0
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'top', -qTop, gp, w);
    end

    % left boundary
    if all(abs(xe([1 4]) - 0) < tol) && isfinite(qLeft) && qLeft ~= 0
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'left', -qLeft, gp, w);
    end

    % right boundary
    if all(abs(xe([2 3]) - xRight) < tol) && isfinite(qRight) && qRight ~= 0
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'right', -qRight, gp, w);
    end
end
end

function fe = edge_load_Q4(xe, ye, edge, qterm, gp, w)
% Edge load: fe += ∫ N^T * qterm dΓ  on the specified edge
% Returns 4x1 vector

fe = zeros(4,1);

for i = 1:2
    s = gp(i); wt = w(i);

    switch edge
        case 'bottom' % eta=-1, xi=s
            xi=s; eta=-1;
            [N, dN_dxi, ~] = shape_Q4(xi, eta);
            dx_dxi = dN_dxi * xe; dy_dxi = dN_dxi * ye;
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

        case 'top'    % eta=+1, xi=s
            xi=s; eta=+1;
            [N, dN_dxi, ~] = shape_Q4(xi, eta);
            dx_dxi = dN_dxi * xe; dy_dxi = dN_dxi * ye;
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

        case 'left'   % xi=-1, eta=s
            xi=-1; eta=s;
            [N, ~, dN_deta] = shape_Q4(xi, eta);
            dx_deta = dN_deta * xe; dy_deta = dN_deta * ye;
            Jedge = sqrt(dx_deta^2 + dy_deta^2);

        case 'right'  % xi=+1, eta=s
            xi=+1; eta=s;
            [N, ~, dN_deta] = shape_Q4(xi, eta);
            dx_deta = dN_deta * xe; dy_deta = dN_deta * ye;
            Jedge = sqrt(dx_deta^2 + dy_deta^2);

        otherwise
            error('Unknown edge type.');
    end

    fe = fe + (N' * qterm) * Jedge * wt;
end
end

function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
dN_dxi  = 0.25 * [-(1-eta), +(1-eta), +(1+eta), -(1+eta)];
dN_deta = 0.25 * [-(1-xi),  -(1+xi), +(1+xi),  +(1-xi)];
end

