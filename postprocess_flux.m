function [qx, qy, xc, yc] = postprocess_flux(xn, yn, conn, T, L1, k1, k2)
% Computes heat flux q = -k * grad(T) at element centers.

ne = size(conn,1);
qx = zeros(ne,1);
qy = zeros(ne,1);
xc = zeros(ne,1);
yc = zeros(ne,1);

xi = 0; eta = 0; % element center
[~, dN_dxi, dN_deta] = shape_Q4(xi, eta);

for e = 1:ne
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);
    Te = T(nodes);

    xce = mean(xe); yce = mean(ye);
    xc(e) = xce; yc(e) = yce;

    if xce <= L1
        k = k1;
    else
        k = k2;
    end

    J = [dN_dxi; dN_deta] * [xe(:) ye(:)];

    dN = J\ [dN_dxi; dN_deta];  % 2x4
    dTdx = dN(1,:) * Te;
    dTdy = dN(2,:) * Te;

    qx(e) = -k * dTdx;
    qy(e) = -k * dTdy;
end
end

function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];

% row vectors (1x4)
dN_dxi  = 0.25 * [-(1-eta), +(1-eta), +(1+eta), -(1+eta)];
dN_deta = 0.25 * [-(1-xi),  -(1+xi), +(1+xi),  +(1-xi)];
end
