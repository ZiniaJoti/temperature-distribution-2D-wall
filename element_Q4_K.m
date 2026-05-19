function Ke = element_Q4_K(xe, ye, k)
% Element stiffness for scalar diffusion (heat conduction)
% Ke = ∫ (B^T * kI * B) dΩ  using 2x2 Gauss quadrature

Ke = zeros(4,4);

gp = [-1/sqrt(3), +1/sqrt(3)];
w  = [1, 1];

for a = 1:2
    for b = 1:2
        xi = gp(a); eta = gp(b); wt = w(a)*w(b);

        [~, dN_dxi, dN_deta] = shape_Q4(xi, eta);
        % Jacobian 2x2
        J = [dN_dxi; dN_deta] * [xe(:) ye(:)]; 
        detJ = det(J);
        dN = J\[dN_dxi; dN_deta];    % 2x4
        dN_dx = dN(1,:);
        dN_dy = dN(2,:);

        B = [dN_dx; dN_dy];                 % 2x4

        Ke = Ke + (B'*(k*eye(2))*B) * detJ * wt;
    end
end
end

function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];

% derivatives are ROW vectors (1x4)
dN_dxi  = 0.25 * [-(1-eta), +(1-eta), +(1+eta), -(1+eta)];
dN_deta = 0.25 * [-(1-xi),  -(1+xi), +(1+xi),  +(1-xi)];
end
