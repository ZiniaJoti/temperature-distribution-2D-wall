function main_composite_wall_FEM_Q4()
clc; clear; close all;

%% ---------------------------
%  USER INPUTS (edit here)
%% ---------------------------
L  = 1.0;          % wall length (x)
H  = 1.0;          % wall height (y)
L1 = 0.5;          % interface location (x = L1)

k1 = 10.0;         % material 1 (brick)
k2 = 1.0;          % material 2 (insulation)

Tin  = 0;          % T at x = L   (inside)
Tout = 100;        % T at x = 0   (outside)
amp = 20;

nx = 30;           % number of elements in x
ny = 30;           % number of elements in y

S  = 50.0;          % source term (volumetric), keep 0 for now

% Neumann flux on boundaries (engineering BC: -k dT/dn = qbar)
qTop    = 0.0;     % at y = H
qBottom = 50.0;     % at y = 0

%% ---------------------------
%  MESH (structured Q4)
%% ---------------------------
[xn, yn, conn] = mesh_Q4_structured(L, H, nx, ny);
nn = numel(xn);            % number of nodes
ne = size(conn,1);         % number of elements

K = sparse(nn, nn);
F = zeros(nn,1);

%% ---------------------------
%  ASSEMBLY
%% ---------------------------
for e = 1:ne
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);

    % choose material based on element center x
    xc = mean(xe);
    if xc <= L1
        k = k1;
    else
        k = k2;
    end

    [Ke, Fe] = element_Q4_Ke_Fe(xe, ye, k, S);

    K(nodes, nodes) = K(nodes, nodes) + Ke;
    F(nodes) = F(nodes) + Fe;
end

%% ---------------------------
%  NEUMANN BC (top and bottom)
%% ---------------------------
% add flux contributions on y=0 and y=H edges
F = apply_neumann_top_bottom(F, xn, yn, conn, qTop, qBottom);

%% ---------------------------
%  DIRICHLET BC (left and right)
%% ---------------------------
tol = 1e-12;
leftNodes  = find(abs(xn - 0) < tol);
rightNodes = find(abs(xn - L) < tol);

dirNodes = [leftNodes; rightNodes];
%dirVals  = [Tout*ones(numel(leftNodes),1); Tin*ones(numel(rightNodes),1)];

Tleft  = Tout + amp * (yn(leftNodes) / H);   % varies linearly with y
Tright = Tin * ones(numel(rightNodes),1);   % keep right constant

dirVals = [Tleft; Tright];


[T, Kmod, Fmod] = apply_dirichlet_elimination(K, F, dirNodes, dirVals);

%% ---------------------------
%  SOLVE
%% ---------------------------
T = Kmod \ Fmod;
Tfree = Kmod \ Fmod;

% reconstruct full T (including Dirichlet values)
Tfull = zeros(nn,1);
Tfull(dirNodes) = dirVals;

free = true(nn,1);
free(dirNodes) = false;
Tfull(free) = Tfree;

T = Tfull;
%% ---------------------------
%  POSTPROCESSING: plots
%% ---------------------------
plot_temperature_contour(xn, yn, conn, T);
plot_midline(xn, yn, T, H/2);
plot_flux_quiver(xn, yn, conn, T, L1, k1, k2);

end

%% ==========================================================
%  Structured Q4 mesh generator
%% ==========================================================
function [xn, yn, conn] = mesh_Q4_structured(L, H, nx, ny)
% nodes on a (nx+1) x (ny+1) grid
x = linspace(0, L, nx+1);
y = linspace(0, H, ny+1);
[X, Y] = meshgrid(x, y);

xn = X(:);
yn = Y(:);

% connectivity: element (i,j) uses nodes:
% n1=(i,j), n2=(i+1,j), n3=(i+1,j+1), n4=(i,j+1)
conn = zeros(nx*ny, 4);
elem = 0;
for j = 1:ny
    for i = 1:nx
        elem = elem + 1;
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n4 = j*(nx+1) + i;
        n3 = n4 + 1;
        conn(elem,:) = [n1 n2 n3 n4];
    end
end
end

%% ==========================================================
%  Q4 element stiffness (2x2 Gauss) and load (source S)
%% ==========================================================
function [Ke, Fe] = element_Q4_Ke_Fe(xe, ye, k, S)
Ke = zeros(4,4);
Fe = zeros(4,1);

% 2x2 Gauss points
gp = [-1/sqrt(3), +1/sqrt(3)];
w  = [1, 1];

for a = 1:2
    for b = 1:2
        xi  = gp(a);
        eta = gp(b);
        wt  = w(a)*w(b);

        [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

        % Jacobian
        J = [dN_dxi; dN_deta] * [xe ye];  % 2x2
        detJ = det(J);
        invJ = inv(J);

        % derivatives w.r.t x,y: [dN/dx; dN/dy] = invJ * [dN/dxi; dN/deta]
        dN = J\[dN_dxi; dN_deta];   % 2x4
        dN_dx = dN(1,:); dN_dy = dN(2,:);

        % B matrix for scalar diffusion: grad(T) = [dN_dx; dN_dy] * Te
        B = [dN_dx; dN_dy]; % 2x4

        Ke = Ke + (B' * (k*eye(2)) * B) * detJ * wt;

        % source term contribution: ∫ N^T * (-S?) depends on your convention.
        % Here we use PDE: d/dx(kTx)+d/dy(kTy)+S=0 => RHS = -∫ N S dΩ moved to RHS?
        % For simplicity with S=0, it doesn't matter.
        Fe = Fe - (N' * S) * detJ * wt;
    end
end
end

%% ==========================================================
%  Q4 shape functions in parent coordinates
%% ==========================================================
function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta), ...
            (1+xi)*(1-eta), ...
            (1+xi)*(1+eta), ...
            (1-xi)*(1+eta)];

dN_dxi  = 0.25 * [-(1-eta), +(1-eta), +(1+eta), -(1+eta)];
dN_deta = 0.25 * [-(1-xi),  -(1+xi), +(1+xi),  +(1-xi)];
end

%% ==========================================================
%  Apply Dirichlet BC by elimination (strong enforcement)
%% ==========================================================
function [T, Kmod, Fmod] = apply_dirichlet_elimination(K, F, dirNodes, dirVals)
nn = size(K,1);
T = zeros(nn,1);

Kmod = K;
Fmod = F;

% set known temperatures into T
T(dirNodes) = dirVals;

% modify RHS: F = F - K(:,dir)*Tdir
Fmod = Fmod - Kmod(:, dirNodes) * T(dirNodes);

% eliminate rows/cols
free = true(nn,1);
free(dirNodes) = false;

Kmod = Kmod(free, free);
Fmod = Fmod(free);

% store mapping back is done in main after solve
% (we return partial system; main solves and inserts)
end

%% ==========================================================
%  Neumann flux on top and bottom edges (qTop, qBottom)
%  Uses 2-point Gauss on each boundary edge
%% ==========================================================
function F = apply_neumann_top_bottom(F, xn, yn, conn, qTop, qBottom)
tol = 1e-12;

% Gauss points for 1D edge integration
gp = [-1/sqrt(3), +1/sqrt(3)];
w  = [1, 1];

for e = 1:size(conn,1)
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);

    % edges in local node order:
    % bottom edge: (1-2), top edge: (4-3) in our connectivity [1 2 3 4]
    % Check if edge lies on y=0 or y=H by looking at y-values of edge nodes
    y12 = ye([1 2]);
    y43 = ye([4 3]);

    % bottom boundary (y=0)
    if all(abs(y12 - 0) < tol) && qBottom ~= 0
        % integrate ∫ N^T * ( -qbar ) dΓ because weak form used -∫ w qbar
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'bottom', -qBottom, gp, w);
    end

    % top boundary (y=H) : detect by max y of the edge
    yTop = max(yn);
    if all(abs(y43 - yTop) < tol) && qTop ~= 0
        F(nodes) = F(nodes) + edge_load_Q4(xe, ye, 'top', -qTop, gp, w);
    end
end
end

function fe_edge = edge_load_Q4(xe, ye, whichEdge, qterm, gp, w)
% returns a 4x1 vector (only the edge nodes get nonzero contributions)
fe_edge = zeros(4,1);

for i = 1:2
    s = gp(i); wt = w(i);

    switch whichEdge
        case 'bottom' % eta = -1, xi = s
            xi = s; eta = -1;
            [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

            % edge length factor: |dx/dxi| along edge
            % For eta fixed, tangent = d(x,y)/d xi
            dx_dxi = dN_dxi * xe;
            dy_dxi = dN_dxi * ye;
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

        case 'top'    % eta = +1, xi = s
            xi = s; eta = +1;
            [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

            dx_dxi = dN_dxi * xe;
            dy_dxi = dN_dxi * ye;
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

        otherwise
            error('Unknown edge');
    end

    fe_edge = fe_edge + (N' * qterm) * Jedge * wt;
end
end

%% ==========================================================
%  Plot temperature contour
%% ==========================================================
function plot_temperature_contour(xn, yn, conn, T)
figure;
patch('Faces', conn, 'Vertices', [xn yn], 'FaceVertexCData', T, ...
      'FaceColor','interp','EdgeColor',[0.6 0.6 0.6]);
colorbar; axis equal tight;
xlabel('x'); ylabel('y');
title('Temperature distribution T(x,y)');
end

%% ==========================================================
%  Plot midline temperature
%% ==========================================================
function plot_midline(xn, yn, T, ymid)
tol = 1e-6;
idx = find(abs(yn - ymid) < tol);
[xx, order] = sort(xn(idx));
TT = T(idx);
TT = TT(order);

figure;
plot(xx, TT, 'o-'); grid on;
xlabel('x'); ylabel('T');
title(sprintf('Midline temperature (y = %.3f)', ymid));
end

%% ==========================================================
%  Compute and plot heat flux quiver (element-centered)
%% ==========================================================
function plot_flux_quiver(xn, yn, conn, T, L1, k1, k2)
ne = size(conn,1);
xc = zeros(ne,1); yc = zeros(ne,1);
qx = zeros(ne,1); qy = zeros(ne,1);

% evaluate flux at element center (xi=0,eta=0)
xi=0; eta=0;
[N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

for e=1:ne
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);
    Te = T(nodes);

    xce = mean(xe); yce = mean(ye);
    xc(e)=xce; yc(e)=yce;

    k = (xce <= L1)*k1 + (xce > L1)*k2;

    J = [dN_dxi; dN_deta] * [xe ye];
    invJ = inv(J);
    dN = J\[dN_dxi; dN_deta];

    dTdx = dN(1,:)*Te;
    dTdy = dN(2,:)*Te;

    % q = -k * grad(T)
    qx(e) = -k*dTdx;
    qy(e) = -k*dTdy;
end

figure;
quiver(xc, yc, qx, qy);
axis equal tight;
xlabel('x'); ylabel('y');
title('Heat flux vectors q = -k [dT/dx, dT/dy]');
grid on;
end
