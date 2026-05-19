function main_composite_wall_2D()
clc; clear; close all;

%% =========================
%  PARAMETERS
%% =========================
L  = 1.0;        % length in x
H  = 1.0;        % height in y
L1 = 0.5;        % interface location (x = L1)

k1 = 10.0;       % brick conductivity
k2 = 1.0;        % insulation conductivity

Tout = 100;      % T at x=0 (outside)
Tin  = 0;        % T at x=L (inside)

nx = 30;         % elements in x
ny = 30;         % elements in y

S  = 0.0;        % source term in PDE: d/dx(kTx)+d/dy(kTy)+S=0

% Neumann BC convention in this project:
% -k * dT/dn = qbar  on Gamma_N
qTop    = 0.0;   % y=H (insulated if 0)
qBottom = 0.0;   % y=0 (insulated if 0)
qLeft   = NaN;   % keep NaN because left is Dirichlet
qRight  = NaN;   % keep NaN because right is Dirichlet

outdir = pwd;    % save figures in current folder

%% =========================
%  MESH
%% =========================
[xn, yn, conn, bnd] = mesh_structured_Q4(L, H, nx, ny);
nn = numel(xn);
ne = size(conn,1);

%% =========================
%  ASSEMBLY
%% =========================
K = sparse(nn, nn);
F = zeros(nn,1);

for e = 1:ne
    nodes = conn(e,:);
    xe = xn(nodes); ye = yn(nodes);

    % material selection by element center
    xc = mean(xe);
    if xc <= L1
        k = k1;
    else
        k = k2;
    end

    Ke = element_Q4_K(xe, ye, k);

    % source load vector (optional)
    Fe = zeros(4,1);
    if S ~= 0
        Fe = element_Q4_source(xe, ye, S);
    end

    K(nodes, nodes) = K(nodes, nodes) + Ke;
    F(nodes) = F(nodes) + Fe;
end

%% =========================
%  NEUMANN BC (adds to F)
%% =========================
F = apply_neumann_edges(F, xn, yn, conn, bnd, qTop, qBottom, qLeft, qRight);

%% =========================
%  DIRICHLET BC (strong)
%% =========================
dirNodes = [bnd.left(:); bnd.right(:)];
dirVals  = [Tout*ones(numel(bnd.left),1); Tin*ones(numel(bnd.right),1)];

T = apply_dirichlet(K, F, dirNodes, dirVals);

%% =========================
%  POSTPROCESS
%% =========================
[qx, qy, xc, yc] = postprocess_flux(xn, yn, conn, T, L1, k1, k2);

%% =========================
%  PLOTS 
%% =========================
plot_results(xn, yn, conn, T, xc, yc, qx, qy, H/2, outdir);

disp('Done.');

disp('Saved: T_contour.png, flux_quiver.png, midline_plot.png');
end

%% ---------------------------------------------------------
% Element consistent source vector for S
% PDE used: d/dx(kTx)+d/dy(kTy)+S=0
% Weak form gives RHS contribution: -∫ N^T * S dΩ
%% ---------------------------------------------------------
function Fe = element_Q4_source(xe, ye, S)
Fe = zeros(4,1);

gp = [-1/sqrt(3), +1/sqrt(3)];
w  = [1, 1];

for a = 1:2
    for b = 1:2
        xi = gp(a); eta = gp(b); wt = w(a)*w(b);

        [N, dN_dxi, dN_deta] = shape_Q4(xi, eta);

        J = [dN_dxi'; dN_deta'] * [xe ye];
        detJ = det(J);

        Fe = Fe - (N' * S) * detJ * wt;
    end
end
end

function [N, dN_dxi, dN_deta] = shape_Q4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
dN_dxi  = 0.25 * [-(1-eta), +(1-eta), +(1+eta), -(1+eta)];
dN_deta = 0.25 * [-(1-xi),  -(1+xi), +(1+xi),  +(1-xi)];
end
