function plot_results(xn, yn, conn, T, xc, yc, qx, qy, ymid, outdir)

%% Temperature contour
fig1 = figure;
patch('Faces', conn, 'Vertices', [xn yn], ...
      'FaceVertexCData', T, 'FaceColor','interp', ...
      'EdgeColor',[0.7 0.7 0.7]);
axis equal tight;
xlabel('x'); ylabel('y');
title('Temperature distribution T(x,y)');
colorbar;

saveas(fig1, fullfile(outdir, 'T_contour.png'));

%% Flux quiver
fig2 = figure;
quiver(xc, yc, qx, qy);
axis equal tight;
xlabel('x'); ylabel('y');
title('Heat flux vectors q = -k [dT/dx, dT/dy]');
grid on;

saveas(fig2, fullfile(outdir, 'flux_quiver.png'));

%% Midline plot (closest nodes to y = ymid)
tol = 1e-6;
idx = find(abs(yn - ymid) < tol);
if isempty(idx)
    [~, j] = min(abs(yn - ymid));
    yline = yn(j);
    idx = find(abs(yn - yline) < tol);
    ymid = yline;
end

[xx, order] = sort(xn(idx));
TT = T(idx);
TT = TT(order);

fig3 = figure;
plot(xx, TT, 'o-'); grid on;
xlabel('x'); ylabel('T');
title(sprintf('Midline temperature, y = %.4f', ymid));

saveas(fig3, fullfile(outdir, 'midline_plot.png'));

end
