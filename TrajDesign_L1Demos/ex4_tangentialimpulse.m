clear; clc; close all;
p = cw_params(400);                         % 400 km circular target orbit
fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
        p.alt/1e3, p.n, p.T/60);

% --- tangential impulse at the target -------------------------------------
% A pure +dVx on V-bar, applied AT the target itself.  The post-impulse
% state is therefore  x0 = y0 = z0 = 0,  xd0 = dVx,  yd0 = zd0 = 0.
% Propagating that state with the CW solution reproduces the note's
%   x(t) = (dVx/n)(4 sin nt - 3 n t)
%   z(t) = (2 dVx/n)(cos nt - 1)
dVx = 0.20;                                 % tangential impulse [m/s]
s0  = [0, 0, 0, dVx, 0, 0];                 % [x y z xd yd zd]

% --- propagate ------------------------------------------------------------
t = linspace(0, 2*p.T, 600);                % two orbital periods
S = cw_propagate3d(s0, t, p.n);

x = S(1,:);                                 % V-bar
y = S(2,:);                                 % H-bar
z = S(3,:);                                 % R-bar

% --- closed-form solution from the note's Eq. (tang-impulse) --------------
% x(t) = (dVx/n)(4 sin nt - 3 n t),   z(t) = (2 dVx/n)(cos nt - 1)
x_formula = s0(1)+ (dVx/p.n) * (4*sin(p.n*t) - 3*p.n*t);
z_formula = s0(3)+ (2*dVx/p.n) * (cos(p.n*t) - 1);

% --- decompose x(t) into its two physical terms ---------------------------
% x(t) = (4 dVx/n) sin(nt)   +   (-3 dVx) t
%        \___ bounded hump ___/   \_ secular drift _/
x_osc    = (4*dVx/p.n) * sin(p.n*t);        % cycloidal hump
x_secular = -3*dVx * t;                      % along-track drift line

% --- plot: trajectory + closed-form overlay -------------------------------
figure; hold on; grid on; box on;
plot(x, z, 'LineWidth', 3.0, 'Color', [0.70 0.82 0.95]);
% closed-form curve plotted on top with a dashed marker line
plot(x_formula, z_formula, '--', 'LineWidth',1.6, 'Color',[0.85 0.1 0.1]);
plot(x(1),   z(1),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
plot(x(end), z(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);
plot(x_secular, zeros(size(t)), ':', 'Color',[.6 .6 .6], 'LineWidth',1.2);
plot(0, 0, 'ks', 'MarkerSize',12, 'MarkerFaceColor','none', 'LineWidth',1.6);
text(0, 0, '  target', 'FontSize',8, 'VerticalAlignment','bottom');
set(gca,'YDir','reverse');                  % +z (R-bar, toward Earth) down
axis equal;
xlabel('V-bar   x   [m]   (+ = ahead)');
ylabel('R-bar   z   [m]   (+ = toward Earth)');
title(sprintf('Tangential impulse  \\DeltaV_x = %.2f m/s', dVx), ...
      'FontWeight','bold');
legend({'propagated (cw\_propagate3d)','closed-form Eq.', ...
        'impulse point','after 2 orbits','secular drift line'}, ...
       'Location','best');

% mark one full orbit
[~,iT] = min(abs(t - p.T));
plot(x(iT), z(iT), 'ks', 'MarkerFaceColor','y', 'MarkerSize',8);
text(x(iT), z(iT), '  after 1 orbit', 'FontSize',8);

% --- second figure: the x(t) decomposition -------------------------------
figure('Color','w'); hold on; grid on; box on;
plot(t/p.T, x,         'LineWidth',2.2, 'Color',[0 0.45 0.74]);
plot(t/p.T, x_secular, '--', 'LineWidth',1.6, 'Color',[0.85 0.1 0.1]);
plot(t/p.T, x_osc,     ':',  'LineWidth',1.6, 'Color',[0.20 0.6 0.2]);
xlabel('time  [orbits]');  ylabel('V-bar  x  [m]');
title('Tangential impulse: x(t) = bounded hump + secular drift', ...
      'FontWeight','bold');
legend({'x(t) total','secular term  -3\DeltaV_x t', ...
        'oscillatory term  (4\DeltaV_x/n) sin nt'}, 'Location','southwest');

% --- third figure: x(t) and z(t) closed-form solutions vs time -----------
% Plots the two formulas from the note's Eq. (tang-impulse) directly, with
% the propagated solution overlaid as markers to confirm they coincide.
figure('Color','w'); hold on; grid on; box on;
plot(t/p.T, x_formula, 'LineWidth',2.2, 'Color',[0 0.45 0.74]);
plot(t/p.T, z_formula, 'LineWidth',2.2, 'Color',[0.85 0.1 0.1]);
% propagator values, sparsely sampled, as open markers on top
idx = 1:30:numel(t);
plot(t(idx)/p.T, x(idx), 'o', 'Color',[0 0.45 0.74], ...
     'MarkerSize',5, 'LineWidth',1.0);
plot(t(idx)/p.T, z(idx), 's', 'Color',[0.85 0.1 0.1], ...
     'MarkerSize',5, 'LineWidth',1.0);
xlabel('time  [orbits]');  ylabel('relative position  [m]');
title('Closed-form solution  x(t), z(t)  vs.  cw\_propagate3d', ...
      'FontWeight','bold');
legend({'x(t) = (\DeltaV_x/n)(4 sin nt - 3 n t)', ...
        'z(t) = (2\DeltaV_x/n)(cos nt - 1)', ...
        'x(t) propagated','z(t) propagated'}, 'Location','southwest');

% --- console summary ------------------------------------------------------
dx_orbit = -6*pi*dVx/p.n;                   % along-track shift per orbit
err_x = max(abs(x - x_formula));
err_z = max(abs(z - z_formula));
fprintf('impulse dVx        = %.3f m/s\n', dVx);
fprintf('secular drift rate = %.4f m/s   (-3 dVx)\n', -3*dVx);
fprintf('along-track shift  = %.1f m per orbit   (-6 pi dVx / n)\n', dx_orbit);
fprintf('after 1 orbit: (x,z) = (%.1f, %.1f) m\n', x(iT), z(iT));
fprintf('after 2 orbits: (x,z) = (%.1f, %.1f) m\n', x(end), z(end));
fprintf('\nclosed-form vs propagator agreement: max|dx| = %.2e m, ', err_x);
fprintf('max|dz| = %.2e m\n', err_z);
fprintf('\nThe chaser does NOT return to the target -- it drifts away each\n');
fprintf('orbit. This secular drift is why a pure tangential impulse is\n');
fprintf('passively unsafe inside the keep-out zone.\n');