clear; clc; close all;
p = cw_params(400);                         % 400 km circular target orbit
fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
        p.alt/1e3, p.n, p.T/60);

% --- radial impulse at the target -----------------------------------------
% A pure +dVz on R-bar, applied AT the target itself.  Convention: +dVz is a
% radially OUTWARD impulse (away from Earth); in the lecture's +z-toward-Earth
% axis that is a post-impulse velocity  zd0 = -dVz.  The post-impulse state
% is therefore  x0 = y0 = z0 = 0,  xd0 = yd0 = 0,  zd0 = -dVz.
% Propagating that state with the CW solution reproduces the note's
%   x(t) =  (2 dVz/n)(cos nt - 1)
%   z(t) = -(  dVz/n) sin nt
dVz = 0.15;                                 % radial impulse [m/s]  (+ = outward)
s0  = [0, 0, 0, 0, 0, -dVz];                % [x y z xd yd zd]

% --- propagate ------------------------------------------------------------
t = linspace(0, 2*p.T, 600);                % two orbital periods
S = cw_propagate3d(s0, t, p.n);

x = S(1,:);                                 % V-bar
y = S(2,:);                                 % H-bar
z = S(3,:);                                 % R-bar

% --- closed-form solution from the note's Eq. (rad-impulse) ---------------
% x(t) = (2 dVz/n)(cos nt - 1),   z(t) = -(dVz/n) sin nt
x_formula = s0(1) +  (2*dVz/p.n) * (cos(p.n*t) - 1);
z_formula = s0(3) -  (  dVz/p.n) * sin(p.n*t);

% --- plot: trajectory + closed-form overlay -------------------------------
figure; hold on; grid on; box on;
plot(x, z, 'LineWidth', 3.0, 'Color', [0.75 0.90 0.78]);
% closed-form curve plotted on top with a dashed line
plot(x_formula, z_formula, '--', 'LineWidth',1.6, 'Color',[0.85 0.1 0.1]);
plot(x(1),   z(1),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
plot(x(end), z(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);
plot(0, 0, 'ks', 'MarkerSize',12, 'MarkerFaceColor','none', 'LineWidth',1.6);
text(0, 0, '  target / impulse point', 'FontSize',8, ...
     'VerticalAlignment','bottom');
set(gca,'YDir','reverse');                  % +z (R-bar, toward Earth) down
axis equal;
xlabel('V-bar   x   [m]   (+ = ahead)');
ylabel('R-bar   z   [m]   (+ = toward Earth)');
title(sprintf('Radial impulse  \\DeltaV_z = %.2f m/s', dVz), ...
      'FontWeight','bold');
legend({'propagated (cw\_propagate3d)','closed-form Eq.', ...
        'impulse point','after 2 orbits'}, 'Location','best');

% mark the half-orbit point: the far end of the ellipse, dx = -4 dVz/n
[~,iH] = min(abs(t - p.T/2));
plot(x(iH), z(iH), 'ks', 'MarkerFaceColor','y', 'MarkerSize',8);
text(x(iH), z(iH), '  half orbit (\Deltax = -4\DeltaV_z/n)', 'FontSize',8);

% --- second figure: x(t) and z(t) closed-form solutions vs time ----------
% Plots the two formulas from the note's Eq. (rad-impulse) directly, with the
% propagated solution overlaid as markers to confirm they coincide.
figure('Color','w'); hold on; grid on; box on;
plot(t/p.T, x_formula, 'LineWidth',2.2, 'Color',[0 0.45 0.74]);
plot(t/p.T, z_formula, 'LineWidth',2.2, 'Color',[0.85 0.1 0.1]);
idx = 1:30:numel(t);
plot(t(idx)/p.T, x(idx), 'o', 'Color',[0 0.45 0.74], ...
     'MarkerSize',5, 'LineWidth',1.0);
plot(t(idx)/p.T, z(idx), 's', 'Color',[0.85 0.1 0.1], ...
     'MarkerSize',5, 'LineWidth',1.0);
xlabel('time  [orbits]');  ylabel('relative position  [m]');
title('Closed-form solution  x(t), z(t)  vs.  cw\_propagate3d', ...
      'FontWeight','bold');
legend({'x(t) = (2\DeltaV_z/n)(cos nt - 1)', ...
        'z(t) = -(\DeltaV_z/n) sin nt', ...
        'x(t) propagated','z(t) propagated'}, 'Location','best');

% --- console summary ------------------------------------------------------
dx_half = -4*dVz/p.n;                       % along-track shift after half orbit
Ar = dVz/p.n;                               % radial semi-axis
err_x = max(abs(x - x_formula));
err_z = max(abs(z - z_formula));
fprintf('impulse dVz        = %.3f m/s   (+ = radially outward)\n', dVz);
fprintf('ellipse semi-axes  = %.1f m (radial) x %.1f m (along-track)  [2:1]\n', ...
        Ar, 2*Ar);
fprintf('half-orbit shift   = %.1f m   (-4 dVz / n)\n', dx_half);
fprintf('at every integer orbit the chaser is exactly back at (0, 0):\n');
fprintf('  x(kT) = (2 dVz/n)(cos 2pi k - 1) = 0,   z(kT) = -(dVz/n) sin 2pi k = 0\n');
fprintf('\nclosed-form vs propagator agreement: max|dx| = %.2e m, ', err_x);
fprintf('max|dz| = %.2e m\n', err_z);
fprintf('\nThe trajectory is a CLOSED 2:1 ellipse -- the chaser returns to\n');
fprintf('the impulse point every orbit, with no secular drift. If a second\n');
fprintf('impulse is missed, the chaser simply coasts back to the start:\n');
fprintf('this is why the radial impulse is passively safe.\n');