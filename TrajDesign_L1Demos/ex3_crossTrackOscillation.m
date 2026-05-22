clear; clc; close all;
addpath('Dynamics')
p = cw_params(400);                         % 400 km circular target orbit
fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
        p.alt/1e3, p.n, p.T/60);

% starting point 
x0 = 0; y0 = 1;  z0 = 0;
y0dot = 1;
z0dot = 0;

x0dot = 0;


% % walkibg ellipse 
% x0dot =1.6*p.n*z0;
% 

s0 = [x0, y0,z0, x0dot, y0dot, z0dot];

 % two orbital periods
t = linspace(0, 10*p.T, 400);

% 6-by-numel(t) state history 
S = cw_propagate3d(s0, t, p.n);            

% pull out the position components for plotting
x = S(1,:);                                 % V-bar
y = S(2,:);                                 % R-bar
z = S(3,:);                                 % H-bar
 
% --- plot -----------------------------------------------------------------
figure; hold on; 
plot3(x,y, z, 'LineWidth', 2.2, 'Color', [0 0.45 0.74]);
plot3(x(1), y(1),  z(1),   'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
plot3(x(end), y(end), z(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
legend({'trajectory', 'start', 'end'}, 'Location', 'best');
 view(3)
fprintf('start (x,y,z) = (%.1f, %.1f, %.1f) m\n', x(1),   y(1),   z(1));
fprintf('end   (x,y,z) = (%.1f, %.1f, %.1f) m\n', x(end), y(end), z(end));
fprintf('drift rate (3/2) n z0 = %.4f m/s  ->  %.1f m over 2 orbits\n', ...
        1.5*p.n*z0, x(end)-x(1));


figure; hold on;
plot(t,y)