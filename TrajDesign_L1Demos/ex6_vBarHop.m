clear; clc; close all;
p = cw_params(400);                         % 400 km circular target orbit
fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
        p.alt/1e3, p.n, p.T/60);

% =========================================================================
%  MULTI-HOP V-BAR HOPPING
%  The chaser descends along V-bar from x = L to the target (x = 0).
%
%  Two independent design choices are illustrated here:
%
%   (1) HOP TYPE -- classical vs advanced.
%       Classical: pure radial impulse, transfer time fixed at half an
%         orbit (psi = n t_f = pi).  Advanced (Sasaki et al. 2021): the
%         transfer time is a free parameter and the impulse has both a
%         radial and a tangential component.  With psi = n t_f and step dx:
%           dVx = n sin(psi) / D(psi) * dx
%           dVz = 4 n sin^2(psi/2) / D(psi) * dx
%           D(psi) = 8(1 - cos psi) - 3 psi sin psi        (PA-safe: D > 0)
%
%   (2) STEP SCHEDULE -- equal vs geometric.
%       Equal:     every hop has the same length dx = -L/N.
%       Geometric: each hop is a fraction g of the previous one,
%           dx_{k+1} = g * dx_k,  g < 1.  Since the radial excursion of a
%           hop scales with its length, the hops get shorter AND flatter
%           as the chaser nears the target -- matching the improving
%           navigation accuracy at close range.  This is the gamma_hop
%           schedule of the advanced-hopping literature.
%
%  Note: hop type and step schedule are SEPARATE choices.  "Advanced" by
%  itself does not change the hop height -- equal-length advanced hops are
%  all the same size.  It is the geometric schedule that shrinks them.
% =========================================================================

% ---- USER SETTINGS ------------------------------------------------------
L        = 2000;        % total V-bar distance to close [m]
N_hops   = 5;           % number of hops
tof_frac = 0.30;        % advanced-hop time of flight per hop, fraction of T
                        %   (0 < tof_frac < 0.5)
g_hop    = 0.70;        % geometric step ratio  (g < 1 -> shrinking hops)
% -------------------------------------------------------------------------

psi_a = 2*pi * tof_frac;                    % advanced-hop transfer phase
if 8*(1-cos(psi_a)) - 3*psi_a*sin(psi_a) <= 0
    warning('tof_frac too large: psi outside the PA-safe range (0, pi).');
end

% per-hop impulse for a given step dx and transfer phase psi
hop_impulse = @(dx,psi) struct( ...
    'dVx', p.n*sin(psi)/(8*(1-cos(psi))-3*psi*sin(psi)) * dx, ...
    'dVz', 4*p.n*sin(psi/2)^2/(8*(1-cos(psi))-3*psi*sin(psi)) * dx, ...
    'tf' , psi/p.n );

% --- build the four hop schedules ----------------------------------------
% equal steps
dx_equal = repmat(-L/N_hops, 1, N_hops);
% geometric steps (first hop longest, shrinking by g each time)
dx1      = L*(1-g_hop)/(1-g_hop^N_hops);
dx_geom  = -dx1 * g_hop.^(0:N_hops-1);

% =========================================================================
%  helper: plot a multi-hop chain, return total dV and the hop list
% =========================================================================
function [dV,hops] = run_chain(p, steps, psi, color, ax)
    xk = sum(abs(steps));                   % start at x = L
    dV = 0;  hops = [];
    for k = 1:numel(steps)
        dx  = steps(k);
        D   = 8*(1-cos(psi)) - 3*psi*sin(psi);
        dVx = p.n*sin(psi)/D * dx;
        dVz = 4*p.n*sin(psi/2)^2/D * dx;
        tf  = psi/p.n;
        S = cw_propagate3d([xk,0,0,dVx,0,dVz], linspace(0,tf,200), p.n);
        plot(ax, S(1,:), S(3,:), 'LineWidth',2.2, 'Color',color);
        plot(ax, xk, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
        dV = dV + 2*hypot(dVx,dVz);          % insertion + nulling
        hops = [hops; abs(dx), max(abs(S(3,:)))];   %#ok<AGROW>
        xk = S(1,end);
    end
    plot(ax, xk, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
    plot(ax, 0,0,'ks','MarkerSize',12,'MarkerFaceColor','none','LineWidth',1.6);
end

% =========================================================================
%  FIGURE 1 -- classical equal-step hopping  (psi = pi)
% =========================================================================
f1 = figure; ax1 = axes(f1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
[dV_classical,~] = run_chain(p, dx_equal, pi, [0 0.45 0.74], ax1);
set(ax1,'YDir','reverse'); axis(ax1,'equal');
xlabel(ax1,'V-bar   x   [m]   (+ = ahead)');
ylabel(ax1,'R-bar   z   [m]   (+ = toward Earth)');
title(ax1,sprintf(['Classical equal-step hopping:  %d hops, L = %.0f m  ' ...
      '(\\psi = \\pi)'], N_hops, L), 'FontWeight','bold');

% =========================================================================
%  FIGURE 2 -- advanced equal-step hopping  (user tof, same size each)
% =========================================================================
f2 = figure; ax2 = axes(f2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
[dV_adv_equal,hops_eq] = run_chain(p, dx_equal, psi_a, [0.85 0.30 0.10], ax2);
set(ax2,'YDir','reverse'); axis(ax2,'equal');
xlabel(ax2,'V-bar   x   [m]   (+ = ahead)');
ylabel(ax2,'R-bar   z   [m]   (+ = toward Earth)');
title(ax2,sprintf(['Advanced equal-step hopping:  %d hops, t_f = %.2f T  ' ...
      '(all hops same size)'], N_hops, tof_frac), 'FontWeight','bold');

% =========================================================================
%  FIGURE 3 -- advanced GEOMETRIC hopping  (hops shrink toward the target)
% =========================================================================
f3 = figure; ax3 = axes(f3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
[dV_adv_geom,hops_gm] = run_chain(p, dx_geom, psi_a, [0.20 0.55 0.20], ax3);
set(ax3,'YDir','reverse'); axis(ax3,'equal');
xlabel(ax3,'V-bar   x   [m]   (+ = ahead)');
ylabel(ax3,'R-bar   z   [m]   (+ = toward Earth)');
title(ax3,sprintf(['Advanced geometric hopping:  %d hops, ratio g = %.2f  ' ...
      '(hops shrink toward target)'], N_hops, g_hop), 'FontWeight','bold');

% =========================================================================
%  console summary
% =========================================================================
fprintf('descent: L = %.0f m closed in %d hops\n\n', L, N_hops);

fprintf('classical equal-step (psi = pi):\n');
fprintf('  total dV = %.4f m/s   (= n L / 2 = %.4f),  time = %.2f T\n\n', ...
        dV_classical, p.n*L/2, N_hops*0.5);

fprintf('advanced equal-step (t_f = %.2f T):\n', tof_frac);
fprintf('  hop lengths : ');  fprintf('%.0f ', hops_eq(:,1));  fprintf('m\n');
fprintf('  hop heights : ');  fprintf('%.0f ', hops_eq(:,2));  fprintf('m\n');
fprintf('  total dV = %.4f m/s,  time = %.2f T\n\n', ...
        dV_adv_equal, N_hops*tof_frac);

fprintf('advanced geometric (t_f = %.2f T, g = %.2f):\n', tof_frac, g_hop);
fprintf('  hop lengths : ');  fprintf('%.0f ', hops_gm(:,1));  fprintf('m\n');
fprintf('  hop heights : ');  fprintf('%.0f ', hops_gm(:,2));  fprintf('m\n');
fprintf('  total dV = %.4f m/s,  time = %.2f T\n\n', ...
        dV_adv_geom, N_hops*tof_frac);

fprintf('The radial excursion of a hop scales with its length, so the\n');
fprintf('geometric schedule shrinks the hops in BOTH length and height\n');
fprintf('as the chaser nears the target -- equal-step hopping does not.\n');