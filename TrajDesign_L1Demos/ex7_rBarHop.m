close all;
p = cw_params(400);                         % 400 km circular target orbit
fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
    p.alt/1e3, p.n, p.T/60);

% =========================================================================
%  R-BAR HOPPING
%  The chaser descends along R-bar from z = L to the target (z = 0) in a
%  sequence of hops between way-points (0, z_k).
%
%  Posing one hop  (0,z_a) -> (0,z_b):  three unknowns -- the impulse
%  (dVx,dVz) and the transfer time t_f -- against two position conditions
%  x(t_f)=0, z(t_f)=z_b.  The third condition is CLEAN ARRIVAL, zd(t_f)=0,
%  so the chaser reaches the way-point with zero radial velocity.  Writing
%  psi = n t_f, the two linear conditions give the impulse
%
%    dVx = n[(4 z_a - z_b) cos psi - 3 z_a] / [2(cos psi - 1)]
%    dVz = n (z_a - z_b) sin psi / (cos psi - 1)
%
%  and x(t_f)=0 reduces to a transcendental equation for psi:
%
%    3 psi (z_a - z_b cos psi) = 4 (z_a - z_b) sin psi
%
%  solved numerically for the root psi in (0,pi).  The chaser arrives with
%  a residual along-track velocity
%
%    xd(t_f) = n[z_a + (3 cos psi - 4) z_b] / [2(cos psi - 1)]
%
%  which a second impulse dVx2 = -xd(t_f) nulls.  Hop cost:
%    dV_hop = sqrt(dVx^2 + dVz^2) + |xd(t_f)|.
%
%  R-bar hopping is NOT passively safe: the way-points are not equilibria,
%  so a missed nulling impulse leaves a residual drift.
% =========================================================================

% ---- USER SETTINGS ------------------------------------------------------
L      = 3000;          % total R-bar distance to descend [m]
N_hops = 5;             % number of hops
g_hop  = 0.75;          % step schedule: g = 1 -> equal hops,
%                g < 1 -> hops shrink toward target
% -------------------------------------------------------------------------

% --- build the R-bar way-point schedule ----------------------------------
if abs(g_hop - 1) < 1e-9
    steps = repmat(L/N_hops, 1, N_hops);            % equal spacing
else
    d1    = L*(1-g_hop)/(1-g_hop^N_hops);           % first (largest) hop
    steps = d1 * g_hop.^(0:N_hops-1);               % geometric grading
end
zwp = [L, L - cumsum(steps)];                       % way-points z_0..z_N

% =========================================================================
%  FIGURE 1 -- single R-bar hop, in detail
% =========================================================================
za1 = zwp(1);  zb1 = zwp(2);
[dVx1,dVz1,tf1,xdr1] = rbar_hop_solve(za1, zb1, p.n);
S1 = cw_propagate3d([0,0,za1,dVx1,0,dVz1], linspace(0,tf1,300), p.n);

figure; hold on; grid on; box on;
plot(S1(1,:), S1(3,:), 'LineWidth',2.6, 'Color',[0.85 0.1 0.1]);
plot(0, za1, 'go','MarkerFaceColor','g','MarkerSize',8);
plot(0, zb1, 'ro','MarkerFaceColor','r','MarkerSize',8);
plot(0, 0,  'ks','MarkerSize',12,'MarkerFaceColor','none','LineWidth',1.6);
text(0,za1,'  z_a','FontSize',9,'FontWeight','bold');
text(0,zb1,'  z_b','FontSize',9,'FontWeight','bold');
text(0,0,  '  target','FontSize',8,'VerticalAlignment','bottom');
set(gca,'YDir','reverse'); axis equal;
xlabel('V-bar   x   [m]   (+ = ahead)');
ylabel('R-bar   z   [m]   (+ = toward Earth)');
title(sprintf(['Single R-bar hop  (0,%.0f) \\rightarrow (0,%.0f):  ' ...
    't_f = %.3f T'], za1, zb1, tf1/p.T), 'FontWeight','bold');
legend({'hop trajectory (clean arrival)','z_a (start)','z_b (arrival)'}, ...
    'Location','best');

% =========================================================================
%  FIGURE 2 -- the full multi-hop R-bar descent
% =========================================================================
figure; hold on; grid on; box on;
dV_total = 0;
hop_data = zeros(N_hops,4);             % [za, zb, tf, dV_hop]
for k = 1:N_hops
    za = zwp(k);  zb = zwp(k+1);
    [dVx,dVz,tf,xdr] = rbar_hop_solve(za, zb, p.n);
    S = cw_propagate3d([0,0,za,dVx,0,dVz], linspace(0,tf,250), p.n);
    plot(S(1,:), S(3,:), 'LineWidth',2.2, 'Color',[0.85 0.1 0.1]);
    plot(0, za, 'ko','MarkerFaceColor','k','MarkerSize',5);
    dV_hop   = hypot(dVx,dVz) + abs(xdr);
    dV_total = dV_total + dV_hop;
    hop_data(k,:) = [za, zb, tf, dV_hop];
end
plot(0, 0, 'ks','MarkerSize',12,'MarkerFaceColor','none','LineWidth',1.6);
set(gca,'YDir','reverse'); axis equal;
xlabel('V-bar   x   [m]   (+ = ahead)');
ylabel('R-bar   z   [m]   (+ = toward Earth)');
if abs(g_hop-1) < 1e-9
    sched = 'equal hops';
else
    sched = sprintf('geometric grading, g = %.2f', g_hop);
end
title(sprintf('R-bar hopping descent:  %d hops over L = %.0f m  (%s)', ...
    N_hops, L, sched), 'FontWeight','bold');

% =========================================================================
%  FIGURE 3 -- transfer phase and hop cost vs. way-point
% =========================================================================
figure;
subplot(2,1,1); hold on; grid on; box on;
plot(1:N_hops, hop_data(:,3)/p.T, 'o-','LineWidth',2, ...
    'Color',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74]);
ylabel('transfer time  t_f / T');
title('R-bar hop: transfer time and cost per hop','FontWeight','bold');
xlim([0.5 N_hops+0.5]);

subplot(2,1,2); hold on; grid on; box on;
bar(1:N_hops, hop_data(:,4)*100, 0.6, 'FaceColor',[0.85 0.30 0.10]);
xlabel('hop number  (1 = first, farthest from target)');
ylabel('hop \DeltaV  [cm/s]');
xlim([0.5 N_hops+0.5]);

% =========================================================================
%  console summary
% =========================================================================
fprintf('R-bar descent: L = %.0f m in %d hops  (%s)\n\n', L, N_hops, sched);
fprintf(' hop    z_a     z_b     t_f/T    dV [cm/s]\n');
fprintf(' ---  ------  ------  -------  ----------\n');
for k = 1:N_hops
    fprintf(' %2d   %6.0f  %6.0f   %.4f   %8.2f\n', k, ...
        hop_data(k,1), hop_data(k,2), hop_data(k,3)/p.T, ...
        hop_data(k,4)*100);
end
fprintf(' ---  ------  ------  -------  ----------\n');
fprintf(' total R-bar hopping dV = %.4f m/s\n', dV_total);
fprintf('\nUnlike V-bar hopping (total dV = nL/2, independent of N), the\n');
fprintf('R-bar total depends on the way-point schedule and is obtained\n');
fprintf('numerically -- each hop solves its own transcendental phase eqn.\n');
fprintf('\nR-bar hopping is NOT passively safe: a missed nulling impulse\n');
fprintf('leaves the residual along-track velocity xd(t_f) uncancelled,\n');
fprintf('and the chaser drifts off with no closed orbit to return on.\n');

 % demo_rbar_hop

% =========================================================================
%  local function: solve one clean-arrival R-bar hop
% =========================================================================
function [dVx, dVz, tf, xd_res] = rbar_hop_solve(za, zb, n)
%RBAR_HOP_SOLVE  Clean-arrival R-bar hop (0,za) -> (0,zb).
%   Solves the transcendental phase equation
%     3 psi (za - zb cos psi) = 4 (za - zb) sin psi
%   for psi in (0,pi), then evaluates the closed-form impulse and the
%   residual along-track velocity at arrival.

% transcendental phase equation, root in (0, pi)
f = @(psi) 3*psi.*(za - zb*cos(psi)) - 4*(za - zb)*sin(psi);
lo = 1e-4;  hi = pi - 1e-4;
if f(lo)*f(hi) > 0
    % scan for a sign change if the default bracket fails
    pg = linspace(lo, hi, 400);
    fg = arrayfun(f, pg);
    idx = find(fg(1:end-1).*fg(2:end) < 0, 1);
    if isempty(idx)
        error('rbar_hop_solve:noroot', ...
            'No phase root for hop (%.0f -> %.0f).', za, zb);
    end
    lo = pg(idx);  hi = pg(idx+1);
end
for it = 1:100                              % bisection
    m = 0.5*(lo+hi);
    if f(lo)*f(m) <= 0,  hi = m;  else,  lo = m;  end
end
psi = 0.5*(lo+hi);

% closed-form impulse and residual along-track velocity
c = cos(psi);
dVx    = n*((4*za - zb)*c - 3*za) / (2*(c - 1));
dVz    = n*(za - zb)*sin(psi)     / (c - 1);
xd_res = n*(za + (3*c - 4)*zb)    / (2*(c - 1));
tf     = psi / n;
end
