function demo_worked_docking()

addpath('Dynamics/')
%DEMO_WORKED_DOCKING  The worked -V-bar docking approach (note Sec. 6).
%
%   Self-contained simulation of the complete five-element approach strategy
%   from the lecture note's worked example: docking on the -V-bar port of an
%   ISS-class target at 400 km.  The chaser is taken from the end of phasing
%   (S0) through to contact (S5):
%
%     S0 -> S1   free drift on a lower orbit      (RGPS filter convergence)
%     S1 -> S2   two-impulse Hohmann transfer     (min-dV altitude change)
%     S2 -> S3   V-bar hold / short drift         (timeline synchronisation)
%     S3 -> S4   radial-boost transfer            (passive safety in the AE)
%     S4 -> S5   forced straight-line approach    (corridor-aligned docking)
%
%   Each leg is built here directly with cw_propagate3d -- the same analytic
%   CW propagator used by the other demos.  An impulse is applied simply by
%   adding the delta-V to the velocity components of the leg's start state.
%
%   Depends only on  cw_params  and  cw_propagate3d.

    close all;
    p = cw_params(400);
    n = p.n;  T = p.T;
    fprintf('Target altitude %.0f km   n = %.4e rad/s   T = %.1f min\n\n', ...
            p.alt/1e3, n, T/60);

    % =====================================================================
    %  BUILD THE TRAJECTORY, LEG BY LEG, WITH cw_propagate3d
    %  State ordering for cw_propagate3d:  [x y z xd yd zd]
    %    x = V-bar,  y = H-bar,  z = R-bar (+z toward Earth)
    % =====================================================================
    z_low = 4000;                       % chaser starts 4 km below the target
    x_S0  = -26000;                     % ... and 26 km behind it

    % ---- leg A : S0 -> S1  free drift on the lower orbit ----------------
    % On a circular lower orbit the chaser drifts forward at xd = (3/2) n z.
    xd_drift = 1.5 * n * z_low;
    dx_hoh   = 0.75 * pi * z_low;        % along-track advance of the Hohmann
    x_S2     = -3000;                    % Hohmann ends here, outside the AE
    x_S1     = x_S2 - dx_hoh;
    t_drift  = (x_S1 - x_S0) / xd_drift;
    tA  = linspace(0, t_drift, 300);
    SA  = cw_propagate3d([x_S0, 0, z_low, xd_drift, 0, 0], tA, n);
    legA = [SA(1,:); SA(3,:)];
    endA = SA(:,end);                    % full state at S1

    % ---- leg B : S1 -> S2  two-impulse Hohmann transfer -----------------
    % First impulse: a tangential burn dVx = n*dz/4 added to the V-bar vel.
    dVx = n * z_low / 4;
    sB0 = endA + [0;0;0; dVx;0;0];       % apply the first Hohmann burn
    tB  = linspace(0, T/2, 300);
    SB  = cw_propagate3d(sB0.', tB, n);
    legB = [SB(1,:); SB(3,:)];
    endB = SB(:,end);                    % state at S2 (2nd burn circularises)

    % ---- leg C : S2 -> S3  V-bar hold / short drift ---------------------
    x_S3 = -2700;
    legC = [linspace(endB(1), x_S3, 50); zeros(1,50)];

    % ---- leg D : S3 -> S4  radial-boost transfer ------------------------
    % A radial impulse from rest on V-bar: net along-track shift after a
    % half orbit is  dx = 4 zd0 / n,  so zd0 = (x_S4 - x_S3) n / 4.
    x_S4 = -500;
    zd0  = (x_S4 - x_S3) * n / 4;
    tD   = linspace(0, T/2, 300);
    SD   = cw_propagate3d([x_S3, 0, 0, 0, 0, zd0], tD, n);
    legD = [SD(1,:); SD(3,:)];

    % ---- leg E : S4 -> S5  forced straight-line approach ----------------
    legE = [linspace(x_S4, 0, 50); zeros(1,50)];

    legs  = {legA, legB, legC, legD, legE};
    names = {'S0->S1 free drift','S1->S2 Hohmann', ...
             'S2->S3 V-bar hold','S3->S4 radial boost', ...
             'S4->S5 forced approach'};
    nodes = struct('S0',[x_S0 z_low], 'S1',[endA(1) endA(3)], ...
                   'S2',[endB(1) endB(3)], 'S3',[x_S3 0], ...
                   'S4',[x_S4 0], 'S5',[0 0]);

    % delta-V budget
    dV_hohmann     = 2 * dVx;            % two equal tangential burns
    dV_radialboost = 2 * abs(zd0);       % insertion + nulling radial burns
    dV_total       = dV_hohmann + dV_radialboost;

    legCols = { [.45 .45 .45], [0 .45 .74], [.60 .50 .10], ...
                [.20 .60 .20], [.75 .20 .65] };
    nodeNames = fieldnames(nodes);
    th = linspace(0, 2*pi, 200);

    figure('Color','w','Position',[80 80 980 760]);

    % =====================================================================
    %  PANEL (a) -- OVERVIEW : S0 through to contact
    % =====================================================================
    ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    for k = 1:numel(legs)
        plot(ax1, legs{k}(1,:), legs{k}(2,:), ...
             'LineWidth',2.6, 'Color',legCols{k});
    end
    for k = 1:numel(nodeNames)
        P = nodes.(nodeNames{k});
        plot(ax1, P(1), P(2), 'ko','MarkerFaceColor','k','MarkerSize',6);
        text(ax1, P(1), P(2), ['  ' nodeNames{k}], ...
             'FontWeight','bold','FontSize',9,'VerticalAlignment','bottom');
    end
    plot(ax1, 0,0,'ks','MarkerSize',12,'MarkerFaceColor','none','LineWidth',1.6);
    plot(ax1, 2000*cos(th), 1000*sin(th), '--', ...
         'Color',[.85 .1 .1],'LineWidth',1.2);
    text(ax1, 1500, -850, 'AE', 'Color',[.85 .1 .1],'FontWeight','bold');
    set(ax1,'YDir','reverse'); axis(ax1,'equal');
    xlabel(ax1,'V-bar   x   [m]   (+ = ahead)');
    ylabel(ax1,'R-bar   z   [m]   (+ = toward Earth)');
    title(ax1,'(a)  Overview:  S_0 (end of phasing)  \rightarrow  contact', ...
          'FontWeight','bold');

    % =====================================================================
    %  PANEL (b) -- CLOSE-RANGE ZOOM : the AE interior, S2 -> S5
    % =====================================================================
    ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
    for k = 3:numel(legs)
        plot(ax2, legs{k}(1,:), legs{k}(2,:), ...
             'LineWidth',2.8, 'Color',legCols{k});
    end
    for nm = {'S2','S3','S4','S5'}
        P = nodes.(nm{1});
        plot(ax2, P(1), P(2), 'ko','MarkerFaceColor','k','MarkerSize',6);
        text(ax2, P(1), P(2), ['  ' nm{1}], ...
             'FontWeight','bold','FontSize',9,'VerticalAlignment','bottom');
    end
    plot(ax2, 0,0,'ks','MarkerSize',12,'MarkerFaceColor','none','LineWidth',1.6);
    plot(ax2, 2000*cos(th), 1000*sin(th), '--', ...
         'Color',[.85 .1 .1],'LineWidth',1.0);
    plot(ax2, 200*cos(th), 200*sin(th), '-', ...
         'Color',[.90 .55 .10],'LineWidth',1.4);
    text(ax2, 250, -300, 'KOZ', 'Color',[.90 .55 .10],'FontWeight','bold');
    half = deg2rad(12);
    plot(ax2, [-500 0], [ tan(half)*500 0], ':', 'Color',[.5 .5 .5]);
    plot(ax2, [-500 0], [-tan(half)*500 0], ':', 'Color',[.5 .5 .5]);
    set(ax2,'YDir','reverse'); axis(ax2,'equal');
    xlim(ax2,[-3200 600]);  ylim(ax2,[-900 900]);
    xlabel(ax2,'V-bar   x   [m]   (+ = ahead)');
    ylabel(ax2,'R-bar   z   [m]   (+ = toward Earth)');
    title(ax2,'(b)  Close-range zoom:  hold, radial boost, forced approach', ...
          'FontWeight','bold');
    legend(ax2, {'S_2\rightarrowS_3 V-bar hold', ...
                 'S_3\rightarrowS_4 radial boost', ...
                 'S_4\rightarrowS_5 forced approach'}, 'Location','northwest');

    % =====================================================================
    %  console summary
    % =====================================================================
    drivers = {'RGPS filter convergence', ...
               'minimum-dV altitude change', ...
               'timeline synchronisation', ...
               'passive safety inside the AE', ...
               'corridor-aligned final approach'};
    fprintf('Worked -V-bar docking approach -- element breakdown:\n');
    fprintf('  %-26s  %-22s\n', 'element', 'driver');
    fprintf('  %s\n', repmat('-',1,52));
    for k = 1:numel(names)
        fprintf('  %-26s  %-22s\n', names{k}, drivers{k});
    end
    fprintf('\n  node positions (V-bar, R-bar) [m]:\n');
    for k = 1:numel(nodeNames)
        P = nodes.(nodeNames{k});
        fprintf('    %-3s = (%8.0f, %7.0f)\n', nodeNames{k}, P(1), P(2));
    end
    fprintf('\n  delta-V budget:\n');
    fprintf('    Hohmann transfer (S1->S2)      = %.4f m/s\n', dV_hohmann);
    fprintf('    radial-boost transfer (S3->S4) = %.4f m/s\n', dV_radialboost);
    fprintf('    %s\n', repmat('-',1,42));
    fprintf('    total                          = %.4f m/s\n', dV_total);
    fprintf('\n  The free drift (S0->S1), V-bar hold (S2->S3), and the\n');
    fprintf('  passive coast inside each transfer cost no dV; the forced\n');
    fprintf('  straight-line approach (S4->S5) is closed-loop controlled.\n');
end