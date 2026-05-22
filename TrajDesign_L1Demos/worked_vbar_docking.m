function out = worked_vbar_docking(p)
%WORKED_VBAR_DOCKING  The worked -V-bar docking approach (note Sec. 6).
%
%   out = WORKED_VBAR_DOCKING(p)
%
%   Assembles the complete five-element approach strategy from the note's
%   worked example: a docking on the -V-bar port of an ISS-class target.
%
%       S0 -> S1   free drift on a lower orbit       (RGPS convergence)
%       S1 -> S2   two-impulse Hohmann transfer      (min-dV altitude change)
%       S2 -> S3   V-bar hold / drift                (timeline sync)
%       S3 -> S4   radial-boost transfer             (passive safety in AE)
%       S4 -> S5   forced straight-line approach     (corridor alignment)
%
%   p   : parameter struct from CW_PARAMS.
%   out : struct with the leg-by-leg trajectories, node positions, and
%         the dV accounting.

    n = p.n;  T = p.T;

    % --- S0 -> S1 : free drift on a lower orbit -------------------------
    z_low  = 4000;                       % start 4 km below the target
    x_S0   = -26000;                     % 26 km behind the target
    xd_d   = 1.5*n*z_low;                % secular drift velocity
    dx_hoh = 0.75*pi*z_low;              % Hohmann along-track advance
    x_S2   = -3000;                      % Hohmann ends here, outside the AE
    x_S1   = x_S2 - dx_hoh;
    t_d    = (x_S1 - x_S0)/xd_d;
    S = cw_propagate([x_S0; z_low; xd_d; 0], linspace(0,t_d,300), n);
    legA = [S(1,:); S(2,:)];
    S1   = [S(1,end), S(2,end)];

    % --- S1 -> S2 : Hohmann transfer up to the target orbit ------------
    dVx = n*z_low/4;
    seg = [T/2, dVx, 0;  1.0, dVx, 0];
    [trajH,~] = cw_chain([S1(1); z_low; xd_d; 0], seg, n);
    legB = trajH(2:3,:);
    S2   = [trajH(2,end), 0];

    % --- S2 -> S3 : V-bar hold / short drift to the AE boundary ---------
    S3   = [-2700, 0];
    legC = [linspace(S2(1),S3(1),50); zeros(1,50)];

    % --- S3 -> S4 : passively safe radial-boost transfer ---------------
    x_S4 = -500;                         % corridor entry behind the port
    zd_b = n*(x_S4 - S3(1))/4;           % net dx = +4 zd/n  -> zd > 0
    segB = [T/2, 0, +zd_b;  1.0, 0, -zd_b];
    [trajB,~] = cw_chain([S3(1);0;0;0], segB, n);
    legD = trajB(2:3,:);
    S4   = [trajB(2,end), 0];

    % --- S4 -> S5 : forced straight-line approach to the port ----------
    legE = [linspace(S4(1),0,50); zeros(1,50)];

    % --- collect --------------------------------------------------------
    out.legs  = {legA, legB, legC, legD, legE};
    out.names = {'S0->S1 free drift','S1->S2 Hohmann', ...
                 'S2->S3 V-bar hold','S3->S4 radial boost', ...
                 'S4->S5 forced approach'};
    out.nodes = struct('S0',[x_S0 z_low],'S1',S1,'S2',S2, ...
                       'S3',S3,'S4',S4,'S5',[0 0]);
    out.dV.hohmann     = 2*dVx;
    out.dV.radialboost = 2*abs(zd_b);
    out.dV.total       = out.dV.hohmann + out.dV.radialboost;
    out.note = sprintf(['worked -V-bar docking: Hohmann %.3f + ' ...
                        'radial-boost %.3f = %.3f m/s'], ...
                        out.dV.hohmann, out.dV.radialboost, out.dV.total);
end
