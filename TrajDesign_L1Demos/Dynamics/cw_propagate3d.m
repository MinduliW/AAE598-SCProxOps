function S = cw_propagate3d(s0, t, n)
%CW_PROPAGATE3D  Analytic 3D Clohessy-Wiltshire state propagation.
%
%   S = CW_PROPAGATE3D(s0, t, n)
%
%   Propagates the full six-element relative state under the unforced
%   Clohessy-Wiltshire dynamics, using the closed-form analytic solution
%   (no numerical integration).
%
%   Inputs:
%     s0 = [x; y; z; xd; yd; zd]   initial relative state
%            x  = V-bar position  (along-track)              [m]
%            y  = H-bar position  (cross-track)              [m]
%            z  = R-bar position  (radial, +z toward Earth)  [m]
%            xd, yd, zd = the corresponding velocities       [m/s]
%     t  = scalar or vector of times [s]
%     n  = mean motion [rad/s]
%
%   Output:
%     S  = 6-by-numel(t) array; rows are [x; y; z; xd; yd; zd].
%
%   Convention (matching the AE 598 note):
%     x = V-bar (along-track, +x along the target velocity)
%     y = H-bar (cross-track)
%     z = R-bar (radial, +z TOWARD Earth -> a lower, faster orbit)
%
%   The in-plane (x, z) motion is coupled through the Coriolis and
%   gravity-gradient terms; the out-of-plane (y) motion is a decoupled
%   simple harmonic oscillator at the orbital rate. Both are propagated
%   here by their exact closed-form solutions.
%
%   Governing equations:
%     xdd - 2 n zd          = 0
%     zdd + 2 n xd - 3 n^2 z = 0
%     ydd + n^2 y           = 0
%
%   Example:
%     n  = 1.131e-3;                       % 400 km orbit
%     s0 = [0; 200; 500; 0; 0.3; 0];       % [x y z xd yd zd]
%     S  = cw_propagate3d(s0, linspace(0,6000,300), n);
%     plot3(S(1,:), S(2,:), S(3,:));       % V-bar, H-bar, R-bar

    % --- unpack initial state ([x y z xd yd zd] ordering) --------------
    x0  = s0(1);   y0  = s0(2);   z0  = s0(3);
    xd0 = s0(4);   yd0 = s0(5);   zd0 = s0(6);

    t  = t(:).';                       % force row vector
    nt = n * t;
    s  = sin(nt);   c = cos(nt);

    % --- in-plane (V-bar / R-bar) closed-form solution -----------------
    x  = x0 + (6*n*z0 - 3*xd0).*t ...
            + (4*xd0 - 6*n*z0)/n.*s ...
            + 2*zd0/n.*(1 - c);
    z  = (4*z0 - 2*xd0/n) ...
            - (3*z0 - 2*xd0/n).*c ...
            + zd0/n.*s;
    xd = (6*n*z0 - 3*xd0) ...
            + (4*xd0 - 6*n*z0).*c ...
            + 2*zd0.*s;
    zd = (3*n*z0 - 2*xd0).*s + zd0.*c;

    % --- out-of-plane (H-bar) decoupled harmonic solution --------------
    y  =  y0.*c + (yd0/n).*s;
    yd = -y0*n.*s + yd0.*c;

    % --- assemble [x; y; z; xd; yd; zd] --------------------------------
    S = [x; y; z; xd; yd; zd];
end