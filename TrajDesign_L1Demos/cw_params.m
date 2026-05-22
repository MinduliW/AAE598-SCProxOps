function p = cw_params(altitude_km)
%CW_PARAMS  Orbital parameters for the Clohessy-Wiltshire model.
%
%   p = CW_PARAMS()             uses a 400 km circular target orbit.
%   p = CW_PARAMS(altitude_km)  uses the given target altitude [km].
%
%   Returns a struct with fields:
%     mu  - Earth gravitational parameter [m^3/s^2]
%     Re  - Earth mean radius             [m]
%     alt - target orbit altitude         [m]
%     a   - target orbit radius           [m]
%     n   - mean motion                   [rad/s]
%     T   - orbital period                [s]
%
%   Convention used throughout this suite (matching the AE 598 note):
%     x = V-bar  (along-track, +x along the target velocity)
%     y = H-bar  (cross-track)
%     z = R-bar  (radial, +z TOWARD Earth -> a lower, faster orbit)

    if nargin < 1 || isempty(altitude_km)
        altitude_km = 400;
    end
    p.mu  = 3.986004418e14;
    p.Re  = 6.378137e6;
    p.alt = altitude_km * 1e3;
    p.a   = p.Re + p.alt;
    p.n   = sqrt(p.mu / p.a^3);
    p.T   = 2*pi / p.n;
end
