function animate_walking_ellipse3d(p, saveGif)
%ANIMATE_WALKING_ELLIPSE3D  3D animation of the CW relative orbit as xd0 is
%   swept through (and away from) the bounded-motion condition xd0 = 2 n z0.
%
%   animate_walking_ellipse3d()           uses a 400 km target orbit
%   animate_walking_ellipse3d(p)          uses the given cw_params struct
%   animate_walking_ellipse3d(p, true)    also writes walking_ellipse3d.gif
%
%   The chaser starts at fixed position (x0,y0,z0) = (0, 300, 300) m and
%   with a fixed cross-track velocity, so the out-of-plane (H-bar) motion
%   is a steady bounded oscillation throughout.  Only the along-track
%   velocity xd0 changes, as a fraction f of the bounded-motion value:
%
%       xd0 = f * (2 n z0)
%
%       f = 1.00 -> in-plane motion closes: a 3D bounded relative orbit
%       f > 1.00 -> in-plane ellipse walks; the 3D orbit becomes a
%                   drifting helix-like trajectory
%
%   With  delta = xd0 - 2 n z0  the in-plane centre drifts at
%   xd_walk = 6 n z0 - 3 xd0 = -3*delta, while the cross-track motion stays
%   bounded -- so a violation of the in-plane condition turns a closed 3D
%   relative orbit into a along-track-drifting one.
%
%   Axes: X = V-bar (along-track), Y = H-bar (cross-track),
%         Z = R-bar (radial); the R-bar axis is drawn pointing down so
%         "+z toward Earth" reads downward, matching the AE 598 note.

    if nargin < 1 || isempty(p),       p = cw_params(400);  end
    if nargin < 2 || isempty(saveGif), saveGif = false;     end

    % --- fixed start state (only xd0 will be swept) ---------------------
    x0 = 0;  y0 = 300;  z0 = 300;
    n  = p.n;  T = p.T;
    y0dot = n * y0;          % a clean bounded cross-track oscillation
    z0dot = 0;

    % --- sweep of xd0 as a fraction of the bounded-motion value ---------
    xd_bound = 2*n*z0;
    fracs = [linspace(0.75,1.00,26), ...     % climb to the bounded value
             linspace(1.00,1.30,55), ...     % violate it progressively
             linspace(1.30,1.00,40)];        % and sweep back

    t = linspace(0, 10*T, 1600);
    figure('Color','w','Position',[100 100 940 560]);

    for k = 1:numel(fracs)
        f     = fracs(k);
        xd0   = f * xd_bound;
        delta = xd0 - xd_bound;              % violation of the condition
        xw    = 6*n*z0 - 3*xd0;              % = -3*delta, the walk rate

        s0 = [x0, y0, z0, xd0, y0dot, z0dot];
        S  = cw_propagate3d(s0, t, n);
        x  = S(1,:);  y = S(2,:);  z = S(3,:);

        clf; hold on; grid on; box on;
        view(35,22);                         % 3D viewing angle
        set(gca,'ZDir','reverse');           % +z (R-bar, toward Earth) down
        xlim([-14000 4000]);
        ylim([-1200 1200]);
        zlim([-1400 1400]);

        % --- trajectory ------------------------------------------------
        plot3(x, y, z, 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);
        plot3(x(1), y(1), z(1), 'go', ...
              'MarkerFaceColor','g','MarkerSize',8);
        plot3(0,0,0,'ks','MarkerSize',12,'MarkerFaceColor','none', ...
              'LineWidth',1.6);
        text(0,0,0,'  target','FontSize',8);

        % --- projections onto the three coordinate planes -------------
        zfloor = 1400;  yfloor = 1200;
        plot3(x, y, zfloor*ones(size(z)), '-', ...     % onto R-bar floor
              'Color',[.7 .7 .85],'LineWidth',0.8);
        plot3(x, yfloor*ones(size(y)), z, '-', ...     % onto H-bar wall
              'Color',[.85 .8 .7],'LineWidth',0.8);

        % --- ellipse-centre markers, one per orbit --------------------
        for m = 0:9
            xc = xw * m * T;
            if xc > -14000 && xc < 4000
                plot3(xc, y0, z0, 'k+','MarkerSize',9,'LineWidth',1.2);
            end
        end

        xlabel('V-bar   x   [m]');
        ylabel('H-bar   y   [m]');
        zlabel('R-bar   z   [m]   (down)');

        % --- regime label ---------------------------------------------
        if abs(delta) < 1e-6
            regime = 'bounded-motion condition: closed 3D relative orbit';
        elseif delta < 0
            regime = 'xd_0 < 2 n z_0';
        else
            regime = sprintf(['in-plane ellipse walks %.0f m/orbit ' ...
                              '-> drifting 3D orbit'], xw*T);
        end
        title({sprintf(['xd_0 = %.2f \\times 2 n z_0     ' ...
                        '(\\delta = %.4f m/s,   xd_{walk} = -3\\delta ' ...
                        '= %.4f m/s)'], f, delta, xw), regime}, ...
              'FontWeight','bold','FontSize',10);

        drawnow;

        if saveGif
            frame = getframe(gcf);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);
            if k == 1
                imwrite(A,map,'walking_ellipse3d.gif','gif', ...
                        'LoopCount',Inf,'DelayTime',0.05);
            else
                imwrite(A,map,'walking_ellipse3d.gif','gif', ...
                        'WriteMode','append','DelayTime',0.05);
            end
        else
            pause(0.04);
        end
    end

    if saveGif
        fprintf('Saved walking_ellipse3d.gif\n');
    end
end