function [f0, F_reg, g, xdot, u_effective] = Dynamics(x, params, theta, u, Delta)
% x=[phi;pitch;psi;p;q;r], theta=[Gam1;Gam2;Gam5;Gam6;Gam7;cd]

    x = x(:);
    if numel(x) ~= 6
        error('Dynamics:StateDimension', 'x must contain six states.');
    end

    if nargin < 3 || isempty(theta)
        theta = params.theta_true;
    end
    if nargin < 4 || isempty(u)
        u = zeros(4,1);
    end
    if nargin < 5 || isempty(Delta)
        Delta = zeros(4,4);
    end

    theta = theta(:);
    u = u(:);
    if numel(theta) ~= 6
        error('Dynamics:ParameterDimension', 'theta must be a 6-by-1 vector.');
    end
    if numel(u) ~= 4
        error('Dynamics:InputDimension', 'u must be a 4-by-1 vector.');
    end
    if ~isequal(size(Delta), [4,4])
        error('Dynamics:DeltaDimension', 'Delta must be a 4-by-4 matrix.');
    end
    if any(~isfinite([x; theta; u; Delta(:)]))
        error('Dynamics:NonFiniteInput', 'All model inputs must be finite.');
    end

    phi = x(1);
    pitch = x(2);
    p = x(4);
    q = x(5);
    r = x(6);
    omega = [p; q; r];

    J = euler_rate_matrix(phi, pitch);
    f0 = [J*omega; zeros(3,1)];
    F_omega = [ p*q, -q*r,  0,    0,             0,   -p;
                0,    0,    p*r, -(p^2 - r^2),   0,   -q;
               -q*r,  0,    0,    0,             p*q, -r];

    F_reg = [zeros(3,6); F_omega];
    g = [zeros(3,4); params.I \ params.B_alloc];
    u_effective = (eye(4) + Delta)*u;
    xdot = f0 + F_reg*theta + g*u_effective;
end


function J = euler_rate_matrix(phi, pitch)

    cphi = cos(phi);
    sphi = sin(phi);
    cpitch = cos(pitch);
    spitch = sin(pitch);
    if abs(cpitch) < 1e-6
        cpitch = sign(cpitch + eps)*1e-6;
    end

    tpitch = spitch/cpitch;
    J = [1, sphi*tpitch,  cphi*tpitch;
         0, cphi,        -sphi;
         0, sphi/cpitch,  cphi/cpitch];
end
