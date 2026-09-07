function reg = Regressors(x, params, x_ref)
% eISS-CLF and relative-degree-two pitch-HOCBF terms

    x = x(:);
    if numel(x) ~= 6
        error('Regressors:StateDimension','x must be a 6-by-1 vector.');
    end

    if nargin < 3 || isempty(x_ref)
        x_ref = params.x_eq;
    end
    x_ref = x_ref(:);
    if numel(x_ref) ~= 6
        error('Regressors:ReferenceDimension', ...
            'x_ref must be a 6-by-1 vector.');
    end

    phi = x(1);
    pitch = x(2);
    body_rates = x(4:6);

    [J, dJ_dphi, dJ_dpitch] = euler_kinematic_terms(phi,pitch);
    euler_rates = J*body_rates;
    Jdot = dJ_dphi*euler_rates(1) + dJ_dpitch*euler_rates(2);
    y = x(1:3)-x_ref(1:3);
    ydot = euler_rates;
    eta = [y; ydot];

    [~,F_reg,g] = Dynamics(x,params);
    F_omega = F_reg(4:6,:);
    G_omega = g(4:6,:);
    alpha0 = Jdot*body_rates;
    alphaF = J*F_omega;
    A_dec = J*G_omega;

    P = 0.5*(params.P+params.P.');
    grad_V_eta = 2*P*eta;

    reg.V = eta.'*P*eta;
    reg.omega_V = grad_V_eta(1:3).'*ydot + ...
        grad_V_eta(4:6).'*alpha0;
    reg.phi_V = alphaF.'*grad_V_eta(4:6);
    reg.LgV = grad_V_eta(4:6).'*A_dec;
    h = params.theta_max-pitch;
    pitch_rate = euler_rates(2);
    hdot = -pitch_rate;
    psi1 = hdot+params.alpha1*h;

    J2 = J(2,:);
    J2dot = Jdot(2,:);
    Lf_psi1 = -J2dot*body_rates-params.alpha1*pitch_rate;
    phi_h = (-J2*F_omega).';
    Lg_psi1 = -J2*G_omega;

    reg.h = h;
    reg.hdot = hdot;
    reg.psi1 = psi1;
    reg.Lf_psi1 = Lf_psi1;
    reg.phi_h = phi_h;
    reg.Lg_psi1 = Lg_psi1;

    reg.y = y;
    reg.ydot = ydot;
    reg.eta = eta;
    reg.J = J;
    reg.Jdot = Jdot;
    reg.alpha0 = alpha0;
    reg.alphaF = alphaF;
    reg.A_dec = A_dec;
    reg.F_reg = F_reg;
    reg.g = g;
end


function [J,dJ_dphi,dJ_dpitch] = euler_kinematic_terms(phi,pitch)

    cphi = cos(phi);
    sphi = sin(phi);
    cpitch = cos(pitch);
    spitch = sin(pitch);
    if abs(cpitch) < 1e-6
        cpitch = sign(cpitch+eps)*1e-6;
    end

    tpitch = spitch/cpitch;
    secpitch = 1/cpitch;

    J = [1, sphi*tpitch,       cphi*tpitch;
         0, cphi,             -sphi;
         0, sphi*secpitch,     cphi*secpitch];

    dJ_dphi = [0, cphi*tpitch,      -sphi*tpitch;
               0, -sphi,            -cphi;
               0, cphi*secpitch,    -sphi*secpitch];

    dJ_dpitch = [0, sphi*secpitch^2, ...
                    cphi*secpitch^2;
                  0, 0, 0;
                  0, sphi*secpitch*tpitch, ...
                    cphi*secpitch*tpitch];
end
