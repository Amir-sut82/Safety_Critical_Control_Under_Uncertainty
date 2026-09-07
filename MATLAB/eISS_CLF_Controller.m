function [u,info] = eISS_CLF_Controller(x,theta_hat,params,x_ref)
% Minimum-norm eISS-CLF QP with input bounds

    if nargin < 4 || isempty(x_ref)
        x_ref = params.x_eq;
    end
    theta_hat = theta_hat(:);
    if numel(theta_hat) ~= params.n_theta
        error('eISS_CLF_Controller:ParameterDimension', ...
            'theta_hat must contain %d parameters.',params.n_theta);
    end
    if params.kappa_iss <= 0
        error('eISS_CLF_Controller:InvalidKappa', ...
            'params.kappa_iss must be strictly positive.');
    end

    reg = Regressors(x,params,x_ref);
    iota = (reg.phi_V.'*reg.phi_V)/params.kappa_iss;

    A_clf = reg.LgV;
    b_clf = -params.c_eiss*reg.V + iota ...
        - reg.omega_V - reg.phi_V.'*theta_hat;

    n_u = size(A_clf,2);
    lower = expand_bound(params.u_min,n_u);
    upper = expand_bound(params.u_max,n_u);
    used_quadprog = false;
    exitflag = 0;
    u = [];
    allow_quadprog = true;
    if isfield(params,'use_quadprog')
        allow_quadprog = params.use_quadprog;
    end
    if allow_quadprog && exist('quadprog','file') == 2
        try
            options = optimoptions('quadprog','Display','off');
            [u,~,exitflag] = quadprog(2*eye(n_u),zeros(n_u,1), ...
                A_clf,b_clf,[],[],lower,upper,[],options);
            used_quadprog = ~isempty(u) && exitflag > 0;
        catch
            used_quadprog = false;
        end
    end
    if ~used_quadprog
        [u,exitflag] = project_origin_to_halfspace_box( ...
            A_clf,b_clf,lower,upper);
    end
    feasibility_tolerance = 1e-12;
    if A_clf*u-b_clf > feasibility_tolerance
        [u,projection_exitflag] = project_origin_to_halfspace_box( ...
            A_clf,b_clf,lower,upper);
        exitflag = projection_exitflag;
        used_quadprog = false;
    end

    estimated_Vdot = reg.omega_V + reg.phi_V.'*theta_hat + A_clf*u;
    constraint_residual = estimated_Vdot + params.c_eiss*reg.V - iota;

    info.V = reg.V;
    info.omega_V = reg.omega_V;
    info.phi_V = reg.phi_V;
    info.LgV = reg.LgV;
    info.iota = iota;
    info.A_clf = A_clf;
    info.b_clf = b_clf;
    info.estimated_Vdot = estimated_Vdot;
    info.constraint_residual = constraint_residual;
    info.exitflag = exitflag;
    info.used_quadprog = used_quadprog;
    info.reg = reg;
end


function bound = expand_bound(value,n)
    if isscalar(value)
        bound = value*ones(n,1);
    else
        bound = value(:);
        if numel(bound) ~= n
            error('eISS_CLF_Controller:BoundDimension', ...
                'Each input bound must be scalar or %d-by-1.',n);
        end
    end
end


function [u,exitflag] = project_origin_to_halfspace_box(a,b,lower,upper)

    tolerance = 1e-12;
    u0 = min(max(zeros(numel(lower),1),lower),upper);
    if a*u0 <= b+tolerance
        u = u0;
        exitflag = 1;
        return;
    end
    u_min_linear = upper;
    u_min_linear(a >= 0) = lower(a >= 0);
    if a*u_min_linear > b+tolerance
        u = u_min_linear;
        exitflag = -2;
        return;
    end
    lambda_low = 0;
    lambda_high = 1;
    u_lambda = min(max(-lambda_high*a.',lower),upper);
    while a*u_lambda > b
        lambda_high = 2*lambda_high;
        u_lambda = min(max(-lambda_high*a.',lower),upper);
        if lambda_high > 1e12
            u = u_min_linear;
            exitflag = -2;
            return;
        end
    end

    for k = 1:80
        lambda_mid = 0.5*(lambda_low+lambda_high);
        u_lambda = min(max(-lambda_mid*a.',lower),upper);
        if a*u_lambda > b
            lambda_low = lambda_mid;
        else
            lambda_high = lambda_mid;
        end
    end

    u = min(max(-lambda_high*a.',lower),upper);
    exitflag = 1;
end
