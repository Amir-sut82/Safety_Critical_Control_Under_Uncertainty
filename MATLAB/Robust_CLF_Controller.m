function [u,info] = Robust_CLF_Controller( ...
    x,theta_known,params,x_ref,delta_bound)
% The l1 epigraph conservatively bounds the spectral uncertainty term

    if nargin < 4 || isempty(x_ref)
        x_ref = params.x_eq;
    end
    if nargin < 5 || isempty(delta_bound)
        delta_bound = params.delta_max;
    end
    if delta_bound < 0
        error('Robust_CLF_Controller:InvalidDelta', ...
            'delta_bound must be nonnegative.');
    end

    theta_known = theta_known(:);
    if numel(theta_known) ~= params.n_theta
        error('Robust_CLF_Controller:ParameterDimension', ...
            'theta_known must contain %d parameters.',params.n_theta);
    end

    reg = Regressors(x,params,x_ref);
    a = reg.LgV;
    drift = reg.omega_V+reg.phi_V.'*theta_known;
    b_clf = -params.robust_gamma*reg.V-drift;
    beta = delta_bound*norm(a,2);

    n_u = size(a,2);
    if isfield(params,'robust_u_min')
        lower_bound = expand_bound(params.robust_u_min,n_u);
    else
        lower_bound = expand_bound(params.u_min,n_u);
    end
    if isfield(params,'robust_u_max')
        upper_bound = expand_bound(params.robust_u_max,n_u);
    else
        upper_bound = expand_bound(params.u_max,n_u);
    end
    if any(lower_bound > upper_bound)
        error('Robust_CLF_Controller:InvalidInputBounds', ...
            'Every robust lower input bound must not exceed its upper bound.');
    end
    z_upper = max(abs(lower_bound),abs(upper_bound));

    used_quadprog = false;
    exitflag = 0;
    u = [];
    z = [];
    solver_name = 'orthant dual-QP';

    allow_quadprog = isfield(params,'robust_use_quadprog') ...
        && params.robust_use_quadprog;
    if allow_quadprog && exist('quadprog','file') == 2
        try
            H = blkdiag(2*eye(n_u),2e-12*eye(n_u));
            A_qp = [a, beta*ones(1,n_u);
                    eye(n_u), -eye(n_u);
                   -eye(n_u), -eye(n_u)];
            b_qp = [b_clf; zeros(2*n_u,1)];
            lb = [lower_bound; zeros(n_u,1)];
            ub = [upper_bound; z_upper];
            options = optimoptions('quadprog','Display','off');
            [decision,~,qp_exitflag] = quadprog( ...
                H,zeros(2*n_u,1),A_qp,b_qp,[],[],lb,ub,[],options);
            if ~isempty(decision) && qp_exitflag > 0
                u = decision(1:n_u);
                z = decision(n_u+1:end);
                exitflag = qp_exitflag;
                used_quadprog = true;
                solver_name = 'quadprog dual-QP';
            end
        catch
            used_quadprog = false;
        end
    end

    tolerance = params.robust_qp_tolerance;
    if used_quadprog
        dual_left = a*u+beta*sum(abs(u));
        if dual_left > b_clf+tolerance
            used_quadprog = false;
            u = [];
            z = [];
        end
    end

    if ~used_quadprog
        [u,z,exitflag] = solve_l1_robust_qp( ...
            a,b_clf,beta,lower_bound,upper_bound,tolerance);
    end
    fallback_attempted = false;
    used_spectral_fallback = false;
    if exitflag <= 0
        fallback_attempted = true;
        [u_spectral,spectral_exitflag] = solve_spectral_projection( ...
            a,b_clf,beta,lower_bound,upper_bound,tolerance);
        u = u_spectral;
        z = abs(u);
        if spectral_exitflag > 0
            exitflag = spectral_exitflag;
            used_spectral_fallback = true;
            solver_name = 'exact spectral fallback';
        else
            exitflag = spectral_exitflag;
            solver_name = 'best-effort spectral control';
        end
    end

    dual_qp_residual = drift+a*u+beta*sum(abs(u)) ...
        + params.robust_gamma*reg.V;
    exact_robust_residual = drift+a*u ...
        + delta_bound*norm(a,2)*norm(u,2) ...
        + params.robust_gamma*reg.V;

    info.V = reg.V;
    info.reg = reg;
    info.drift = drift;
    info.LgV = a;
    info.delta_bound = delta_bound;
    info.beta = beta;
    info.z = z;
    info.dual_qp_residual = dual_qp_residual;
    info.exact_robust_residual = exact_robust_residual;
    info.conservative_qp_satisfied = dual_qp_residual <= tolerance;
    info.hard_robust_constraint_satisfied = ...
        exact_robust_residual <= tolerance;
    info.exitflag = exitflag;
    info.used_quadprog = used_quadprog;
    info.fallback_attempted = fallback_attempted;
    info.used_spectral_fallback = used_spectral_fallback;
    info.best_effort_control = exitflag <= 0;
    info.lower_bound = lower_bound;
    info.upper_bound = upper_bound;
    info.solver = solver_name;
end


function bound = expand_bound(value,n)
    if isscalar(value)
        bound = value*ones(n,1);
    else
        bound = value(:);
        if numel(bound) ~= n
            error('Robust_CLF_Controller:BoundDimension', ...
                'Each input bound must be scalar or %d-by-1.',n);
        end
    end
end


function [best_u,best_z,exitflag] = solve_l1_robust_qp( ...
    a,b,beta,lower_bound,upper_bound,tolerance)

    n = numel(lower_bound);
    best_cost = inf;
    best_u = min(max(zeros(n,1),lower_bound),upper_bound);
    best_z = abs(best_u);
    exitflag = -2;
    fully_unbounded = all(isinf(lower_bound) & lower_bound < 0) ...
        && all(isinf(upper_bound) & upper_bound > 0);

    for mask = 0:(2^n-1)
        signs = 2*bitget(mask,1:n)-1;
        signs = signs(:);
        orthant_lower = lower_bound;
        orthant_upper = upper_bound;

        positive = signs > 0;
        orthant_lower(positive) = max(orthant_lower(positive),0);
        orthant_upper(~positive) = min(orthant_upper(~positive),0);
        if any(orthant_lower > orthant_upper)
            continue;
        end

        effective_row = a+beta*signs.';
        if fully_unbounded
            if b >= -tolerance
                candidate = zeros(n,1);
                candidate_exitflag = 1;
            else
                direction = -effective_row.';
                direction(signs.*direction < 0) = 0;
                denominator = effective_row*direction;
                if denominator < -eps
                    candidate = (b/denominator)*direction;
                    candidate_exitflag = 1;
                else
                    candidate = zeros(n,1);
                    candidate_exitflag = -2;
                end
            end
        else
            [candidate,candidate_exitflag] = project_origin_halfspace_box( ...
                effective_row,b,orthant_lower,orthant_upper,tolerance);
        end
        if candidate_exitflag <= 0
            continue;
        end

        if a*candidate+beta*sum(abs(candidate)) > b+tolerance
            continue;
        end
        candidate_cost = candidate.'*candidate;
        if candidate_cost < best_cost
            best_cost = candidate_cost;
            best_u = candidate;
            best_z = abs(candidate);
            exitflag = 1;
        end
    end
end


function [u,exitflag] = project_origin_halfspace_box( ...
    a,b,lower_bound,upper_bound,tolerance)

    u0 = min(max(zeros(numel(lower_bound),1),lower_bound),upper_bound);
    if a*u0 <= b+tolerance
        u = u0;
        exitflag = 1;
        return;
    end
    u_min_linear = u0;
    positive_row = a > 0;
    negative_row = a < 0;
    u_min_linear(positive_row) = lower_bound(positive_row);
    u_min_linear(negative_row) = upper_bound(negative_row);
    if a*u_min_linear > b+tolerance
        u = u_min_linear;
        exitflag = -2;
        return;
    end

    multiplier_low = 0;
    multiplier_high = 1;
    u_multiplier = min(max(-multiplier_high*a.', ...
        lower_bound),upper_bound);
    while a*u_multiplier > b
        multiplier_high = 2*multiplier_high;
        u_multiplier = min(max(-multiplier_high*a.', ...
            lower_bound),upper_bound);
        if multiplier_high > 1e12
            u = u_min_linear;
            exitflag = -2;
            return;
        end
    end

    for k = 1:80
        multiplier_mid = 0.5*(multiplier_low+multiplier_high);
        u_multiplier = min(max(-multiplier_mid*a.', ...
            lower_bound),upper_bound);
        if a*u_multiplier > b
            multiplier_low = multiplier_mid;
        else
            multiplier_high = multiplier_mid;
        end
    end

    u = min(max(-multiplier_high*a.',lower_bound),upper_bound);
    exitflag = 1;
end


function [u,exitflag] = solve_spectral_projection( ...
    a,b,beta,lower_bound,upper_bound,tolerance)

    u0 = min(max(zeros(numel(lower_bound),1),lower_bound),upper_bound);
    if a*u0+beta*norm(u0,2) <= b+tolerance
        u = u0;
        exitflag = 1;
        return;
    end
    if norm(a,2) <= eps
        u = u0;
        exitflag = -2;
        return;
    end

    scale_low = 0;
    scale_high = 1;
    previous = nan(size(u0));
    while true
        candidate = min(max(-scale_high*a.',lower_bound),upper_bound);
        if a*candidate+beta*norm(candidate,2) <= b+tolerance
            break;
        end
        if all(isfinite(previous)) && norm(candidate-previous,inf) <= eps
            u = candidate;
            exitflag = -2;
            return;
        end
        previous = candidate;
        scale_high = 2*scale_high;
        if scale_high > 1e12
            u = candidate;
            exitflag = -2;
            return;
        end
    end

    for k = 1:80
        scale_mid = 0.5*(scale_low+scale_high);
        candidate = min(max(-scale_mid*a.',lower_bound),upper_bound);
        if a*candidate+beta*norm(candidate,2) > b
            scale_low = scale_mid;
        else
            scale_high = scale_mid;
        end
    end
    u = min(max(-scale_high*a.',lower_bound),upper_bound);
    exitflag = 1;
end
