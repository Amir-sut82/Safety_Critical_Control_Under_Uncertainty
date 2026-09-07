function [u,info] = Robust_CBF_SafetyFilter( ...
    x,u_nom,theta_known,params,delta_bound)
% Robust relative-degree-two pitch safety filter

    if nargin < 5 || isempty(delta_bound)
        delta_bound = params.delta_max;
    end
    if delta_bound < 0
        error('Robust_CBF_SafetyFilter:InvalidDelta', ...
            'delta_bound must be nonnegative.');
    end

    u_nom = u_nom(:);
    theta_known = theta_known(:);
    if numel(u_nom) ~= 4
        error('Robust_CBF_SafetyFilter:InputDimension', ...
            'u_nom must be a 4-by-1 vector.');
    end
    if numel(theta_known) ~= params.n_theta
        error('Robust_CBF_SafetyFilter:ParameterDimension', ...
            'theta_known must contain %d parameters.',params.n_theta);
    end

    reg = Regressors(x,params,params.x_eq);
    a = reg.Lg_psi1;
    drift = reg.Lf_psi1+reg.phi_h.'*theta_known;
    beta = delta_bound*norm(a,2);
    b_cbf = drift+params.alpha2*reg.psi1;

    n_u = numel(u_nom);
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
        error('Robust_CBF_SafetyFilter:InvalidInputBounds', ...
            'Every robust lower input bound must not exceed its upper bound.');
    end
    z_upper = max(abs(lower_bound),abs(upper_bound));

    tolerance = params.robust_qp_tolerance;
    used_quadprog = false;
    exitflag = 0;
    u = [];
    z = [];
    solver_name = 'orthant robust-CBF QP';

    allow_quadprog = isfield(params,'robust_cbf_use_quadprog') ...
        && params.robust_cbf_use_quadprog;
    if allow_quadprog && exist('quadprog','file') == 2
        try
            H = blkdiag(2*eye(n_u),2e-12*eye(n_u));
            f_qp = [-2*u_nom; zeros(n_u,1)];
            A_qp = [-a, beta*ones(1,n_u);
                     eye(n_u), -eye(n_u);
                    -eye(n_u), -eye(n_u)];
            rhs_qp = [b_cbf; zeros(2*n_u,1)];
            lb = [lower_bound; zeros(n_u,1)];
            ub = [upper_bound; z_upper];
            options = optimoptions('quadprog','Display','off');
            [decision,~,qp_exitflag] = quadprog( ...
                H,f_qp,A_qp,rhs_qp,[],[],lb,ub,[],options);
            if ~isempty(decision) && qp_exitflag > 0
                u = decision(1:n_u);
                z = decision(n_u+1:end);
                exitflag = qp_exitflag;
                used_quadprog = true;
                solver_name = 'quadprog robust-CBF QP';
            end
        catch
            used_quadprog = false;
        end
    end

    if used_quadprog
        conservative_margin_test = b_cbf+a*u-beta*sum(abs(u));
        if conservative_margin_test < -tolerance
            used_quadprog = false;
            u = [];
            z = [];
        end
    end

    if ~used_quadprog
        [u,z,exitflag] = solve_l1_robust_cbf_qp( ...
            u_nom,a,b_cbf,beta,lower_bound,upper_bound,tolerance);
    end

    best_effort_control = false;
    if exitflag <= 0
        u = minimum_residual_control( ...
            u_nom,a,beta,lower_bound,upper_bound);
        z = abs(u);
        best_effort_control = true;
        solver_name = 'best-effort robust-CBF control';
    end

    nominal_margin = drift+a*u+params.alpha2*reg.psi1;
    exact_robust_margin = nominal_margin-beta*norm(u,2);
    conservative_margin = nominal_margin-beta*sum(abs(u));

    info.h = reg.h;
    info.hdot = reg.hdot;
    info.psi1 = reg.psi1;
    info.Lf_psi1 = reg.Lf_psi1;
    info.phi_h = reg.phi_h;
    info.Lg_psi1 = a;
    info.drift = drift;
    info.beta = beta;
    info.delta_bound = delta_bound;
    info.nominal_margin = nominal_margin;
    info.exact_robust_margin = exact_robust_margin;
    info.conservative_margin = conservative_margin;
    info.hard_robust_constraint_satisfied = ...
        exact_robust_margin >= -tolerance;
    info.conservative_qp_satisfied = conservative_margin >= -tolerance;
    info.exitflag = exitflag;
    info.used_quadprog = used_quadprog;
    info.best_effort_control = best_effort_control;
    info.lower_bound = lower_bound;
    info.upper_bound = upper_bound;
    info.z = z;
    info.reg = reg;
    info.solver = solver_name;
end


function bound = expand_bound(value,n)
    if isscalar(value)
        bound = value*ones(n,1);
    else
        bound = value(:);
        if numel(bound) ~= n
            error('Robust_CBF_SafetyFilter:BoundDimension', ...
                'Each input bound must be scalar or %d-by-1.',n);
        end
    end
end


function [best_u,best_z,exitflag] = solve_l1_robust_cbf_qp( ...
    u_nom,a,b,beta,lower_bound,upper_bound,tolerance)

    n = numel(u_nom);
    best_cost = inf;
    best_u = min(max(u_nom,lower_bound),upper_bound);
    best_z = abs(best_u);
    exitflag = -2;

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

        effective_row = -a+beta*signs.';
        [candidate,candidate_exitflag] = project_point_halfspace_box( ...
            u_nom,effective_row,b,orthant_lower,orthant_upper,tolerance);
        if candidate_exitflag <= 0
            continue;
        end
        if -a*candidate+beta*sum(abs(candidate)) > b+tolerance
            continue;
        end

        candidate_cost = sum((candidate-u_nom).^2);
        if candidate_cost < best_cost
            best_cost = candidate_cost;
            best_u = candidate;
            best_z = abs(candidate);
            exitflag = 1;
        end
    end
end


function [u,exitflag] = project_point_halfspace_box( ...
    v,a,b,lower_bound,upper_bound,tolerance)

    u0 = min(max(v,lower_bound),upper_bound);
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

    lambda_low = 0;
    lambda_high = 1;
    u_lambda = min(max(v-lambda_high*a.',lower_bound),upper_bound);
    while a*u_lambda > b
        lambda_high = 2*lambda_high;
        u_lambda = min(max(v-lambda_high*a.',lower_bound),upper_bound);
        if lambda_high > 1e12
            u = u_min_linear;
            exitflag = -2;
            return;
        end
    end

    for k = 1:70
        lambda_mid = 0.5*(lambda_low+lambda_high);
        u_lambda = min(max(v-lambda_mid*a.',lower_bound),upper_bound);
        if a*u_lambda > b
            lambda_low = lambda_mid;
        else
            lambda_high = lambda_mid;
        end
    end
    u = min(max(v-lambda_high*a.',lower_bound),upper_bound);
    exitflag = 1;
end


function u = minimum_residual_control( ...
    u_nom,a,beta,lower_bound,upper_bound)

    if any(~isfinite([lower_bound;upper_bound]))
        u = min(max(u_nom,lower_bound),upper_bound);
        return;
    end

    n = numel(u_nom);
    u = zeros(n,1);
    for index = 1:n
        zero_candidate = min(max(0,lower_bound(index)),upper_bound(index));
        candidates = [lower_bound(index);zero_candidate;upper_bound(index)];
        scores = a(index)*candidates-beta*abs(candidates);
        best_score = max(scores);
        best_indices = find(abs(scores-best_score) <= 1e-14);
        [~,nearest] = min(abs(candidates(best_indices)-u_nom(index)));
        u(index) = candidates(best_indices(nearest));
    end
end
