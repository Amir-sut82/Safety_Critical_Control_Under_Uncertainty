function [state,info] = SMID(state,x,u,xdot,theta_known,params)
% Identify the active 3x3 uncertainty block on channels 2:4

    x = x(:);
    u = u(:);
    xdot = xdot(:);
    theta_known = theta_known(:);
    active = params.smid_active_channels(:).';

    if numel(x) ~= 6 || numel(xdot) ~= 6
        error('SMID:StateDimension','x and xdot must be 6-by-1 vectors.');
    end
    if numel(u) ~= 4
        error('SMID:InputDimension','u must be a 4-by-1 vector.');
    end
    if numel(theta_known) ~= params.n_theta
        error('SMID:ParameterDimension', ...
            'theta_known must contain %d parameters.',params.n_theta);
    end
    if numel(active) ~= 3 || any(active < 1) || any(active > 4)
        error('SMID:ActiveChannels', ...
            'params.smid_active_channels must contain three input indices.');
    end

    if isempty(state)
        state.U = zeros(numel(active),0);
        state.Y = zeros(numel(active),0);
        state.Delta_hat_active = zeros(numel(active));
        state.delta_bound = params.delta_max;
        state.informative_count = 0;
    end

    [f0,F_reg,g] = Dynamics(x,params);
    G_active = g(4:6,active);
    if rcond(G_active) < 1e-12
        error('SMID:AllocationRank', ...
            'The active angular-acceleration allocation matrix is singular.');
    end

    u_active = u(active);
    nominal_omega_dot = f0(4:6)+F_reg(4:6,:)*theta_known ...
        + G_active*u_active;
    y = G_active\(xdot(4:6)-nominal_omega_dot);

    accepted = norm(u_active,2) >= params.smid_min_input_norm;
    if accepted
        state.U(:,end+1) = u_active;
        state.Y(:,end+1) = y;
        state.informative_count = state.informative_count+1;
    end

    sample_count = size(state.U,2);
    singular_values = svd(state.U,'econ');
    sigma_min = 0;
    lambda_min = 0;
    full_rank = false;
    radius = inf;
    residual_max = 0;

    if numel(singular_values) >= numel(active)
        sigma_min = singular_values(end);
        lambda_min = sigma_min^2;
        full_rank = sigma_min >= params.smid_sigma_threshold;
    end

    if full_rank
        Delta_hat = state.Y*pinv(state.U);
        residual = state.Y-Delta_hat*state.U;
        if ~isempty(residual)
            residual_max = max(sqrt(sum(residual.^2,1)));
        end
        radius = sqrt(sample_count)*params.smid_noise_bound/sigma_min;
        candidate_bound = norm(Delta_hat,2)+radius ...
            + params.smid_numerical_margin;
        candidate_bound = min(params.delta_max,candidate_bound);
        rate_limited_bound = state.delta_bound ...
            - params.smid_max_bound_decrease;
        state.delta_bound = min(state.delta_bound, ...
            max(candidate_bound,rate_limited_bound));
        state.Delta_hat_active = Delta_hat;
    end

    info.accepted = accepted;
    info.sample_count = sample_count;
    info.full_rank = full_rank;
    info.sigma_min = sigma_min;
    info.lambda_min = lambda_min;
    info.radius = radius;
    info.delta_bound = state.delta_bound;
    info.Delta_hat_active = state.Delta_hat_active;
    info.residual_max = residual_max;
    info.u_active = u_active;
    info.y = y;
end
