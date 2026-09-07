function result = simulate_smid_case( ...
    use_online_bound,label,params,x0,x_ref,theta_known,Delta_true)
% Compare a fixed uncertainty bound with the online SMID bound

    dt = params.smid_dt;
    t = (0:dt:params.smid_t_final).';
    n_samples = numel(t);

    x = zeros(n_samples,6);
    u_clf = zeros(n_samples,4);
    u_nom = zeros(n_samples,4);
    u = zeros(n_samples,4);
    h = zeros(n_samples,1);
    psi1 = zeros(n_samples,1);
    delta_used = params.delta_max*ones(n_samples,1);
    delta_certified = params.delta_max*ones(n_samples,1);
    sigma_min = zeros(n_samples,1);
    lambda_min = zeros(n_samples,1);
    Delta_error = NaN(n_samples,1);
    actual_cbf_margin = zeros(n_samples,1);
    exact_set_cbf_margin = zeros(n_samples,1);
    conservative_cbf_margin = zeros(n_samples,1);
    safety_inflation = zeros(n_samples,1);
    clf_set_residual = zeros(n_samples,1);
    filter_intervention = zeros(n_samples,1);
    exitflag = ones(n_samples,1);
    best_effort = false(n_samples,1);
    full_rank = false(n_samples,1);
    smid_residual = zeros(n_samples,1);

    active = params.smid_active_channels(:).';
    Delta_active_true = Delta_true(active,active);
    true_delta_norm = norm(Delta_true,2);
    x(1,:) = x0(:).';
    smid_state = [];

    for k = 1:n_samples
        xk = x(k,:).';
        if use_online_bound && ~isempty(smid_state)
            controller_delta = smid_state.delta_bound;
        else
            controller_delta = params.delta_max;
        end

        [u_clf_k,clf_info] = Robust_CLF_Controller( ...
            xk,theta_known,params,x_ref,controller_delta);
        probe_k = smid_probe(t(k),params);
        u_nom_k = u_clf_k+probe_k;
        [uk,cbf_info] = Robust_CBF_SafetyFilter( ...
            xk,u_nom_k,theta_known,params,controller_delta);
        [~,~,~,xdot_k] = Dynamics( ...
            xk,params,theta_known,uk,Delta_true);

        [smid_state,smid_info] = SMID( ...
            smid_state,xk,uk,xdot_k,theta_known,params);

        reg = cbf_info.reg;
        barrier_drift = reg.Lf_psi1+reg.phi_h.'*theta_known;
        nominal_cbf_margin = barrier_drift+reg.Lg_psi1*uk ...
            + params.alpha2*reg.psi1;

        u_clf(k,:) = u_clf_k.';
        u_nom(k,:) = u_nom_k.';
        u(k,:) = uk.';
        h(k) = reg.h;
        psi1(k) = reg.psi1;
        delta_used(k) = controller_delta;
        delta_certified(k) = smid_info.delta_bound;
        sigma_min(k) = smid_info.sigma_min;
        lambda_min(k) = smid_info.lambda_min;
        Delta_error(k) = norm( ...
            smid_info.Delta_hat_active-Delta_active_true,2);
        actual_cbf_margin(k) = barrier_drift ...
            + reg.Lg_psi1*(eye(4)+Delta_true)*uk ...
            + params.alpha2*reg.psi1;
        exact_set_cbf_margin(k) = nominal_cbf_margin ...
            - controller_delta*norm(reg.Lg_psi1,2)*norm(uk,2);
        conservative_cbf_margin(k) = cbf_info.conservative_margin;
        safety_inflation(k) = controller_delta ...
            *norm(reg.Lg_psi1,2)*norm(uk,2);
        clf_set_residual(k) = clf_info.exact_robust_residual;
        filter_intervention(k) = norm(uk-u_nom_k,2);
        exitflag(k) = cbf_info.exitflag;
        best_effort(k) = cbf_info.best_effort_control;
        full_rank(k) = smid_info.full_rank;
        smid_residual(k) = smid_info.residual_max;

        if k < n_samples
            x(k+1,:) = rk4_uncertain_sample_hold( ...
                xk,uk,dt,params,theta_known,Delta_true).';
        end
    end

    result.label = label;
    result.use_online_bound = use_online_bound;
    result.t = t;
    result.x = x;
    result.u_clf = u_clf;
    result.u_nom = u_nom;
    result.u = u;
    result.h = h;
    result.psi1 = psi1;
    result.delta_used = delta_used;
    result.delta_certified = delta_certified;
    result.sigma_min = sigma_min;
    result.lambda_min = lambda_min;
    result.Delta_error = Delta_error;
    result.actual_cbf_margin = actual_cbf_margin;
    result.exact_set_cbf_margin = exact_set_cbf_margin;
    result.conservative_cbf_margin = conservative_cbf_margin;
    result.safety_inflation = safety_inflation;
    result.clf_set_residual = clf_set_residual;
    result.filter_intervention = filter_intervention;
    result.exitflag = exitflag;
    result.best_effort = best_effort;
    result.full_rank = full_rank;
    result.smid_residual = smid_residual;
    result.true_delta_norm = true_delta_norm;
    result.final_Delta_hat_active = smid_state.Delta_hat_active;
    result.final_delta_bound = delta_used(end);
    result.final_certified_bound = delta_certified(end);
    result.min_containment_margin = min(delta_used-true_delta_norm);
    result.min_h = min(h);
    result.min_psi1 = min(psi1);
    result.min_actual_cbf_margin = min(actual_cbf_margin);
    result.min_exact_set_cbf_margin = min(exact_set_cbf_margin);
    result.control_energy = trapz(t,sum(u.^2,2));
    result.intervention_energy = trapz(t,filter_intervention.^2);
    result.safety_margin_area = trapz(t,max(h,0));
    result.safety_inflation_area = trapz(t,safety_inflation);
    result.infeasible_count = sum(exitflag <= 0);
    result.best_effort_count = sum(best_effort);

    first_full_rank = find(full_rank,1,'first');
    if isempty(first_full_rank)
        result.full_rank_time = inf;
    else
        result.full_rank_time = t(first_full_rank);
    end
end


function probe = smid_probe(t,params)

    probe = zeros(4,1);
    if t > params.smid_probe_end
        return;
    end

    a = params.smid_probe_amplitude;
    envelope = sin(0.5*pi*min(t/0.25,1))^2;
    probe(2) = a*envelope*(sin(2.1*t)+0.30*sin(4.9*t+0.2));
    probe(3) = a*envelope*(sin(2.8*t+0.7)+0.25*cos(5.7*t));
    probe(4) = a*envelope*(cos(3.6*t-0.4)+0.20*sin(6.3*t));
end


function x_next = rk4_uncertain_sample_hold( ...
    x,u,dt,params,theta_known,Delta_true)

    k1 = plant_rhs(x,u,params,theta_known,Delta_true);
    k2 = plant_rhs(x+0.5*dt*k1,u,params,theta_known,Delta_true);
    k3 = plant_rhs(x+0.5*dt*k2,u,params,theta_known,Delta_true);
    k4 = plant_rhs(x+dt*k3,u,params,theta_known,Delta_true);
    x_next = x+(dt/6)*(k1+2*k2+2*k3+k4);
end


function dx = plant_rhs(x,u,params,theta_known,Delta_true)
    [~,~,~,dx] = Dynamics(x,params,theta_known,u,Delta_true);
end
