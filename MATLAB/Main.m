%% Mini-Project 3 - System model and open-loop validation
clear; clc; close all;

project_dir = fileparts(mfilename('fullpath'));
addpath(genpath(project_dir));

params = System_Parameters();

%% Actuation-uncertainty set and one admissible true matrix
delta_max = params.delta_max;
Delta_raw = [0, 0,     0,     0;
             0, -0.18, 0.03,  0;
             0, 0.02,  0.12, -0.02;
             0, 0,     0.02, -0.10];

if norm(Delta_raw,2) > delta_max
    Delta_true = Delta_raw*(delta_max/norm(Delta_raw,2));
else
    Delta_true = Delta_raw;
end

Delta_zero = zeros(4,4);
assert(norm(Delta_true,2) <= delta_max + 1e-12, ...
    'The selected Delta_true is outside the prescribed uncertainty set.');
assert(min(svd(eye(4) + Delta_true)) > 0, ...
    'I_4 + Delta_true must remain nonsingular for this validation.');

%% Validation 1: original closed-form model versus affine decomposition
t_equiv = (0:0.02:params.t_final).';
x0_equiv = [deg2rad(5); deg2rad(3); deg2rad(10); 0.5; -0.3; 0.2];
u_zero = @(t,x) zeros(4,1);

rhs_original = @(t,x) original_no_damping_rhs(x, u_zero(t,x), params);
rhs_affine = @(t,x) affine_rhs(x, u_zero(t,x), params, ...
    params.theta_nodamp, Delta_zero);

[t_original, x_original] = ode45(rhs_original, t_equiv, x0_equiv, ...
    params.ode_opts);
[t_affine, x_affine] = ode45(rhs_affine, t_equiv, x0_equiv, ...
    params.ode_opts);
x_affine_on_original_grid = interp1(t_affine, x_affine, t_original, 'pchip');
equivalence_error = x_affine_on_original_grid - x_original;
max_equivalence_error = max(abs(equivalence_error), [], 1);

%% Validation 2: isolate parametric and actuation uncertainty
t_uncertainty = (0:0.01:10).';
x0_uncertainty = x0_equiv;
u_open = @(t,x) 0.012*exp(-0.12*t)*[ ...
    0;
    sin(0.8*t);
    0.75*sin(1.1*t + 0.4);
    0.60*cos(0.7*t - 0.2)];

theta_nominal = params.theta_hat0;
theta_true = params.theta_true;

rhs_nominal = @(t,x) affine_rhs(x, u_open(t,x), params, ...
    theta_nominal, Delta_zero);
rhs_parametric = @(t,x) affine_rhs(x, u_open(t,x), params, ...
    theta_true, Delta_zero);
rhs_actuation = @(t,x) affine_rhs(x, u_open(t,x), params, ...
    theta_nominal, Delta_true);
rhs_combined = @(t,x) affine_rhs(x, u_open(t,x), params, ...
    theta_true, Delta_true);

[t_uncertainty, x_nominal] = ode45(rhs_nominal, t_uncertainty, ...
    x0_uncertainty, params.ode_opts);
[~, x_parametric] = ode45(rhs_parametric, t_uncertainty, ...
    x0_uncertainty, params.ode_opts);
[~, x_actuation] = ode45(rhs_actuation, t_uncertainty, ...
    x0_uncertainty, params.ode_opts);
[~, x_combined] = ode45(rhs_combined, t_uncertainty, ...
    x0_uncertainty, params.ode_opts);

U_commanded = zeros(numel(t_uncertainty),4);
U_effective = zeros(numel(t_uncertainty),4);
for k = 1:numel(t_uncertainty)
    uk = u_open(t_uncertainty(k), x_combined(k,:).');
    U_commanded(k,:) = uk.';
    U_effective(k,:) = ((eye(4) + Delta_true)*uk).';
end

%% Task 1.1: eISS-CLF minimum-norm QP (Delta = 0)
t_eiss = (0:0.01:15).';
x0_eiss = [deg2rad(10); 0; deg2rad(30); 0.5; 0.4; 0];
x_ref_eiss = params.x_eq;
u_reg_test = [0; 0.01; -0.008; 0.006];
reg_test = Regressors(x0_eiss,params,x_ref_eiss);
[~,~,~,dx_reg_test] = Dynamics(x0_eiss,params,params.theta_true, ...
    u_reg_test,Delta_zero);
epsilon_reg = 1e-6;
V_plus = Regressors(x0_eiss+epsilon_reg*dx_reg_test, ...
    params,x_ref_eiss).V;
V_minus = Regressors(x0_eiss-epsilon_reg*dx_reg_test, ...
    params,x_ref_eiss).V;
Vdot_finite_difference = (V_plus-V_minus)/(2*epsilon_reg);
Vdot_regressor = reg_test.omega_V ...
    + reg_test.phi_V.'*params.theta_true + reg_test.LgV*u_reg_test;
regressor_validation_error = abs(Vdot_finite_difference-Vdot_regressor);
assert(regressor_validation_error < 1e-7, ...
    'The CLF regressor decomposition failed its finite-difference check.');

theta_hat_cases = {params.theta_true, params.theta_hat0};
case_labels = {'Perfect estimate: $\hat{\theta}=\theta^*$', ...
    'Fixed poor estimate: $\hat{\theta}=\hat{\theta}(0)$'};
eiss_clf = repmat(struct(),numel(theta_hat_cases),1);

for ic = 1:numel(theta_hat_cases)
    theta_hat_fixed = theta_hat_cases{ic};
    rhs_eiss = @(t,x) eiss_closed_loop_rhs(x,theta_hat_fixed, ...
        params,x_ref_eiss,Delta_zero);

    [t_case,x_case] = ode45(rhs_eiss,t_eiss,x0_eiss,params.ode_opts);
    n_case = numel(t_case);
    u_case = zeros(n_case,4);
    V_case = zeros(n_case,1);
    iota_case = zeros(n_case,1);
    residual_case = zeros(n_case,1);
    true_Vdot_case = zeros(n_case,1);
    exitflag_case = zeros(n_case,1);

    for k = 1:n_case
        [uk,controller_info] = eISS_CLF_Controller( ...
            x_case(k,:).',theta_hat_fixed,params,x_ref_eiss);
        reg_k = controller_info.reg;

        u_case(k,:) = uk.';
        V_case(k) = controller_info.V;
        iota_case(k) = controller_info.iota;
        residual_case(k) = controller_info.constraint_residual;
        true_Vdot_case(k) = reg_k.omega_V ...
            + reg_k.phi_V.'*params.theta_true + reg_k.LgV*uk;
        exitflag_case(k) = controller_info.exitflag;
    end

    eiss_clf(ic).label = case_labels{ic};
    eiss_clf(ic).theta_hat = theta_hat_fixed;
    eiss_clf(ic).t = t_case;
    eiss_clf(ic).x = x_case;
    eiss_clf(ic).u = u_case;
    eiss_clf(ic).V = V_case;
    eiss_clf(ic).iota = iota_case;
    eiss_clf(ic).constraint_residual = residual_case;
    eiss_clf(ic).true_Vdot = true_Vdot_case;
    eiss_clf(ic).exitflag = exitflag_case;

    if any(exitflag_case <= 0)
        warning('Main:eISSQPInfeasible', ...
            'The eISS-CLF QP was infeasible in case %d.',ic);
    end
end

%% Task 1.2: compare four history-stack parameter estimators
estimator_functions = {@Gradient_Estimator,@RLS_Estimator, ...
    @RLS_ConstForget,@RLS_VariableForget};
estimator_labels = {'Gradient + history stack','RLS + history stack', ...
    'RLS + constant forgetting','RLS + variable forgetting'};

estimator_results = struct([]);
for method_index = 1:numel(estimator_functions)
    method_result = simulate_estimator_case( ...
        estimator_functions{method_index},estimator_labels{method_index}, ...
        params,x0_eiss,x_ref_eiss);
    if method_index == 1
        estimator_results = repmat(method_result,numel(estimator_functions),1);
    end
    estimator_results(method_index) = method_result;
end

%% Task 1.3: ISSf-CBF safety filter (Delta = 0)
t_cbf = (0:0.01:15).';
x0_cbf = [deg2rad(-5); deg2rad(15); deg2rad(-20); 0.2; -0.1; 0.3];
x_ref_cbf = params.x_ref;
theta_hat_cbf = params.theta_hat0;
cbf_sim_params = params;
cbf_sim_params.use_quadprog = false;
u_cbf_test = [0; -0.006; 0.009; -0.004];
reg_cbf_test = Regressors(x0_cbf,cbf_sim_params,x_ref_cbf);
[~,~,~,dx_cbf_test] = Dynamics(x0_cbf,cbf_sim_params, ...
    params.theta_true,u_cbf_test,Delta_zero);
epsilon_cbf = 1e-6;
psi1_plus = Regressors(x0_cbf+epsilon_cbf*dx_cbf_test, ...
    cbf_sim_params,x_ref_cbf).psi1;
psi1_minus = Regressors(x0_cbf-epsilon_cbf*dx_cbf_test, ...
    cbf_sim_params,x_ref_cbf).psi1;
psi1dot_finite_difference = (psi1_plus-psi1_minus)/(2*epsilon_cbf);
psi1dot_regressor = reg_cbf_test.Lf_psi1 ...
    + reg_cbf_test.phi_h.'*params.theta_true ...
    + reg_cbf_test.Lg_psi1*u_cbf_test;
cbf_regressor_validation_error = abs( ...
    psi1dot_finite_difference-psi1dot_regressor);
assert(cbf_regressor_validation_error < 1e-7, ...
    'The HOCBF regressor decomposition failed its finite-difference check.');

cbf_modes = {'none','naive','issf'};
cbf_labels = {'Unsafe nominal eISS-CLF', ...
    'Naive certainty-equivalence CBF', ...
    'ISSf-CBF safety filter'};
task1_3 = repmat(struct(),numel(cbf_modes),1);

for ic = 1:numel(cbf_modes)
    mode = cbf_modes{ic};
    rhs_cbf = @(t,x) issf_closed_loop_rhs(x,mode,theta_hat_cbf, ...
        cbf_sim_params,x_ref_cbf,Delta_zero);
    [t_case,x_case] = ode45(rhs_cbf,t_cbf,x0_cbf, ...
        cbf_sim_params.ode_opts);

    n_case = numel(t_case);
    u_nom_case = zeros(n_case,4);
    u_case = zeros(n_case,4);
    h_case = zeros(n_case,1);
    psi1_case = zeros(n_case,1);
    iota_h_case = zeros(n_case,1);
    estimated_margin_case = nan(n_case,1);
    true_margin_case = zeros(n_case,1);
    exitflag_case = ones(n_case,1);

    for k = 1:n_case
        xk = x_case(k,:).';
        [u_nom_k,~] = eISS_CLF_Controller( ...
            xk,theta_hat_cbf,cbf_sim_params,x_ref_cbf);
        reg_k = Regressors(xk,cbf_sim_params,x_ref_cbf);

        if strcmp(mode,'none')
            uk = u_nom_k;
            iota_h_k = 0;
            estimated_margin_k = reg_k.Lf_psi1 ...
                + reg_k.phi_h.'*theta_hat_cbf + reg_k.Lg_psi1*uk ...
                + cbf_sim_params.alpha2*reg_k.psi1;
            exitflag_k = 1;
        else
            [uk,filter_info] = ISSf_CBF_SafetyFilter( ...
                xk,u_nom_k,theta_hat_cbf,cbf_sim_params,mode);
            iota_h_k = filter_info.iota_h;
            estimated_margin_k = filter_info.estimated_margin;
            exitflag_k = filter_info.exitflag;
        end

        true_margin_k = reg_k.Lf_psi1 ...
            + reg_k.phi_h.'*params.theta_true + reg_k.Lg_psi1*uk ...
            + cbf_sim_params.alpha2*reg_k.psi1;

        u_nom_case(k,:) = u_nom_k.';
        u_case(k,:) = uk.';
        h_case(k) = reg_k.h;
        psi1_case(k) = reg_k.psi1;
        iota_h_case(k) = iota_h_k;
        estimated_margin_case(k) = estimated_margin_k;
        true_margin_case(k) = true_margin_k;
        exitflag_case(k) = exitflag_k;
    end

    task1_3(ic).mode = mode;
    task1_3(ic).label = cbf_labels{ic};
    task1_3(ic).t = t_case;
    task1_3(ic).x = x_case;
    task1_3(ic).u_nom = u_nom_case;
    task1_3(ic).u = u_case;
    task1_3(ic).h = h_case;
    task1_3(ic).psi1 = psi1_case;
    task1_3(ic).iota_h = iota_h_case;
    task1_3(ic).estimated_margin = estimated_margin_case;
    task1_3(ic).true_margin = true_margin_case;
    task1_3(ic).exitflag = exitflag_case;
    task1_3(ic).intervention_norm = sqrt(sum((u_case-u_nom_case).^2,2));
    task1_3(ic).min_h = min(h_case);
    task1_3(ic).control_energy = trapz(t_case,sum(u_case.^2,2));

    if ~strcmp(mode,'none') && any(exitflag_case <= 0)
        warning('Main:CBFQPInfeasible', ...
            'The %s safety-filter QP was infeasible at %d samples.', ...
            mode,sum(exitflag_case <= 0));
    end
end

if task1_3(3).min_h < -1e-6
    warning('Main:ISSfSafetyViolation', ...
        'The ISSf-CBF trajectory violated h>=0 by %.3e rad.', ...
        -task1_3(3).min_h);
end

%% Task 1.4: estimator comparison with eISS-CLF + ISSf-CBF
task1_4 = struct([]);
for method_index = 1:numel(estimator_functions)
    method_result = simulate_estimator_issf_case( ...
        estimator_functions{method_index},estimator_labels{method_index}, ...
        params,x0_cbf,x_ref_cbf);
    if method_index == 1
        task1_4 = repmat(method_result,numel(estimator_functions),1);
    end
    task1_4(method_index) = method_result;
end

%% Task 2.1: robust CLF under multiplicative actuation uncertainty
t_robust = (0:0.01:params.robust_t_final).';
x0_robust = x0_eiss;
x_ref_robust = params.x_eq;
theta_robust = params.theta_true;

robust_case_labels = { ...
    'Nominal CLF, $\Delta=0$', ...
    'Nominal CLF, adversarial $\Delta_{\mathrm{wc}}$', ...
    'Robust CLF, $\Delta=0$', ...
    'Robust CLF, fixed $\Delta^*$', ...
    'Robust CLF, adversarial $\Delta_{\mathrm{wc}}$'};
controller_delta_bounds = [0,0,delta_max,delta_max,delta_max];
uncertainty_modes = {'zero','worst','zero','fixed','worst'};
robust_clf = repmat(struct(),numel(robust_case_labels),1);

for ic = 1:numel(robust_case_labels)
    controller_delta = controller_delta_bounds(ic);
    uncertainty_mode = uncertainty_modes{ic};
    rhs_robust = @(t,x) robust_clf_closed_loop_rhs( ...
        x,controller_delta,uncertainty_mode,theta_robust,params, ...
        x_ref_robust,Delta_true);
    [t_case,x_case] = ode45(rhs_robust,t_robust,x0_robust,params.ode_opts);

    n_case = numel(t_case);
    u_case = zeros(n_case,4);
    V_case = zeros(n_case,1);
    actual_decay_residual = zeros(n_case,1);
    exact_set_residual = zeros(n_case,1);
    dual_set_residual = zeros(n_case,1);
    Delta_norm_case = zeros(n_case,1);
    exitflag_case = zeros(n_case,1);
    conservative_qp_case = false(n_case,1);
    fallback_attempted_case = false(n_case,1);
    spectral_fallback_case = false(n_case,1);
    best_effort_case = false(n_case,1);

    for k = 1:n_case
        xk = x_case(k,:).';
        [uk,controller_info] = Robust_CLF_Controller( ...
            xk,theta_robust,params,x_ref_robust,controller_delta);
        Delta_k = select_actuation_uncertainty( ...
            uncertainty_mode,controller_info.LgV,uk,delta_max,Delta_true);
        reg_k = controller_info.reg;

        actual_Vdot = reg_k.omega_V+reg_k.phi_V.'*theta_robust ...
            + reg_k.LgV*(eye(4)+Delta_k)*uk;
        exact_worst_Vdot = reg_k.omega_V+reg_k.phi_V.'*theta_robust ...
            + reg_k.LgV*uk ...
            + delta_max*norm(reg_k.LgV,2)*norm(uk,2);
        dual_worst_Vdot = reg_k.omega_V+reg_k.phi_V.'*theta_robust ...
            + reg_k.LgV*uk ...
            + delta_max*norm(reg_k.LgV,2)*sum(abs(uk));

        u_case(k,:) = uk.';
        V_case(k) = reg_k.V;
        actual_decay_residual(k) = actual_Vdot ...
            + params.robust_gamma*reg_k.V;
        exact_set_residual(k) = exact_worst_Vdot ...
            + params.robust_gamma*reg_k.V;
        dual_set_residual(k) = dual_worst_Vdot ...
            + params.robust_gamma*reg_k.V;
        Delta_norm_case(k) = norm(Delta_k,2);
        exitflag_case(k) = controller_info.exitflag;
        conservative_qp_case(k) = controller_info.conservative_qp_satisfied;
        fallback_attempted_case(k) = controller_info.fallback_attempted;
        spectral_fallback_case(k) = controller_info.used_spectral_fallback;
        best_effort_case(k) = controller_info.best_effort_control;
    end

    robust_clf(ic).label = robust_case_labels{ic};
    robust_clf(ic).controller_delta = controller_delta;
    robust_clf(ic).uncertainty_mode = uncertainty_mode;
    robust_clf(ic).t = t_case;
    robust_clf(ic).x = x_case;
    robust_clf(ic).state_norm = sqrt(sum(x_case.^2,2));
    robust_clf(ic).u = u_case;
    robust_clf(ic).V = V_case;
    robust_clf(ic).actual_decay_residual = actual_decay_residual;
    robust_clf(ic).exact_set_residual = exact_set_residual;
    robust_clf(ic).dual_set_residual = dual_set_residual;
    robust_clf(ic).Delta_norm = Delta_norm_case;
    robust_clf(ic).exitflag = exitflag_case;
    robust_clf(ic).conservative_qp_satisfied = conservative_qp_case;
    robust_clf(ic).fallback_attempted = fallback_attempted_case;
    robust_clf(ic).used_spectral_fallback = spectral_fallback_case;
    robust_clf(ic).best_effort_control = best_effort_case;
    robust_clf(ic).control_energy = trapz(t_case,sum(u_case.^2,2));
    robust_clf(ic).final_V_ratio = V_case(end)/max(V_case(1),eps);
    robust_clf(ic).max_actual_decay_residual = max(actual_decay_residual);
    robust_clf(ic).max_exact_set_residual = max(exact_set_residual);
    robust_clf(ic).fallback_attempt_count = sum(fallback_attempted_case);
    robust_clf(ic).spectral_fallback_count = sum(spectral_fallback_case);
    robust_clf(ic).best_effort_count = sum(best_effort_case);

    if controller_delta > 0
        if any(exitflag_case <= 0)
            warning('Main:RobustCLFInfeasible', ...
                'The robust CLF QP was infeasible in case %d.',ic);
        end
        if max(exact_set_residual) > 1e-8
            warning('Main:RobustCLFCertificateViolation', ...
                'The exact spectral robust residual reached %.3e in case %d.', ...
                max(exact_set_residual),ic);
        end
        if max(actual_decay_residual) > 1e-8
            warning('Main:RobustCLFDecayViolation', ...
                'The actual CLF decay residual reached %.3e in case %d.', ...
                max(actual_decay_residual),ic);
        end
        if any(spectral_fallback_case)
            warning('Main:RobustCLFSpectralFallback', ...
                ['The conservative l1 QP was infeasible and the exact ', ...
                 'spectral fallback was used successfully at %d samples ', ...
                 'in case %d.'],sum(spectral_fallback_case),ic);
        end
        if any(best_effort_case)
            warning('Main:RobustCLFBestEffort', ...
                ['No hard robust solution existed at %d samples in case ', ...
                 '%d; minimum-residual control was applied.'], ...
                sum(best_effort_case),ic);
        end
    end
end

%% Task 2.2: robust HOCBF safety filter under actuation uncertainty
t_robust_cbf = (0:0.01:params.robust_cbf_t_final).';
x0_robust_cbf = x0_cbf;
x_ref_robust_cbf = params.x_ref;
theta_robust_cbf = params.theta_true;

robust_cbf_labels = { ...
    'Unsafe nominal CLF, adversarial $\Delta_{h,\mathrm{wc}}$', ...
    'Nominal CBF, adversarial $\Delta_{h,\mathrm{wc}}$', ...
    'Robust CBF, $\Delta=0$', ...
    'Robust CBF, fixed $\Delta^*$', ...
    'Robust CBF, adversarial $\Delta_{h,\mathrm{wc}}$'};
robust_cbf_filter_modes = {'none','filter','filter','filter','filter'};
robust_cbf_delta_bounds = [0,0,delta_max,delta_max,delta_max];
robust_cbf_uncertainty_modes = {'worst','worst','zero','fixed','worst'};
robust_cbf = repmat(struct(),numel(robust_cbf_labels),1);

for ic = 1:numel(robust_cbf_labels)
    filter_mode = robust_cbf_filter_modes{ic};
    filter_delta = robust_cbf_delta_bounds(ic);
    uncertainty_mode = robust_cbf_uncertainty_modes{ic};
    rhs_robust_cbf = @(t,x) robust_cbf_closed_loop_rhs( ...
        x,filter_mode,filter_delta,uncertainty_mode,theta_robust_cbf, ...
        params,x_ref_robust_cbf,Delta_true);
    [t_case,x_case] = ode45(rhs_robust_cbf,t_robust_cbf, ...
        x0_robust_cbf,params.ode_opts);

    n_case = numel(t_case);
    u_nom_case = zeros(n_case,4);
    u_case = zeros(n_case,4);
    h_case = zeros(n_case,1);
    psi1_case = zeros(n_case,1);
    actual_margin_case = zeros(n_case,1);
    exact_set_margin_case = zeros(n_case,1);
    dual_set_margin_case = zeros(n_case,1);
    Delta_norm_case = zeros(n_case,1);
    exitflag_case = ones(n_case,1);
    best_effort_case = false(n_case,1);

    for k = 1:n_case
        xk = x_case(k,:).';
        [u_nom_k,~] = Robust_CLF_Controller( ...
            xk,theta_robust_cbf,params,x_ref_robust_cbf,0);
        reg_k = Regressors(xk,params,x_ref_robust_cbf);
        if strcmp(filter_mode,'none')
            uk = u_nom_k;
            exitflag_k = 1;
            best_effort_k = false;
        else
            [uk,filter_info] = Robust_CBF_SafetyFilter( ...
                xk,u_nom_k,theta_robust_cbf,params,filter_delta);
            exitflag_k = filter_info.exitflag;
            best_effort_k = filter_info.best_effort_control;
        end

        Delta_k = select_cbf_actuation_uncertainty( ...
            uncertainty_mode,reg_k.Lg_psi1,uk,delta_max,Delta_true);
        barrier_drift = reg_k.Lf_psi1 ...
            + reg_k.phi_h.'*theta_robust_cbf;
        nominal_margin = barrier_drift+reg_k.Lg_psi1*uk ...
            + params.alpha2*reg_k.psi1;
        actual_margin = barrier_drift ...
            + reg_k.Lg_psi1*(eye(4)+Delta_k)*uk ...
            + params.alpha2*reg_k.psi1;
        exact_set_margin = nominal_margin ...
            - delta_max*norm(reg_k.Lg_psi1,2)*norm(uk,2);
        dual_set_margin = nominal_margin ...
            - delta_max*norm(reg_k.Lg_psi1,2)*sum(abs(uk));

        u_nom_case(k,:) = u_nom_k.';
        u_case(k,:) = uk.';
        h_case(k) = reg_k.h;
        psi1_case(k) = reg_k.psi1;
        actual_margin_case(k) = actual_margin;
        exact_set_margin_case(k) = exact_set_margin;
        dual_set_margin_case(k) = dual_set_margin;
        Delta_norm_case(k) = norm(Delta_k,2);
        exitflag_case(k) = exitflag_k;
        best_effort_case(k) = best_effort_k;
    end

    robust_cbf(ic).label = robust_cbf_labels{ic};
    robust_cbf(ic).filter_mode = filter_mode;
    robust_cbf(ic).filter_delta = filter_delta;
    robust_cbf(ic).uncertainty_mode = uncertainty_mode;
    robust_cbf(ic).t = t_case;
    robust_cbf(ic).x = x_case;
    robust_cbf(ic).u_nom = u_nom_case;
    robust_cbf(ic).u = u_case;
    robust_cbf(ic).h = h_case;
    robust_cbf(ic).psi1 = psi1_case;
    robust_cbf(ic).actual_margin = actual_margin_case;
    robust_cbf(ic).exact_set_margin = exact_set_margin_case;
    robust_cbf(ic).dual_set_margin = dual_set_margin_case;
    robust_cbf(ic).Delta_norm = Delta_norm_case;
    robust_cbf(ic).exitflag = exitflag_case;
    robust_cbf(ic).best_effort_control = best_effort_case;
    robust_cbf(ic).intervention_norm = ...
        sqrt(sum((u_case-u_nom_case).^2,2));
    robust_cbf(ic).min_h = min(h_case);
    robust_cbf(ic).min_psi1 = min(psi1_case);
    robust_cbf(ic).min_actual_margin = min(actual_margin_case);
    robust_cbf(ic).min_exact_set_margin = min(exact_set_margin_case);
    robust_cbf(ic).control_energy = trapz(t_case,sum(u_case.^2,2));
    robust_cbf(ic).best_effort_count = sum(best_effort_case);

    if filter_delta > 0
        if any(exitflag_case <= 0)
            warning('Main:RobustCBFInfeasible', ...
                'The robust CBF QP was infeasible in case %d.',ic);
        end
        if min(exact_set_margin_case) < -1e-8
            warning('Main:RobustCBFCertificateViolation', ...
                ['The exact spectral robust-CBF margin reached %.3e ', ...
                 'in case %d.'],min(exact_set_margin_case),ic);
        end
        if min(actual_margin_case) < -1e-8
            warning('Main:RobustCBFMarginViolation', ...
                'The actual HOCBF margin reached %.3e in case %d.', ...
                min(actual_margin_case),ic);
        end
        if min(h_case) < -1e-6
            warning('Main:RobustCBFSafetyViolation', ...
                'The robust-CBF trajectory violated h>=0 in case %d.',ic);
        end
    end
end

if robust_cbf(1).min_h >= 0
    warning('Main:RobustCBFScenarioNotActive', ...
        ['The selected nominal-CLF/worst-case scenario did not violate ', ...
         'safety; adjust the unsafe reference or initial state.']);
end

%% Task 2.3: set-membership identification for uncertainty reduction
active_smid = params.smid_active_channels(:).';
inactive_smid = setdiff(1:4,active_smid);
if norm(Delta_true(inactive_smid,:),'fro') > 1e-12 ...
        || norm(Delta_true(:,inactive_smid),'fro') > 1e-12
    error('Main:SMIDStructureMismatch', ...
        ['Task 2.3 assumes the unallocated channel has zero uncertainty ', ...
         'row and column.']);
end

smid_labels = { ...
    'Fixed initial set, $\delta=\delta_{\max}$', ...
    'Online SMID-reduced set, $\delta=\hat\delta(t)$'};
task2_3 = struct([]);
x_ref_smid = x_ref_robust_cbf;
x_ref_smid(2) = params.smid_theta_ref;
for ic = 1:numel(smid_labels)
    use_online_bound = ic == 2;
    smid_result = simulate_smid_case( ...
        use_online_bound,smid_labels{ic},params,x0_robust_cbf, ...
        x_ref_smid,theta_robust_cbf,Delta_true);
    if ic == 1
        task2_3 = repmat(smid_result,numel(smid_labels),1);
    end
    task2_3(ic) = smid_result;
end

online_smid = task2_3(2);
if ~isfinite(online_smid.full_rank_time)
    warning('Main:SMIDInsufficientExcitation', ...
        'The Task 2.3 data matrix never became full rank.');
end
if online_smid.final_delta_bound >= delta_max-1e-6
    warning('Main:SMIDNoContraction', ...
        'The online set-membership bound did not contract.');
end
if online_smid.min_containment_margin < -1e-8
    warning('Main:SMIDContainmentViolation', ...
        'The SMID bound excluded the true uncertainty by %.3e.', ...
        -online_smid.min_containment_margin);
end
if online_smid.min_h < -1e-6
    warning('Main:SMIDSafetyViolation', ...
        'The SMID robust loop violated h>=0 by %.3e rad.', ...
        -online_smid.min_h);
end
if online_smid.min_exact_set_cbf_margin < -1e-8
    warning('Main:SMIDCertificateViolation', ...
        'The SMID robust-CBF margin reached %.3e.', ...
        online_smid.min_exact_set_cbf_margin);
end
if online_smid.control_energy >= task2_3(1).control_energy
    warning('Main:SMIDNoControlEffortReduction', ...
        'The selected SMID scenario did not reduce integrated control effort.');
end
if online_smid.safety_inflation_area >= task2_3(1).safety_inflation_area
    warning('Main:SMIDNoInflationReduction', ...
        'The selected SMID scenario did not reduce robust-CBF inflation.');
end

all_states = [x_original; x_affine; x_nominal; x_parametric; ...
              x_actuation; x_combined];
for ic = 1:numel(task1_3)
    all_states = [all_states; task1_3(ic).x];
end
for method_index = 1:numel(task1_4)
    all_states = [all_states; task1_4(method_index).x];
end
for ic = 1:numel(robust_clf)
    all_states = [all_states; robust_clf(ic).x];
end
for ic = 1:numel(robust_cbf)
    all_states = [all_states; robust_cbf(ic).x];
end
for ic = 1:numel(task2_3)
    all_states = [all_states; task2_3(ic).x];
end
assert(all(isfinite(all_states(:))), ...
    'A validation trajectory contains NaN or Inf.');

%% Collect, print, save, and plot results
results.t_equiv = t_original;
results.x_original = x_original;
results.x_affine = x_affine_on_original_grid;
results.equivalence_error = equivalence_error;
results.max_equivalence_error = max_equivalence_error;

results.t_uncertainty = t_uncertainty;
results.x_nominal = x_nominal;
results.x_parametric = x_parametric;
results.x_actuation = x_actuation;
results.x_combined = x_combined;
results.U_commanded = U_commanded;
results.U_effective = U_effective;
results.Delta_true = Delta_true;
results.delta_max = delta_max;
results.theta_nominal = theta_nominal;
results.theta_true = theta_true;
results.output_dir = project_dir;
results.eiss_clf = eiss_clf;
results.estimator_comparison = estimator_results;
results.task1_3 = task1_3;
results.task1_4 = task1_4;
results.robust_clf = robust_clf;
results.robust_cbf = robust_cbf;
results.task2_3 = task2_3;
results.regressor_validation.Vdot_finite_difference = ...
    Vdot_finite_difference;
results.regressor_validation.Vdot_regressor = Vdot_regressor;
results.regressor_validation.absolute_error = regressor_validation_error;
results.cbf_regressor_validation.psi1dot_finite_difference = ...
    psi1dot_finite_difference;
results.cbf_regressor_validation.psi1dot_regressor = psi1dot_regressor;
results.cbf_regressor_validation.absolute_error = ...
    cbf_regressor_validation_error;
results.max_parametric_effect = max(abs(x_parametric-x_nominal),[],1);
results.max_actuation_effect = max(abs(x_actuation-x_nominal),[],1);
results.max_combined_effect = max(abs(x_combined-x_nominal),[],1);

fprintf('\nMini-Project 3 system validation\n');
fprintf('--------------------------------\n');
fprintf('||Delta_true||_2       = %.6f (bound %.6f)\n', ...
    norm(Delta_true,2), delta_max);
fprintf('min sigma(I_4+Delta)  = %.6f\n', ...
    min(svd(eye(4) + Delta_true)));
fprintf('Maximum model mismatch by state:\n');
state_names = {'phi','pitch','psi','p','q','r'};
for i = 1:6
    fprintf('  %-6s : %.3e\n', state_names{i}, max_equivalence_error(i));
end
fprintf('Overall maximum mismatch = %.3e\n\n', ...
    max(max_equivalence_error));

fprintf('Maximum open-loop trajectory change relative to nominal:\n');
fprintf('  State      Parametric      Actuation       Combined\n');
for i = 1:6
    fprintf('  %-6s   %11.3e   %11.3e   %11.3e\n',state_names{i}, ...
        results.max_parametric_effect(i), ...
        results.max_actuation_effect(i), ...
        results.max_combined_effect(i));
end
fprintf('\n');

fprintf('Task 1.2 estimator comparison:\n');
fprintf(['  Method                         final ||theta_tilde||', ...
    '   t_5%% [s]   int ||u||^2 dt   final lambda_min\n']);
for method_index = 1:numel(estimator_results)
    item = estimator_results(method_index);
    fprintf('  %-30s   %14.3e   %8.3f   %14.3e   %14.3e\n', ...
        item.label,item.parameter_error_norm(end),item.convergence_time_5pct, ...
        item.control_energy,item.lambda_min(end));
end
fprintf('\n');

fprintf('Task 1.3 ISSf-CBF safety-filter summary:\n');
fprintf('  HOCBF psi1dot check error: %.3e\n', ...
    cbf_regressor_validation_error);
fprintf('  Case       min h [deg]   max pitch [deg]   energy       infeasible\n');
for ic = 1:numel(task1_3)
    fprintf('  %-8s   %11.4f   %15.4f   %9.3e   %6d\n', ...
        task1_3(ic).mode,rad2deg(task1_3(ic).min_h), ...
        max(rad2deg(task1_3(ic).x(:,2))), ...
        task1_3(ic).control_energy,sum(task1_3(ic).exitflag <= 0));
end
fprintf('\n');

fprintf('Task 1.4 estimator comparison with eISS-CLF + ISSf-CBF:\n');
fprintf(['  Method                         final ||theta_tilde||', ...
    '   t_5%% [s]   lambda_min    min h [deg]   energy      infeasible\n']);
for method_index = 1:numel(task1_4)
    item = task1_4(method_index);
    fprintf('  %-30s   %14.3e   %8.3f   %10.3e   %11.4f   %9.3e   %6d\n', ...
        item.label,item.parameter_error_norm(end), ...
        item.convergence_time_5pct,item.lambda_min(end), ...
        rad2deg(item.min_h),item.control_energy,item.infeasible_count);
end
fprintf('\n');

fprintf('Task 2.1 robust-CLF summary:\n');
fprintf(['  Case                               V(tf)/V(0)', ...
    '   max actual r   max set r     energy      infeasible  fallback  best-effort\n']);
for ic = 1:numel(robust_clf)
    item = robust_clf(ic);
    fprintf(['  %-34s   %10.3e   %12.3e   %10.3e', ...
        '   %9.3e   %6d   %6d   %6d\n'], ...
        sprintf('case %d',ic),item.final_V_ratio, ...
        item.max_actual_decay_residual,item.max_exact_set_residual, ...
        item.control_energy,sum(item.exitflag <= 0), ...
        item.spectral_fallback_count,item.best_effort_count);
end
fprintf('\n');

fprintf('Task 2.2 robust-CBF safety-filter summary:\n');
fprintf(['  Case                         min h [deg]   min psi1', ...
    '   min actual m   min set m     energy      infeasible  best-effort\n']);
for ic = 1:numel(robust_cbf)
    item = robust_cbf(ic);
    fprintf(['  %-28s   %11.4f   %9.3e   %12.3e', ...
        '   %10.3e   %9.3e   %6d   %6d\n'], ...
        sprintf('case %d',ic),rad2deg(item.min_h),item.min_psi1, ...
        item.min_actual_margin,item.min_exact_set_margin, ...
        item.control_energy,sum(item.exitflag <= 0), ...
        item.best_effort_count);
end
fprintf('\n');

fprintf('Task 2.3 SMID uncertainty-reduction summary:\n');
fprintf(['  Case                         final bound   true ||Delta||', ...
    '   rank time   min h [deg]   control E   intervention E   CBF inflation\n']);
for ic = 1:numel(task2_3)
    item = task2_3(ic);
    fprintf(['  %-28s   %11.6f   %14.6f   %9.3f', ...
        '   %11.5f   %9.3e   %14.3e   %9.3e\n'], ...
        sprintf('case %d',ic),item.final_delta_bound, ...
        item.true_delta_norm,item.full_rank_time,rad2deg(item.min_h), ...
        item.control_energy,item.intervention_energy, ...
        item.safety_inflation_area);
end
fprintf(['  SMID final active-block error = %.3e, ', ...
    'minimum containment margin = %.3e\n\n'], ...
    online_smid.Delta_error(end),online_smid.min_containment_margin);

fprintf('Task 1.1 eISS-CLF summary:\n');
fprintf('  Regressor Vdot check error: %.3e\n',regressor_validation_error);
fprintf('  Case                       V(0)       V(tf)      max|u|    max QP residual\n');
for ic = 1:numel(eiss_clf)
    fprintf('  %-24s  %9.3e  %9.3e  %9.3e  %13.3e\n', ...
        sprintf('case %d',ic),eiss_clf(ic).V(1),eiss_clf(ic).V(end), ...
        max(abs(eiss_clf(ic).u(:))), ...
        max(eiss_clf(ic).constraint_residual));
end
fprintf('\n');

save(fullfile(project_dir, 'Validation_Results.mat'), 'results', 'params');
Plot_Results(results);


function dx = affine_rhs(x, u, params, theta, Delta)
    [~,~,~,dx] = Dynamics(x, params, theta, u, Delta);
end


function dx = original_no_damping_rhs(x, u, params)

    [f0,~,g] = Dynamics(x, params);
    p = x(4);
    q = x(5);
    r = x(6);

    omega_dot = [params.Gam1*p*q - params.Gam2*q*r;
                 params.Gam5*p*r - params.Gam6*(p^2-r^2);
                 params.Gam7*p*q - params.Gam1*q*r];

    dx = f0 + [zeros(3,1); omega_dot] + g*u;
end


function dx = eiss_closed_loop_rhs(x,theta_hat,params,x_ref,Delta)
    [u,~] = eISS_CLF_Controller(x,theta_hat,params,x_ref);
    [~,~,~,dx] = Dynamics(x,params,params.theta_true,u,Delta);
end


function dx = issf_closed_loop_rhs(x,mode,theta_hat,params,x_ref,Delta)
    [u_nom,~] = eISS_CLF_Controller(x,theta_hat,params,x_ref);
    if strcmp(mode,'none')
        u = u_nom;
    else
        [u,~] = ISSf_CBF_SafetyFilter(x,u_nom,theta_hat,params,mode);
    end
    [~,~,~,dx] = Dynamics(x,params,params.theta_true,u,Delta);
end


function dx = robust_clf_closed_loop_rhs( ...
    x,controller_delta,uncertainty_mode,theta_known,params,x_ref,Delta_fixed)
    [u,controller_info] = Robust_CLF_Controller( ...
        x,theta_known,params,x_ref,controller_delta);
    Delta = select_actuation_uncertainty( ...
        uncertainty_mode,controller_info.LgV,u, ...
        params.delta_max,Delta_fixed);
    [~,~,~,dx] = Dynamics(x,params,theta_known,u,Delta);
end


function Delta = select_actuation_uncertainty( ...
    uncertainty_mode,LgV,u,delta_max,Delta_fixed)

    if strcmp(uncertainty_mode,'zero')
        Delta = zeros(4,4);
    elseif strcmp(uncertainty_mode,'fixed')
        Delta = Delta_fixed;
    elseif strcmp(uncertainty_mode,'worst')
        row_norm = norm(LgV,2);
        input_norm = norm(u,2);
        if row_norm <= eps || input_norm <= eps
            Delta = zeros(4,4);
        else
            Delta = delta_max*(LgV.'/row_norm)*(u.'/input_norm);
        end
    else
        error('Main:UnknownUncertaintyMode', ...
            'Unknown uncertainty mode: %s',uncertainty_mode);
    end
end


function dx = robust_cbf_closed_loop_rhs( ...
    x,filter_mode,filter_delta,uncertainty_mode,theta_known,params, ...
    x_ref,Delta_fixed)

    [u_nom,~] = Robust_CLF_Controller(x,theta_known,params,x_ref,0);
    reg = Regressors(x,params,x_ref);
    if strcmp(filter_mode,'none')
        u = u_nom;
    else
        [u,~] = Robust_CBF_SafetyFilter( ...
            x,u_nom,theta_known,params,filter_delta);
    end
    Delta = select_cbf_actuation_uncertainty( ...
        uncertainty_mode,reg.Lg_psi1,u,params.delta_max,Delta_fixed);
    [~,~,~,dx] = Dynamics(x,params,theta_known,u,Delta);
end


function Delta = select_cbf_actuation_uncertainty( ...
    uncertainty_mode,Lg_psi1,u,delta_max,Delta_fixed)

    if strcmp(uncertainty_mode,'zero')
        Delta = zeros(4,4);
    elseif strcmp(uncertainty_mode,'fixed')
        Delta = Delta_fixed;
    elseif strcmp(uncertainty_mode,'worst')
        row_norm = norm(Lg_psi1,2);
        input_norm = norm(u,2);
        if row_norm <= eps || input_norm <= eps
            Delta = zeros(4,4);
        else
            Delta = -delta_max*(Lg_psi1.'/row_norm)*(u.'/input_norm);
        end
    else
        error('Main:UnknownCBFUncertaintyMode', ...
            'Unknown CBF uncertainty mode: %s',uncertainty_mode);
    end
end


function result = simulate_estimator_case( ...
    estimator_function,label,params,x0,x_ref)

    dt = params.estimator_dt;
    t = (0:dt:params.estimator_t_final).';
    n_samples = numel(t);
    n_theta = params.n_theta;

    simulation_params = params;
    simulation_params.use_quadprog = false;
    rng(params.estimator_rng_seed,'twister');

    x = zeros(n_samples,6);
    theta_hat = zeros(n_samples,n_theta);
    u = zeros(n_samples,4);
    V = zeros(n_samples,1);
    parameter_error_norm = zeros(n_samples,1);
    prediction_error_norm = zeros(n_samples,1);
    lambda_min = zeros(n_samples,1);
    stack_count = zeros(n_samples,1);
    forgetting_factor = ones(n_samples,1);
    covariance_trace = NaN(n_samples,1);
    qp_margin = zeros(n_samples,1);

    x(1,:) = x0.';
    theta_hat(1,:) = params.theta_hat0.';
    stack = [];
    estimator_state = [];

    for k = 1:n_samples
        xk = x(k,:).';
        theta_k = theta_hat(k,:).';
        [uk,controller_info] = eISS_CLF_Controller( ...
            xk,theta_k,simulation_params,x_ref);
        [f0,L,g,xdot] = Dynamics(xk,simulation_params, ...
            params.theta_true,uk,zeros(4));

        Y = xdot-f0-g*uk;
        if params.estimator_noise_std > 0
            Y = Y+params.estimator_noise_std*randn(size(Y));
        end
        [stack,stack_info] = HistoryStack(stack,L,Y,params);

        u(k,:) = uk.';
        V(k) = controller_info.V;
        parameter_error_norm(k) = norm(params.theta_true-theta_k);
        prediction_error_norm(k) = norm(Y-L*theta_k);
        lambda_min(k) = stack_info.lambda_min;
        stack_count(k) = stack_info.count;
        qp_margin(k) = -controller_info.constraint_residual;

        if k < n_samples
            [theta_next,estimator_state,estimator_info] = ...
                estimator_function(theta_k,L,Y,stack,estimator_state, ...
                params,dt);
            theta_hat(k+1,:) = theta_next.';
            forgetting_factor(k) = estimator_info.forgetting_factor;
            covariance_trace(k) = estimator_info.covariance_trace;
            x(k+1,:) = rk4_eiss_step(xk,theta_k,simulation_params,x_ref,dt).';
        elseif k > 1
            forgetting_factor(k) = forgetting_factor(k-1);
            covariance_trace(k) = covariance_trace(k-1);
        end
    end

    if any(~isfinite(x(:))) || any(~isfinite(theta_hat(:)))
        error('Main:EstimatorDivergence', ...
            'Non-finite state or parameter estimate in method: %s',label);
    end

    error_threshold = 0.05*parameter_error_norm(1);
    convergence_index = find(parameter_error_norm <= error_threshold,1,'first');
    if isempty(convergence_index)
        convergence_time = NaN;
    else
        convergence_time = t(convergence_index);
    end

    result.label = label;
    result.t = t;
    result.x = x;
    result.theta_hat = theta_hat;
    result.u = u;
    result.V = V;
    result.parameter_error_norm = parameter_error_norm;
    result.prediction_error_norm = prediction_error_norm;
    result.lambda_min = lambda_min;
    result.stack_count = stack_count;
    result.forgetting_factor = forgetting_factor;
    result.covariance_trace = covariance_trace;
    result.qp_margin = qp_margin;
    result.control_energy = trapz(t,sum(u.^2,2));
    result.convergence_time_5pct = convergence_time;
    result.final_stack = stack;
end


function result = simulate_estimator_issf_case( ...
    estimator_function,label,params,x0,x_ref)

    dt = params.estimator_dt;
    t = (0:dt:params.estimator_t_final).';
    n_samples = numel(t);
    n_theta = params.n_theta;

    simulation_params = params;
    simulation_params.use_quadprog = false;
    rng(params.estimator_rng_seed,'twister');

    x = zeros(n_samples,6);
    theta_hat = zeros(n_samples,n_theta);
    u_nom = zeros(n_samples,4);
    u = zeros(n_samples,4);
    h = zeros(n_samples,1);
    psi1 = zeros(n_samples,1);
    iota_h = zeros(n_samples,1);
    estimated_margin = zeros(n_samples,1);
    true_margin = zeros(n_samples,1);
    parameter_error_norm = zeros(n_samples,1);
    prediction_error_norm = zeros(n_samples,1);
    lambda_min = zeros(n_samples,1);
    stack_count = zeros(n_samples,1);
    forgetting_factor = ones(n_samples,1);
    covariance_trace = NaN(n_samples,1);
    clf_exitflag = ones(n_samples,1);
    cbf_exitflag = ones(n_samples,1);

    x(1,:) = x0.';
    theta_hat(1,:) = params.theta_hat0.';
    stack = [];
    estimator_state = [];

    for k = 1:n_samples
        xk = x(k,:).';
        theta_k = theta_hat(k,:).';

        [u_nom_k,clf_info] = eISS_CLF_Controller( ...
            xk,theta_k,simulation_params,x_ref);
        [uk,cbf_info] = ISSf_CBF_SafetyFilter( ...
            xk,u_nom_k,theta_k,simulation_params,'issf');
        [f0,L,g,xdot] = Dynamics(xk,simulation_params, ...
            params.theta_true,uk,zeros(4));

        Y = xdot-f0-g*uk;
        if params.estimator_noise_std > 0
            Y = Y+params.estimator_noise_std*randn(size(Y));
        end
        [stack,stack_info] = HistoryStack(stack,L,Y,params);

        reg_k = Regressors(xk,simulation_params,x_ref);
        u_nom(k,:) = u_nom_k.';
        u(k,:) = uk.';
        h(k) = reg_k.h;
        psi1(k) = reg_k.psi1;
        iota_h(k) = cbf_info.iota_h;
        estimated_margin(k) = cbf_info.estimated_margin;
        true_margin(k) = reg_k.Lf_psi1 ...
            + reg_k.phi_h.'*params.theta_true + reg_k.Lg_psi1*uk ...
            + simulation_params.alpha2*reg_k.psi1;
        parameter_error_norm(k) = norm(params.theta_true-theta_k);
        prediction_error_norm(k) = norm(Y-L*theta_k);
        lambda_min(k) = stack_info.lambda_min;
        stack_count(k) = stack_info.count;
        clf_exitflag(k) = clf_info.exitflag;
        cbf_exitflag(k) = cbf_info.exitflag;

        if k < n_samples
            [theta_next,estimator_state,estimator_info] = ...
                estimator_function(theta_k,L,Y,stack,estimator_state, ...
                params,dt);
            theta_hat(k+1,:) = theta_next.';
            forgetting_factor(k) = estimator_info.forgetting_factor;
            covariance_trace(k) = estimator_info.covariance_trace;
            x(k+1,:) = rk4_eiss_issf_step( ...
                xk,theta_k,simulation_params,x_ref,dt).';
        elseif k > 1
            forgetting_factor(k) = forgetting_factor(k-1);
            covariance_trace(k) = covariance_trace(k-1);
        end
    end

    if any(~isfinite(x(:))) || any(~isfinite(theta_hat(:)))
        error('Main:Task14Divergence', ...
            'Non-finite Task 1.4 result in method: %s',label);
    end

    error_threshold = 0.05*parameter_error_norm(1);
    convergence_index = find(parameter_error_norm <= error_threshold,1,'first');
    if isempty(convergence_index)
        convergence_time = NaN;
    else
        convergence_time = t(convergence_index);
    end

    result.label = label;
    result.t = t;
    result.x = x;
    result.state_norm = sqrt(sum(x.^2,2));
    result.theta_hat = theta_hat;
    result.u_nom = u_nom;
    result.u = u;
    result.intervention_norm = sqrt(sum((u-u_nom).^2,2));
    result.h = h;
    result.psi1 = psi1;
    result.iota_h = iota_h;
    result.estimated_margin = estimated_margin;
    result.true_margin = true_margin;
    result.parameter_error_norm = parameter_error_norm;
    result.prediction_error_norm = prediction_error_norm;
    result.lambda_min = lambda_min;
    result.stack_count = stack_count;
    result.forgetting_factor = forgetting_factor;
    result.covariance_trace = covariance_trace;
    result.clf_exitflag = clf_exitflag;
    result.cbf_exitflag = cbf_exitflag;
    result.control_energy = trapz(t,sum(u.^2,2));
    result.convergence_time_5pct = convergence_time;
    result.min_h = min(h);
    result.min_psi1 = min(psi1);
    result.min_true_margin = min(true_margin);
    result.min_qp_margin = min(estimated_margin);
    result.max_abs_input = max(abs(u(:)));
    result.clf_infeasible_count = sum(clf_exitflag <= 0);
    result.cbf_infeasible_count = sum(cbf_exitflag <= 0);
    result.infeasible_count = sum(clf_exitflag <= 0 | cbf_exitflag <= 0);
    result.final_stack = stack;

    if result.min_h < -1e-6
        warning('Main:Task14SafetyViolation', ...
            '%s violated h>=0 by %.3e rad.',label,-result.min_h);
    end
    if result.infeasible_count > 0
        warning('Main:Task14QPInfeasible', ...
            '%s had %d infeasible CLF/CBF samples.', ...
            label,result.infeasible_count);
    end
    if result.min_qp_margin < -1e-8
        warning('Main:Task14CBFResidual', ...
            '%s had a minimum implemented CBF margin of %.3e.', ...
            label,result.min_qp_margin);
    end
    if result.min_true_margin < -1e-6
        warning('Main:Task14TrueCBFViolation', ...
            '%s had a minimum true HOCBF margin of %.3e.', ...
            label,result.min_true_margin);
    end

    if isscalar(params.u_min)
        lower_limit = params.u_min*ones(1,size(u,2));
    else
        lower_limit = params.u_min(:).';
    end
    if isscalar(params.u_max)
        upper_limit = params.u_max*ones(1,size(u,2));
    else
        upper_limit = params.u_max(:).';
    end
    lower_violation = any(any(u < repmat(lower_limit,size(u,1),1)-1e-10));
    upper_violation = any(any(u > repmat(upper_limit,size(u,1),1)+1e-10));
    if lower_violation || upper_violation
        error('Main:Task14InputViolation', ...
            '%s exceeded the actuator limit.',label);
    end
end


function x_next = rk4_eiss_step(x,theta_hat,params,x_ref,dt)
    k1 = eiss_vector_field(x,theta_hat,params,x_ref);
    k2 = eiss_vector_field(x+0.5*dt*k1,theta_hat,params,x_ref);
    k3 = eiss_vector_field(x+0.5*dt*k2,theta_hat,params,x_ref);
    k4 = eiss_vector_field(x+dt*k3,theta_hat,params,x_ref);
    x_next = x+(dt/6)*(k1+2*k2+2*k3+k4);
end


function x_next = rk4_eiss_issf_step(x,theta_hat,params,x_ref,dt)
    k1 = eiss_issf_vector_field(x,theta_hat,params,x_ref);
    k2 = eiss_issf_vector_field(x+0.5*dt*k1,theta_hat,params,x_ref);
    k3 = eiss_issf_vector_field(x+0.5*dt*k2,theta_hat,params,x_ref);
    k4 = eiss_issf_vector_field(x+dt*k3,theta_hat,params,x_ref);
    x_next = x+(dt/6)*(k1+2*k2+2*k3+k4);
end


function dx = eiss_vector_field(x,theta_hat,params,x_ref)
    [u,~] = eISS_CLF_Controller(x,theta_hat,params,x_ref);
    [~,~,~,dx] = Dynamics(x,params,params.theta_true,u,zeros(4));
end


function dx = eiss_issf_vector_field(x,theta_hat,params,x_ref)
    [u_nom,~] = eISS_CLF_Controller(x,theta_hat,params,x_ref);
    [u,~] = ISSf_CBF_SafetyFilter(x,u_nom,theta_hat,params,'issf');
    [~,~,~,dx] = Dynamics(x,params,params.theta_true,u,zeros(4));
end
