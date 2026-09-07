function params = System_Parameters()

    params.m     = 0.486;
    params.g_acc = 9.81;

    params.Ixx = 0.0043;
    params.Iyy = 0.0142;
    params.Izz = 0.0176;
    params.Ixz = 0.0001;

    params.I = [ params.Ixx, 0, -params.Ixz;
                 0, params.Iyy, 0;
                -params.Ixz, 0, params.Izz];

    params.dcg = 0.175;
    params.hcg = 0.085;
    params.b   = 3.13e-5;
    params.d   = 1.2e-6;
    params.d_over_b = params.d / params.b;

    params.B_alloc = [0, 0,          params.dcg,      params.d_over_b;
                      0, params.hcg, 0,               0;
                      0, 0,          params.d_over_b, -params.dcg];

    Ixx = params.Ixx;
    Iyy = params.Iyy;
    Izz = params.Izz;
    Ixz = params.Ixz;
    Gam = Ixx*Izz - Ixz^2;

    Gam1 = Ixz*(Ixx - Iyy + Izz)/Gam;
    Gam2 = (Izz*(Izz - Iyy) + Ixz^2)/Gam;
    Gam5 = (Izz - Ixx)/Iyy;
    Gam6 = Ixz/Iyy;
    Gam7 = (Ixx*(Ixx - Iyy) + Ixz^2)/Gam;

    params.Gam  = Gam;
    params.Gam1 = Gam1;
    params.Gam2 = Gam2;
    params.Gam5 = Gam5;
    params.Gam6 = Gam6;
    params.Gam7 = Gam7;

    params.c_d = 0.01;
    params.theta_true = [Gam1; Gam2; Gam5; Gam6; Gam7; params.c_d];
    params.n_theta = 6;
    params.theta_nodamp = [Gam1; Gam2; Gam5; Gam6; Gam7; 0];
    params.theta_hat0 = [0; 0.5; 0.5; 0; -0.3; 0];
    params.theta_bar = 1.5*norm(params.theta_true - params.theta_hat0);

    params.Gamma_adapt = 10*eye(params.n_theta);
    params.gamma_clf = 0.5;
    params.c_eiss = 0.5;
    params.kappa_iss = 10;
    params.gamma_cbf = 5;
    params.gamma_cl = 5;
    params.N_stack = 20;
    % Part 1 estimation
    params.estimator_dt = 0.01;
    params.estimator_t_final = 30;
    params.estimator_noise_std = 0;
    params.estimator_rng_seed = 7;
    params.stack_regressor_threshold = 1e-5;
    params.stack_replacement_tolerance = 1e-10;
    params.stack_lambda_threshold = 1e-4;
    params.Gamma_gradient = 10*eye(params.n_theta);
    params.gradient_history_gain = 5;
    params.rls_P0 = 100;
    params.rls_cov_min = 1e-8;
    params.rls_cov_max = 1e5;
    params.rls_lambda_constant = 0.995;
    params.vff_lambda_min = 0.97;
    params.vff_error_gain = 20;
    params.vff_epsilon = 1e-8;
    params.vff_excitation_threshold = 1e-4;
    params.vff_cov_min = 1e-6;
    params.vff_cov_max = 1e3;
    params.use_quadprog = true;

    params.theta_max = deg2rad(20);
    params.alpha1 = 5;
    params.alpha2 = 5;
    params.kappa_issf = 0.5;
    % Part 2 robust control
    params.delta_max = 0.25;
    params.robust_gamma = 0.5;
    params.robust_t_final = 15;
    params.robust_use_quadprog = false;
    params.robust_cbf_use_quadprog = false;
    params.robust_cbf_t_final = 15;
    params.robust_qp_tolerance = 1e-10;
    % SMID uses only the identifiable actuator channels
    params.smid_active_channels = 2:4;
    params.smid_dt = 0.01;
    params.smid_t_final = 15;
    params.smid_noise_bound = 1e-9;
    params.smid_numerical_margin = 1e-9;
    params.smid_min_input_norm = 1e-6;
    params.smid_sigma_threshold = 5e-3;
    params.smid_max_bound_decrease = 5e-4;
    params.smid_probe_amplitude = 3.5e-3;
    params.smid_probe_end = 4;
    params.smid_theta_ref = deg2rad(19.5);
    params.robust_u_min = -inf;
    params.robust_u_max = inf;

    params.u_min = -0.1;
    params.u_max = 0.1;

    params.Kp = 4*eye(3);
    params.Kd = 4*eye(3);
    F_mat = [zeros(3), eye(3); zeros(3), zeros(3)];
    G_mat = [zeros(3); eye(3)];
    K = [params.Kp, params.Kd];
    Acl = F_mat - G_mat*K;
    P = lyap(Acl', 0.1*eye(6));

    params.P = P;
    params.F_mat = F_mat;
    params.G_mat = G_mat;
    params.K = K;
    params.Acl = Acl;
    params.Q = 0.1*eye(6);

    params.p_penalty = 1e4;

    params.theta_ref = deg2rad(30);
    params.x_ref = [0; params.theta_ref; 0; 0; 0; 0];
    params.x_eq = zeros(6,1);

    params.t_final = 60;
    params.ode_opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, ...
        'MaxStep',0.01);
end
