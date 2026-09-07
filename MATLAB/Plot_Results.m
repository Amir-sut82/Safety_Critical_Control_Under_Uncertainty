function figures = Plot_Results(results)

    required = {'t_equiv','x_original','x_affine','equivalence_error', ...
        'max_equivalence_error','t_uncertainty','x_nominal', ...
        'x_parametric','x_actuation','x_combined','U_commanded', ...
        'U_effective','output_dir'};
    for k = 1:numel(required)
        if ~isfield(results, required{k})
            error('Plot_Results:MissingField', ...
                'results.%s is required.', required{k});
        end
    end

    state_titles = {'Roll angle $\phi$','Pitch angle $\theta$', ...
        'Yaw angle $\psi$','Roll rate $p$','Pitch rate $q$', ...
        'Yaw rate $r$'};
    y_labels = {'$\phi~[\mathrm{deg}]$','$\theta~[\mathrm{deg}]$', ...
        '$\psi~[\mathrm{deg}]$','$p~[\mathrm{rad/s}]$', ...
        '$q~[\mathrm{rad/s}]$','$r~[\mathrm{rad/s}]$'};
    state_symbols = {'$\phi$','$\theta$','$\psi$','$p$','$q$','$r$'};
    angle_scale = [180/pi, 180/pi, 180/pi, 1, 1, 1];
    time_label = '$t~[\mathrm{s}]$';
    figures(1) = figure('Color','w','Name','Validation Model Equivalence', ...
        'Position',[80 60 1050 720]);
    for i = 1:6
        subplot(3,2,i);
        plot(results.t_equiv, results.x_original(:,i)*angle_scale(i), ...
            'b-', 'LineWidth',1.7); hold on;
        plot(results.t_equiv, results.x_affine(:,i)*angle_scale(i), ...
            'r--', 'LineWidth',1.25);
        latex_style(y_labels{i}, time_label, state_titles{i});
        if i == 1
            legend({'Original closed-form model', ...
                'Parameter-affine model $f_0+F\theta+gu$'}, ...
                'Interpreter','latex','Location','best','FontSize',9);
        end
    end
    latex_super_title(['Open-loop validation: original model versus ', ...
        '$\dot{x}=f_0(x)+F(x)\theta+g(x)u$']);
    figures(2) = figure('Color','w','Name','Validation Maximum Mismatch', ...
        'Position',[180 120 760 470]);

    max_error = max(results.max_equivalence_error(:), 1e-16);
    numerical_floor = 1e-16;
    for i = 1:6
        plot([i i], [numerical_floor max_error(i)], '-', ...
            'Color',[0.15 0.45 0.75], 'LineWidth',2.0); hold on;
    end
    scatter(1:6, max_error, 65, [0.05 0.30 0.65], 'filled');
    yline(1e-10,'k--','LineWidth',1.2);
    text(3.45,1.35e-10,'$10^{-10}$ validation tolerance', ...
        'Interpreter','latex','FontSize',10, ...
        'HorizontalAlignment','center');

    set(gca,'YScale','log','XTick',1:6,'XTickLabel',state_symbols);
    ylim([1e-16 1e-9]);
    xlim([0.5 6.5]);
    latex_style('$\max_t\!\left|x_{\mathrm{aff}}-x_{\mathrm{orig}}\right|$', ...
        '', 'Maximum absolute model mismatch');

    overall_error = max(max_error);
    text(0.03,0.91, ...
        ['$\max_{t,i}|e_i(t)| = ' sprintf('%.2e',overall_error) '$'], ...
        'Units','normalized','Interpreter','latex','FontSize',12, ...
        'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75]);
    figures(3) = figure('Color','w','Name','Validation Uncertainty Effects', ...
        'Position',[80 60 1050 720]);
    case_styles = {'k-','b--','m-.','r-'};
    case_names = {'Nominal: $\hat{\theta}(0),\;\Delta=0$', ...
        'Parametric: $\theta^*,\;\Delta=0$', ...
        'Actuation: $\hat{\theta}(0),\;\Delta=\Delta^*$', ...
        'Combined: $\theta^*,\;\Delta=\Delta^*$'};
    states = {results.x_nominal, results.x_parametric, ...
        results.x_actuation, results.x_combined};

    for i = 1:6
        subplot(3,2,i); hold on;
        for j = 1:4
            line_width = 1.15;
            if j == 4
                line_width = 1.7;
            end
            plot(results.t_uncertainty, ...
                states{j}(:,i)*angle_scale(i), case_styles{j}, ...
                'LineWidth',line_width);
        end
        latex_style(y_labels{i}, time_label, state_titles{i});
        if i == 1
            legend(case_names,'Interpreter','latex', ...
                'Location','best','FontSize',8.5);
        end
    end
    latex_super_title(['Open-loop response under parametric and ', ...
        'multiplicative actuation uncertainty']);
    active_channels = 2:4;
    figures(4) = figure('Color','w','Name','Validation Actuation Inputs', ...
        'Position',[80 160 1200 390]);
    for j = 1:numel(active_channels)
        i = active_channels(j);
        subplot(1,3,j);
        plot(results.t_uncertainty, results.U_commanded(:,i), ...
            'k-', 'LineWidth',1.55); hold on;
        plot(results.t_uncertainty, results.U_effective(:,i), ...
            'r--', 'LineWidth',1.55);
        latex_style(sprintf('$u_%d$',i), time_label, ...
            sprintf('Active input channel $%d$',i));
        if j == 1
            legend({'Commanded input $u$', ...
                'Effective input $(I_4+\Delta)u$'}, ...
                'Interpreter','latex','Location','best','FontSize',9);
        end
    end
    latex_super_title('Multiplicative actuation uncertainty: $u_{\mathrm{eff}}=(I_4+\Delta)u$');
    if isfield(results,'eiss_clf') && ~isempty(results.eiss_clf)
        eiss = results.eiss_clf;
        n_cases = numel(eiss);
        colors = lines(n_cases);

        figures(5) = figure('Color','w','Name','Task 1.1 eISS CLF States', ...
            'Position',[80 60 1050 720]);
        line_styles = {'-','--',':','-.'};
        for i = 1:6
            subplot(3,2,i); hold on;
            for ic = 1:n_cases
                plot(eiss(ic).t,eiss(ic).x(:,i)*angle_scale(i), ...
                    'Color',colors(ic,:), ...
                    'LineStyle',line_styles{min(ic,numel(line_styles))}, ...
                    'LineWidth',1.55);
            end
            latex_style(y_labels{i},time_label,state_titles{i});
            if i == 1
                legend({eiss.label},'Interpreter','latex', ...
                    'Location','best','FontSize',8.5);
            end
        end
        latex_super_title('Task 1.1: closed-loop response with the eISS-CLF-QP');

        figures(6) = figure('Color','w','Name','Task 1.1 eISS CLF Diagnostics', ...
            'Position',[100 80 1050 690]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cases
            semilogy(eiss(ic).t,max(eiss(ic).V,1e-16), ...
                'Color',colors(ic,:), ...
                'LineStyle',line_styles{min(ic,numel(line_styles))}, ...
                'LineWidth',1.55);
        end
        latex_style('$V(x)$',time_label,'eISS control Lyapunov function');
        set(gca,'YScale','log');
        legend({eiss.label},'Interpreter','latex','Location','best','FontSize',8.5);

        subplot(2,2,2); hold on;
        for ic = 1:n_cases
            control_norm = sqrt(sum(eiss(ic).u.^2,2));
            plot(eiss(ic).t,control_norm,'Color',colors(ic,:), ...
                'LineStyle',line_styles{min(ic,numel(line_styles))}, ...
                'LineWidth',1.55);
        end
        latex_style('$\|u(t)\|_2$',time_label,'Control effort');

        subplot(2,2,3); hold on;
        for ic = 1:n_cases
            semilogy(eiss(ic).t,max(eiss(ic).iota,1e-16), ...
                'Color',colors(ic,:), ...
                'LineStyle',line_styles{min(ic,numel(line_styles))}, ...
                'LineWidth',1.55);
        end
        latex_style('$\iota(x)=\|\varphi(x)\|_2^2/\kappa$', ...
            time_label,'ISS compensation term');

        subplot(2,2,4); hold on;
        for ic = 1:n_cases
            plot(eiss(ic).t,eiss(ic).constraint_residual, ...
                'Color',colors(ic,:), ...
                'LineStyle',line_styles{min(ic,numel(line_styles))}, ...
                'LineWidth',1.35);
        end
        yline(0,'k--','LineWidth',1.1);
        latex_style(['$\omega+\varphi^\top\hat{\theta}+L_gVu', ...
            '+cV-\iota$'],time_label,'QP constraint residual ($\leq0$)');

        latex_super_title('Task 1.1: eISS-CLF-QP diagnostics');
    end
    if isfield(results,'estimator_comparison') ...
            && ~isempty(results.estimator_comparison)
        comparison = results.estimator_comparison;
        n_methods = numel(comparison);
        method_colors = lines(n_methods);
        method_styles = {'-','--','-.',':'};
        method_labels = {comparison.label};

        figures(7) = figure('Color','w','Name', ...
            'Task 1.2 Estimator Comparison','Position',[90 60 1080 720]);

        subplot(2,2,1); hold on;
        for method = 1:n_methods
            semilogy(comparison(method).t, ...
                max(comparison(method).parameter_error_norm,1e-12), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.6);
        end
        latex_style('$\|\widetilde{\theta}(t)\|_2$',time_label, ...
            'Parameter-estimation error');
        legend(method_labels,'Interpreter','latex','Location','best','FontSize',8);

        subplot(2,2,2); hold on;
        for method = 1:n_methods
            plot(comparison(method).t,comparison(method).lambda_min, ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        latex_style('$\lambda_{\min}(\mathcal{H})$',time_label, ...
            'History-stack information level');

        subplot(2,2,3); hold on;
        for method = 1:n_methods
            state_norm = sqrt(sum(comparison(method).x.^2,2));
            semilogy(comparison(method).t,max(state_norm,1e-12), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        latex_style('$\|x(t)\|_2$',time_label,'Closed-loop state convergence');

        subplot(2,2,4); hold on;
        for method = 1:n_methods
            control_norm = sqrt(sum(comparison(method).u.^2,2));
            plot(comparison(method).t,control_norm, ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        latex_style('$\|u(t)\|_2$',time_label,'Instantaneous control effort');
        latex_super_title('Task 1.2: comparison of history-stack estimators');

        figures(8) = figure('Color','w','Name', ...
            'Task 1.2 Parameter Estimates','Position',[70 50 1120 760]);
        parameter_labels = {'$\Gamma_1$','$\Gamma_2$','$\Gamma_5$', ...
            '$\Gamma_6$','$\Gamma_7$','$c_d$'};
        for parameter = 1:6
            subplot(3,2,parameter); hold on;
            for method = 1:n_methods
                plot(comparison(method).t, ...
                    comparison(method).theta_hat(:,parameter), ...
                    'Color',method_colors(method,:), ...
                    'LineStyle',method_styles{method},'LineWidth',1.4);
            end
            yline(results.theta_true(parameter),'k--','LineWidth',1.25);
            latex_style(['$\widehat{\theta}_{' num2str(parameter) '}$'], ...
                time_label,parameter_labels{parameter});
            if parameter == 1
                legend([method_labels,{'True value'}], ...
                    'Interpreter','latex','Location','best','FontSize',7.5);
            end
        end
        latex_super_title('Task 1.2: estimated parameters versus true values');

        figures(9) = figure('Color','w','Name', ...
            'Task 1.2 Forgetting and Margins','Position',[90 60 1080 720]);

        subplot(2,2,1); hold on;
        for method = 1:n_methods
            plot(comparison(method).t,comparison(method).forgetting_factor, ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.5);
        end
        ylim([0.96 1.002]);
        latex_style('$\lambda_k$',time_label,'Forgetting-factor evolution');
        legend(method_labels,'Interpreter','latex','Location','best','FontSize',8);

        subplot(2,2,2); hold on;
        for method = 1:n_methods
            if any(isfinite(comparison(method).covariance_trace))
                semilogy(comparison(method).t, ...
                    max(comparison(method).covariance_trace,1e-12), ...
                    'Color',method_colors(method,:), ...
                    'LineStyle',method_styles{method},'LineWidth',1.5);
            end
        end
        latex_style('$\mathrm{tr}(P_k)$',time_label,'RLS covariance trace');

        subplot(2,2,3); hold on;
        for method = 1:n_methods
            semilogy(comparison(method).t, ...
                max(comparison(method).qp_margin,1e-16), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.5);
        end
        latex_style('$-r_{\mathrm{CLF}}(t)$',time_label, ...
            'eISS-CLF QP feasibility margin');

        subplot(2,2,4);
        energies = arrayfun(@(item)item.control_energy,comparison);
        bar(energies,0.65,'FaceColor',[0.25 0.55 0.80]);
        set(gca,'XTick',1:n_methods, ...
            'XTickLabel',{'Gradient','RLS','RLS-CF','RLS-VFF'});
        latex_style('$\int_0^T\|u(t)\|_2^2\,dt$','', ...
            'Integrated control effort');
        latex_super_title('Task 1.2: forgetting, covariance, and stability margins');
    end
    if isfield(results,'task1_3') && ~isempty(results.task1_3)
        cbf = results.task1_3;
        n_cbf = numel(cbf);
        cbf_colors = [0.20 0.20 0.20; 0.85 0.33 0.10; 0.00 0.45 0.74];
        cbf_styles = {'--','-.','-'};
        first_cbf_figure = numel(figures)+1;

        figures(first_cbf_figure) = figure('Color','w', ...
            'Name','Task 1.3 ISSf CBF Comparison', ...
            'Position',[80 60 1050 700]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cbf
            plot(cbf(ic).t,rad2deg(cbf(ic).x(:,2)),cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.65);
        end
        yline(20,'r--','LineWidth',1.25);
        yline(30,'k:','LineWidth',1.15);
        latex_style('$\theta~[\mathrm{deg}]$',time_label, ...
            'Pitch response');
        legend([{cbf.label},{'$\theta_{\max}=20^\circ$', ...
            '$\theta_{\mathrm{des}}=30^\circ$'}], ...
            'Interpreter','latex','Location','best','FontSize',8.2);

        subplot(2,2,2); hold on;
        for ic = 2:n_cbf
            plot(cbf(ic).t,rad2deg(cbf(ic).h),cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.65);
        end
        yline(0,'r--','LineWidth',1.2);
        latex_style('$h(x)=\theta_{\max}-\theta~[\mathrm{deg}]$', ...
            time_label,'Filtered safety margin ($h\geq0$)');
        legend({cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',8.5);

        subplot(2,2,3); hold on;
        for ic = 1:n_cbf
            control_norm = sqrt(sum(cbf(ic).u.^2,2));
            plot(cbf(ic).t,control_norm,cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.55);
        end
        latex_style('$\|u(t)\|_2$',time_label,'Applied control effort');

        subplot(2,2,4); hold on;
        for ic = 2:n_cbf
            plot(cbf(ic).t,cbf(ic).intervention_norm,cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.55);
        end
        latex_style('$\|u-u_{\mathrm{nom}}\|_2$',time_label, ...
            'Safety-filter intervention');
        legend({cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',8.5);

        latex_super_title(['Task 1.3: certainty-equivalence CBF versus ', ...
            'ISSf-CBF safety filter']);

        figures(first_cbf_figure+1) = figure('Color','w', ...
            'Name','Task 1.3 ISSf CBF Diagnostics', ...
            'Position',[100 80 1050 700]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cbf
            plot(cbf(ic).t,cbf(ic).psi1,cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.55);
        end
        yline(0,'r--','LineWidth',1.15);
        latex_style('$\psi_1=\dot h+\alpha_1h$',time_label, ...
            'High-order barrier state');

        subplot(2,2,2); hold on;
        for ic = 2:n_cbf
            semilogy(cbf(ic).t,max(cbf(ic).iota_h,1e-16), ...
                cbf_styles{ic},'Color',cbf_colors(ic,:), ...
                'LineWidth',1.55);
        end
        latex_style('$\iota_h(x)=\|\varphi_h(x)\|_2^2/\kappa_h$', ...
            time_label,'ISSf compensation term');
        legend({cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',8.5);

        subplot(2,2,3); hold on;
        for ic = 1:n_cbf
            plot(cbf(ic).t,cbf(ic).true_margin,cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.45);
        end
        yline(0,'r--','LineWidth',1.15);
        latex_style(['$L_f\psi_1+\varphi_h^\top\theta^*+', ...
            'L_g\psi_1u+\alpha_2\psi_1$'],time_label, ...
            'True HOCBF condition');

        subplot(2,2,4); hold on;
        for ic = 2:n_cbf
            plot(cbf(ic).t,cbf(ic).estimated_margin,cbf_styles{ic}, ...
                'Color',cbf_colors(ic,:),'LineWidth',1.45);
        end
        yline(0,'k--','LineWidth',1.1);
        latex_style('QP feasibility margin',time_label, ...
            'Implemented CBF constraint ($\geq0$)');

        latex_super_title('Task 1.3: ISSf-HOCBF diagnostic quantities');
    end
    if isfield(results,'task1_4') && ~isempty(results.task1_4)
        task14 = results.task1_4;
        n_methods = numel(task14);
        method_colors = lines(n_methods);
        method_styles = {'-','--','-.',':'};
        method_labels = {task14.label};
        short_labels = {'Gradient','RLS','RLS-CF','RLS-VFF'};
        first_task14_figure = numel(figures)+1;

        figures(first_task14_figure) = figure('Color','w','Name', ...
            'Task 1.4 Safe Closed Loop Estimator Comparison', ...
            'Position',[80 50 1100 740]);

        subplot(2,2,1); hold on;
        for method = 1:n_methods
            semilogy(task14(method).t, ...
                max(task14(method).parameter_error_norm,1e-12), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.65);
        end
        latex_style('$\|\widetilde{\theta}(t)\|_2$',time_label, ...
            'Parameter-estimation error');
        legend(method_labels,'Interpreter','latex', ...
            'Location','best','FontSize',8);

        subplot(2,2,2); hold on;
        for method = 1:n_methods
            plot(task14(method).t,task14(method).lambda_min, ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        latex_style('$\lambda_{\min}(\mathcal{H}(t))$',time_label, ...
            'History-stack excitation level');

        subplot(2,2,3); hold on;
        for method = 1:n_methods
            semilogy(task14(method).t,max(task14(method).state_norm,1e-12), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        latex_style('$\|x(t)\|_2$',time_label,'State convergence');

        subplot(2,2,4); hold on;
        for method = 1:n_methods
            plot(task14(method).t,rad2deg(task14(method).h), ...
                'Color',method_colors(method,:), ...
                'LineStyle',method_styles{method},'LineWidth',1.55);
        end
        yline(0,'r--','LineWidth',1.2);
        latex_style('$h(x)=\theta_{\max}-\theta~[\mathrm{deg}]$', ...
            time_label,'Safety-margin evolution');

        latex_super_title(['Task 1.4: four estimators with ', ...
            'eISS-CLF + ISSf-CBF']);

        figures(first_task14_figure+1) = figure('Color','w','Name', ...
            'Task 1.4 Estimator Tradeoff Metrics', ...
            'Position',[120 170 1050 430]);

        subplot(1,2,1);
        energies = arrayfun(@(item)item.control_energy,task14);
        bar(energies,0.65,'FaceColor',[0.25 0.55 0.80]);
        set(gca,'XTick',1:n_methods,'XTickLabel',short_labels);
        latex_style('$\int_0^T\|u(t)\|_2^2\,dt$','', ...
            'Integrated control effort');

        subplot(1,2,2);
        minimum_h_deg = rad2deg(arrayfun(@(item)item.min_h,task14));
        bar(minimum_h_deg,0.65,'FaceColor',[0.25 0.68 0.45]); hold on;
        yline(0,'r--','LineWidth',1.2);
        set(gca,'XTick',1:n_methods,'XTickLabel',short_labels);
        latex_style('$\min_t h(x(t))~[\mathrm{deg}]$','', ...
            'Minimum safety margin');

        latex_super_title(['Task 1.4: control-effort and ', ...
            'safety trade-off']);
    end
    if isfield(results,'robust_clf') && ~isempty(results.robust_clf)
        robust = results.robust_clf;
        n_cases = numel(robust);
        robust_colors = lines(n_cases);
        robust_styles = {'-','--','-.',':','-'};
        robust_labels = {robust.label};
        short_labels = {'Nom-0','Nom-WC','Rob-0','Rob-Fixed','Rob-WC'};
        first_robust_figure = numel(figures)+1;

        figures(first_robust_figure) = figure('Color','w','Name', ...
            'Task 2.1 Robust CLF Closed Loop', ...
            'Position',[80 50 1100 740]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cases
            semilogy(robust(ic).t,max(robust(ic).state_norm,1e-14), ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.55);
        end
        latex_style('$\|x(t)\|_2$',time_label,'State convergence');
        legend(robust_labels,'Interpreter','latex', ...
            'Location','best','FontSize',7.5);

        subplot(2,2,2); hold on;
        for ic = 1:n_cases
            semilogy(robust(ic).t,max(robust(ic).V,1e-16), ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.55);
        end
        latex_style('$V(x(t))$',time_label,'Control Lyapunov function');

        subplot(2,2,3); hold on;
        for ic = 1:n_cases
            control_norm = sqrt(sum(robust(ic).u.^2,2));
            plot(robust(ic).t,control_norm, ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.5);
        end
        latex_style('$\|u(t)\|_2$',time_label,'Control effort');

        subplot(2,2,4); hold on;
        for ic = 1:n_cases
            plot(robust(ic).t,robust(ic).actual_decay_residual, ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.4);
        end
        yline(0,'k--','LineWidth',1.15);
        latex_style('$\dot V+\gamma V$',time_label, ...
            'Actual CLF decay residual ($\leq0$)');

        latex_super_title(['Task 2.1: robust CLF under ', ...
            'multiplicative actuation uncertainty']);

        figures(first_robust_figure+1) = figure('Color','w','Name', ...
            'Task 2.1 Robust CLF Certificate', ...
            'Position',[95 55 1100 720]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cases
            plot(robust(ic).t,robust(ic).exact_set_residual, ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.4);
        end
        yline(0,'k--','LineWidth',1.15);
        latex_style(['$r_{2}=\dot V_{\mathrm{nom}}+', ...
            '\delta_{\max}\|L_gV\|_2\|u\|_2+\gamma V$'], ...
            time_label,'Exact spectral-ball certificate');
        legend(robust_labels,'Interpreter','latex', ...
            'Location','best','FontSize',7.5);

        subplot(2,2,2); hold on;
        for ic = 1:n_cases
            plot(robust(ic).t,robust(ic).dual_set_residual, ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.4);
        end
        yline(0,'k--','LineWidth',1.15);
        latex_style(['$r_{1}=\dot V_{\mathrm{nom}}+', ...
            '\delta_{\max}\|L_gV\|_2\|u\|_1+\gamma V$'], ...
            time_label,'Conservative dual-QP certificate');

        subplot(2,2,3); hold on;
        for ic = 1:n_cases
            plot(robust(ic).t,robust(ic).Delta_norm, ...
                'Color',robust_colors(ic,:), ...
                'LineStyle',robust_styles{ic},'LineWidth',1.45);
        end
        yline(results.delta_max,'r--','LineWidth',1.2);
        latex_style('$\|\Delta(t)\|_2$',time_label, ...
            'Applied uncertainty norm');

        subplot(2,2,4);
        robust_energies = arrayfun(@(item)item.control_energy,robust);
        bar(robust_energies,0.65,'FaceColor',[0.30 0.58 0.78]);
        set(gca,'XTick',1:n_cases,'XTickLabel',short_labels);
        latex_style('$\int_0^T\|u(t)\|_2^2\,dt$','', ...
            'Integrated control effort');

        latex_super_title(['Task 2.1: dual-QP robustness ', ...
            'certificate and uncertainty verification']);
    end
    if isfield(results,'robust_cbf') && ~isempty(results.robust_cbf)
        robust_cbf = results.robust_cbf;
        n_cases = numel(robust_cbf);
        cbf_colors = lines(n_cases);
        cbf_styles = {'--','-.','-','--',':'};
        cbf_labels = {robust_cbf.label};
        short_labels = {'Nom-WC','CBF-WC','Rob-0','Rob-Fixed','Rob-WC'};
        first_robust_cbf_figure = numel(figures)+1;

        figures(first_robust_cbf_figure) = figure('Color','w','Name', ...
            'Task 2.2 Robust CBF Comparison', ...
            'Position',[80 50 1100 740]);

        subplot(2,2,1); hold on;
        for ic = 1:n_cases
            plot(robust_cbf(ic).t,rad2deg(robust_cbf(ic).x(:,2)), ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.55);
        end
        yline(20,'r--','LineWidth',1.2);
        yline(30,'k:','LineWidth',1.1);
        latex_style('$\theta~[\mathrm{deg}]$',time_label,'Pitch response');
        legend([cbf_labels,{'$\theta_{\max}=20^\circ$', ...
            '$\theta_{\mathrm{des}}=30^\circ$'}], ...
            'Interpreter','latex','Location','best','FontSize',7.2);

        subplot(2,2,2); hold on;
        for ic = 2:n_cases
            plot(robust_cbf(ic).t,rad2deg(robust_cbf(ic).h), ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.55);
        end
        yline(0,'r--','LineWidth',1.2);
        latex_style('$h=\theta_{\max}-\theta~[\mathrm{deg}]$', ...
            time_label,'Filtered safety margin ($h\geq0$)');
        legend({robust_cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',7.5);

        subplot(2,2,3); hold on;
        for ic = 1:n_cases
            control_norm = sqrt(sum(robust_cbf(ic).u.^2,2));
            plot(robust_cbf(ic).t,control_norm, ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.5);
        end
        latex_style('$\|u(t)\|_2$',time_label,'Applied control effort');

        subplot(2,2,4); hold on;
        for ic = 2:n_cases
            plot(robust_cbf(ic).t,robust_cbf(ic).intervention_norm, ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.5);
        end
        latex_style('$\|u-u_{\mathrm{nom}}\|_2$',time_label, ...
            'Safety-filter intervention');
        legend({robust_cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',7.5);

        latex_super_title(['Task 2.2: robust HOCBF safety filtering ', ...
            'under actuation uncertainty']);

        figures(first_robust_cbf_figure+1) = figure('Color','w','Name', ...
            'Task 2.2 Robust CBF Certificate', ...
            'Position',[95 55 1100 720]);

        subplot(2,2,1); hold on;
        for ic = 2:n_cases
            plot(robust_cbf(ic).t,robust_cbf(ic).actual_margin, ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.45);
        end
        yline(0,'k--','LineWidth',1.1);
        latex_style('$\dot\psi_1+\alpha_2\psi_1$',time_label, ...
            'Actual HOCBF margin ($\geq0$)');
        legend({robust_cbf(2:end).label},'Interpreter','latex', ...
            'Location','best','FontSize',7.5);

        subplot(2,2,2); hold on;
        for ic = 2:n_cases
            plot(robust_cbf(ic).t,robust_cbf(ic).exact_set_margin, ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.45);
        end
        yline(0,'k--','LineWidth',1.1);
        latex_style(['$m_2=\dot\psi_{1,\mathrm{nom}}-', ...
            '\delta_{\max}\|L_g\psi_1\|_2\|u\|_2+', ...
            '\alpha_2\psi_1$'],time_label, ...
            'Exact spectral-ball certificate');

        subplot(2,2,3); hold on;
        for ic = 1:n_cases
            plot(robust_cbf(ic).t,robust_cbf(ic).Delta_norm, ...
                'Color',cbf_colors(ic,:), ...
                'LineStyle',cbf_styles{ic},'LineWidth',1.45);
        end
        yline(results.delta_max,'r--','LineWidth',1.2);
        latex_style('$\|\Delta(t)\|_2$',time_label, ...
            'Applied uncertainty norm');

        subplot(2,2,4);
        minimum_h_deg = rad2deg(arrayfun( ...
            @(item)item.min_h,robust_cbf(2:end)));
        bar(minimum_h_deg,0.65,'FaceColor',[0.25 0.68 0.45]); hold on;
        yline(0,'r--','LineWidth',1.2);
        set(gca,'XTick',1:(n_cases-1),'XTickLabel',short_labels(2:end));
        latex_style('$\min_t h(x(t))~[\mathrm{deg}]$','', ...
            'Minimum filtered safety margin');

        latex_super_title(['Task 2.2: robust-CBF certificate and ', ...
            'worst-case safety verification']);
    end
    if isfield(results,'task2_3') && ~isempty(results.task2_3)
        smid = results.task2_3;
        n_smid = numel(smid);
        smid_colors = [0.20 0.45 0.75; 0.85 0.33 0.10];
        smid_styles = {'--','-'};
        smid_labels = {smid.label};
        online = smid(2);
        first_smid_figure = numel(figures)+1;

        figures(first_smid_figure) = figure('Color','w','Name', ...
            'Task 2.3 SMID Set Contraction', ...
            'Position',[105 55 1100 720]);

        subplot(2,2,1); hold on;
        plot(online.t,online.delta_used,'Color',smid_colors(2,:), ...
            'LineWidth',1.8);
        yline(results.delta_max,'k--','LineWidth',1.15);
        yline(online.true_delta_norm,'r:','LineWidth',1.45);
        latex_style('$\hat\delta(t)$',time_label, ...
            'Certified uncertainty-bound contraction');
        legend({'Online SMID bound','$\delta_{\max}$', ...
            '$\|\Delta^*\|_2$'},'Interpreter','latex', ...
            'Location','best','FontSize',8.5);

        subplot(2,2,2); hold on;
        semilogy(online.t,max(online.lambda_min,1e-16), ...
            'Color',[0.25 0.60 0.35],'LineWidth',1.65);
        yline(1e-16,'k:','LineWidth',0.9);
        latex_style('$\lambda_{\min}(UU^\top)$',time_label, ...
            'SMID data excitation level');

        subplot(2,2,3); hold on;
        semilogy(online.t,max(online.Delta_error,1e-16), ...
            'Color',[0.55 0.25 0.65],'LineWidth',1.65);
        latex_style('$\|\widehat\Delta_a-\Delta_a^*\|_2$', ...
            time_label,'Active-block identification error');

        subplot(2,2,4); hold on;
        plot(online.t,online.exact_set_cbf_margin, ...
            'Color',[0.10 0.45 0.75],'LineWidth',1.55);
        plot(online.t,online.actual_cbf_margin,'--', ...
            'Color',[0.85 0.33 0.10],'LineWidth',1.45);
        yline(0,'k:','LineWidth',1.05);
        latex_style('$m_h(t)$',time_label, ...
            'Robust HOCBF verification ($m_h\geq0$)');
        legend({'Certified set margin','Actual margin'}, ...
            'Interpreter','latex','Location','best','FontSize',8.5);

        latex_super_title(['Task 2.3: set-membership identification ', ...
            'and certified uncertainty reduction']);

        figures(first_smid_figure+1) = figure('Color','w','Name', ...
            'Task 2.3 SMID Robust Control Tradeoff', ...
            'Position',[120 60 1100 720]);

        subplot(2,2,1); hold on;
        for ic = 1:n_smid
            plot(smid(ic).t,rad2deg(smid(ic).h), ...
                'Color',smid_colors(ic,:), ...
                'LineStyle',smid_styles{ic},'LineWidth',1.6);
        end
        yline(0,'r:','LineWidth',1.1);
        latex_style('$h=\theta_{\max}-\theta~[\mathrm{deg}]$', ...
            time_label,'Safety-margin evolution');
        legend(smid_labels,'Interpreter','latex','Location','best', ...
            'FontSize',8.2);

        subplot(2,2,2); hold on;
        for ic = 1:n_smid
            cumulative_energy = cumtrapz( ...
                smid(ic).t,sum(smid(ic).u.^2,2));
            plot(smid(ic).t,cumulative_energy, ...
                'Color',smid_colors(ic,:), ...
                'LineStyle',smid_styles{ic},'LineWidth',1.6);
        end
        latex_style('$\int_0^t\|u(\tau)\|_2^2d\tau$', ...
            time_label,'Cumulative control effort');

        subplot(2,2,3); hold on;
        for ic = 1:n_smid
            plot(smid(ic).t,smid(ic).filter_intervention, ...
                'Color',smid_colors(ic,:), ...
                'LineStyle',smid_styles{ic},'LineWidth',1.55);
        end
        latex_style('$\|u-u_{\mathrm{nom}}\|_2$',time_label, ...
            'Robust safety-filter intervention');

        subplot(2,2,4);
        fixed_metrics = [smid(1).control_energy, ...
            smid(1).intervention_energy,smid(1).safety_inflation_area];
        smid_metrics = [smid(2).control_energy, ...
            smid(2).intervention_energy,smid(2).safety_inflation_area];
        normalized_metrics = [ones(1,3); ...
            smid_metrics./max(fixed_metrics,eps)];
        bar(normalized_metrics.');
        yline(1,'k:','LineWidth',1.0);
        set(gca,'XTick',1:3,'XTickLabel', ...
            {'Control energy','Intervention','CBF inflation'});
        latex_style('Ratio to fixed-set design','', ...
            'Conservatism reduction metrics');
        legend({'Fixed bound','Online SMID'},'Interpreter','latex', ...
            'Location','best','FontSize',8.5);

        latex_super_title(['Task 2.3: robust-control performance ', ...
            'before and after online set contraction']);
    end
    export_dir = fullfile(results.output_dir,'plots_eps');
    if ~exist(export_dir,'dir')
        mkdir(export_dir);
    end

    drawnow;
    all_figures = findall(0,'Type','figure');
    all_figures = flipud(all_figures(:));
    used_names = {};
    for k = 1:numel(all_figures)
        fig = all_figures(k);
        if isempty(fig.Name)
            base_name = sprintf('Figure_%02d',fig.Number);
        else
            base_name = regexprep(fig.Name,'[^a-zA-Z0-9]+','_');
            base_name = regexprep(base_name,'^_+|_+$','');
        end
        if isempty(base_name)
            base_name = sprintf('Figure_%02d',fig.Number);
        end
        candidate_name = base_name;
        duplicate_index = 2;
        while any(strcmp(used_names,candidate_name))
            candidate_name = sprintf('%s_%02d',base_name,duplicate_index);
            duplicate_index = duplicate_index+1;
        end
        used_names{end+1} = candidate_name;

        eps_path = fullfile(export_dir,[candidate_name '.eps']);
        png_path = fullfile(export_dir,[candidate_name '.png']);
        try
            exportgraphics(fig,eps_path, ...
                'ContentType','vector','BackgroundColor','w');
            exportgraphics(fig,png_path, ...
                'Resolution',300,'BackgroundColor','w');
        catch
            print(fig,fullfile(export_dir,candidate_name), ...
                '-depsc','-r300');
            print(fig,png_path,'-dpng','-r300');
        end
    end
    fprintf('All %d open figures exported as EPS and PNG to: %s\n', ...
        numel(all_figures),export_dir);
end


function latex_style(y_label, x_label, title_text)
    grid on;
    grid minor;
    box on;
    set(gca,'FontSize',11,'LineWidth',1.1, ...
        'TickLabelInterpreter','latex');
    ylabel(y_label,'Interpreter','latex','FontSize',12);
    if ~isempty(x_label)
        xlabel(x_label,'Interpreter','latex','FontSize',12);
    end
    if nargin >= 3 && ~isempty(title_text)
        title(title_text,'Interpreter','latex','FontSize',13, ...
            'FontWeight','bold');
    end
end


function latex_super_title(title_text)
    if exist('sgtitle','file') == 2
        sgtitle(title_text,'Interpreter','latex','FontSize',15, ...
            'FontWeight','bold');
    else
        annotation('textbox',[0,0.955,1,0.04], ...
            'String',title_text,'Interpreter','latex', ...
            'EdgeColor','none','HorizontalAlignment','center', ...
            'FontSize',15,'FontWeight','bold');
    end
end
