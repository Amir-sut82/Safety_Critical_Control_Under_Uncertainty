function [u,info] = ISSf_CBF_SafetyFilter(x,u_nom,theta_hat,params,mode)
% mode is naive or issf; the pitch barrier has relative degree two

    if nargin < 5 || isempty(mode)
        mode = 'issf';
    end
    mode = lower(char(mode));
    if ~ismember(mode,{'naive','issf'})
        error('ISSf_CBF_SafetyFilter:InvalidMode', ...
            'mode must be ''naive'' or ''issf''.');
    end

    u_nom = u_nom(:);
    theta_hat = theta_hat(:);
    if numel(u_nom) ~= 4
        error('ISSf_CBF_SafetyFilter:InputDimension', ...
            'u_nom must be a 4-by-1 vector.');
    end
    if numel(theta_hat) ~= params.n_theta
        error('ISSf_CBF_SafetyFilter:ParameterDimension', ...
            'theta_hat must contain %d parameters.',params.n_theta);
    end

    reg = Regressors(x,params,params.x_eq);
    if strcmp(mode,'issf')
        if params.kappa_issf <= 0
            error('ISSf_CBF_SafetyFilter:InvalidKappa', ...
                'params.kappa_issf must be strictly positive.');
        end
        iota_h = (reg.phi_h.'*reg.phi_h)/params.kappa_issf;
    else
        iota_h = 0;
    end
    A_cbf = -reg.Lg_psi1;
    b_cbf = reg.Lf_psi1 + reg.phi_h.'*theta_hat ...
        + params.alpha2*reg.psi1 - iota_h;

    n_u = numel(u_nom);
    lower_bound = expand_bound(params.u_min,n_u);
    upper_bound = expand_bound(params.u_max,n_u);

    used_quadprog = false;
    exitflag = 0;
    u = [];
    use_quadprog = isfield(params,'use_quadprog') && params.use_quadprog;
    if use_quadprog && exist('quadprog','file') == 2
        try
            options = optimoptions('quadprog','Display','off');
            [u,~,exitflag] = quadprog(2*eye(n_u),-2*u_nom, ...
                A_cbf,b_cbf,[],[],lower_bound,upper_bound,[],options);
            used_quadprog = ~isempty(u) && exitflag > 0;
        catch
            used_quadprog = false;
        end
    end

    if ~used_quadprog
        [u,exitflag] = project_to_halfspace_box( ...
            u_nom,A_cbf,b_cbf,lower_bound,upper_bound);
    end
    if A_cbf*u-b_cbf > 1e-12
        [u,exitflag] = project_to_halfspace_box( ...
            u_nom,A_cbf,b_cbf,lower_bound,upper_bound);
        used_quadprog = false;
    end

    estimated_margin = reg.Lf_psi1 + reg.phi_h.'*theta_hat ...
        + reg.Lg_psi1*u + params.alpha2*reg.psi1 - iota_h;

    info.mode = mode;
    info.h = reg.h;
    info.hdot = reg.hdot;
    info.psi1 = reg.psi1;
    info.Lf_psi1 = reg.Lf_psi1;
    info.phi_h = reg.phi_h;
    info.Lg_psi1 = reg.Lg_psi1;
    info.iota_h = iota_h;
    info.A_cbf = A_cbf;
    info.b_cbf = b_cbf;
    info.estimated_margin = estimated_margin;
    info.exitflag = exitflag;
    info.used_quadprog = used_quadprog;
end


function bound = expand_bound(value,n)
    if isscalar(value)
        bound = value*ones(n,1);
    else
        bound = value(:);
        if numel(bound) ~= n
            error('ISSf_CBF_SafetyFilter:BoundDimension', ...
                'Each input bound must be scalar or %d-by-1.',n);
        end
    end
end


function [u,exitflag] = project_to_halfspace_box( ...
    v,a,b,lower_bound,upper_bound)

    tolerance = 1e-10;
    u0 = min(max(v,lower_bound),upper_bound);
    if a*u0 <= b+tolerance
        u = u0;
        exitflag = 1;
        return;
    end

    u_min_linear = upper_bound;
    u_min_linear(a >= 0) = lower_bound(a >= 0);
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

    for k = 1:80
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
