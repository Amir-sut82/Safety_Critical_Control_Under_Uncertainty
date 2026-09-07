function [theta_next,state,info] = Gradient_Estimator( ...
    theta_hat,L,Y,stack,state,params,dt)
% Current-data gradient plus history-stack correction

    if nargin < 5 || isempty(state)
        state.initialized = true;
    end

    prediction_error = Y-L*theta_hat;
    history_gradient = stack.b-stack.H*theta_hat;
    theta_dot = params.Gamma_gradient*(L.'*prediction_error ...
        + params.gradient_history_gain*history_gradient);
    theta_next = theta_hat+dt*theta_dot;

    info.prediction_error_norm = norm(prediction_error);
    info.forgetting_factor = 1;
    info.covariance_trace = NaN;
    info.theta_dot = theta_dot;
end
