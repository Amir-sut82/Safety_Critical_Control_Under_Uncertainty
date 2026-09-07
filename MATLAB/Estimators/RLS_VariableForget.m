function [theta_next,state,info] = RLS_VariableForget( ...
    theta_hat,L,Y,stack,state,params,dt)
% Error-dependent forgetting with covariance clipping

    state = initialize_state(state,params);
    prediction_error = Y-L*theta_hat;

    normalized_error = norm(prediction_error)^2/( ...
        params.vff_epsilon+norm(Y)^2);
    rho = params.vff_error_gain*normalized_error;
    lambda = 1-(1-params.vff_lambda_min)*rho/(1+rho);

    if norm(L,'fro') < params.vff_excitation_threshold
        lambda = 1;
    end
    lambda = min(1,max(params.vff_lambda_min,lambda));

    [theta_next,state.P] = rls_update(theta_hat,state.P,L,Y,lambda,params);
    [theta_next,state] = replay_history(theta_next,state,stack,params);
    state.P = clip_covariance(state.P,params.vff_cov_min,params.vff_cov_max);

    info.prediction_error_norm = norm(prediction_error);
    info.normalized_error = normalized_error;
    info.forgetting_factor = lambda;
    info.covariance_trace = trace(state.P);
end


function state = initialize_state(state,params)
    if isempty(state)
        state.P = params.rls_P0*eye(params.n_theta);
        state.replay_index = 1;
    end
end


function [theta,state] = replay_history(theta,state,stack,params)
    if stack.count == 0
        return;
    end
    index = min(state.replay_index,stack.count);
    [theta,state.P] = rls_update(theta,state.P,stack.L(:,:,index), ...
        stack.Y(:,index),1,params);
    state.replay_index = mod(index,stack.count)+1;
end


function [theta_new,P_new] = rls_update(theta,P,L,Y,lambda,params)
    S = lambda*eye(size(L,1))+L*P*L.';
    K = (P*L.')/S;
    theta_new = theta+K*(Y-L*theta);
    P_new = (P-K*L*P)/lambda;
    P_new = clip_covariance(P_new,params.rls_cov_min,params.rls_cov_max);
end


function P = clip_covariance(P,lower,upper)
    [vectors,values] = eig(0.5*(P+P.'));
    diagonal = min(max(real(diag(values)),lower),upper);
    P = vectors*diag(diagonal)*vectors.';
    P = 0.5*(P+P.');
end
