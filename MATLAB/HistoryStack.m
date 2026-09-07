function [stack,info] = HistoryStack(stack,L,Y,params)
% Replace stored data only when lambda_min increases

    if isempty(stack)
        [n_y,n_theta] = size(L);
        stack.capacity = params.N_stack;
        stack.count = 0;
        stack.L = zeros(n_y,n_theta,stack.capacity);
        stack.Y = zeros(n_y,stack.capacity);
        stack.H = zeros(n_theta,n_theta);
        stack.b = zeros(n_theta,1);
        stack.lambda_min = 0;
    end

    if size(L,1) ~= size(stack.Y,1) || size(L,2) ~= size(stack.H,1)
        error('HistoryStack:RegressorDimension', ...
            'The new regression sample has incompatible dimensions.');
    end
    Y = Y(:);
    if numel(Y) ~= size(L,1)
        error('HistoryStack:OutputDimension', ...
            'Y must have one entry for each row of L.');
    end

    info.added = false;
    info.replaced_index = 0;
    info.previous_lambda_min = stack.lambda_min;

    if norm(L,'fro') < params.stack_regressor_threshold
        info.lambda_min = stack.lambda_min;
        info.count = stack.count;
        return;
    end

    if stack.count < stack.capacity
        stack.count = stack.count+1;
        index = stack.count;
        stack.L(:,:,index) = L;
        stack.Y(:,index) = Y;
        stack.H = stack.H+L.'*L;
        stack.b = stack.b+L.'*Y;
        stack.lambda_min = minimum_eigenvalue(stack.H);
        info.added = true;
        info.replaced_index = index;
    else
        current_lambda = stack.lambda_min;
        best_lambda = current_lambda;
        best_index = 0;
        best_H = stack.H;

        for index = 1:stack.capacity
            L_old = stack.L(:,:,index);
            H_candidate = stack.H-L_old.'*L_old+L.'*L;
            candidate_lambda = minimum_eigenvalue(H_candidate);
            if candidate_lambda > best_lambda+params.stack_replacement_tolerance
                best_lambda = candidate_lambda;
                best_index = index;
                best_H = H_candidate;
            end
        end

        if best_index > 0
            L_old = stack.L(:,:,best_index);
            Y_old = stack.Y(:,best_index);
            stack.L(:,:,best_index) = L;
            stack.Y(:,best_index) = Y;
            stack.H = best_H;
            stack.b = stack.b-L_old.'*Y_old+L.'*Y;
            stack.lambda_min = best_lambda;
            info.added = true;
            info.replaced_index = best_index;
        end
    end

    info.lambda_min = stack.lambda_min;
    info.count = stack.count;
    info.full_rank = stack.lambda_min >= params.stack_lambda_threshold;
end


function value = minimum_eigenvalue(H)
    eigenvalues = eig(0.5*(H+H.'));
    value = max(0,min(real(eigenvalues)));
end
