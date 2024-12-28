function X_hat_k = extract_multiple_target_state(w, m)
    X_hat_k = {};
    Jk = length(w);
    for i = 1:Jk
        disp(w{i})
        if w{i} > 0.5
            for j = 1:round(w{i})
                X_hat_k = [X_hat_k, m{i}];
            end
        end
    end
end