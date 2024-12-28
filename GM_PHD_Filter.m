function estimated_states = GM_PHD_Filter(detections,time_steps,P_d,P_s,beta_FA,A,G,C,Q,R,V)

    estimated_states = cell(1,length(time_steps));

    birth_means = {};
    birth_covariances = {};
    birth_weights = {};

    birth_means{1} = [500,500,0,0]';
    birth_covariances{1} = diag([100,100,50,50]);
    birth_weights{1} = 0.1;
   
    birth_means{2} = [250,400,0,0]';
    birth_covariances{2} = diag([100,100,50,50]);
    birth_weights{2} = 0.1;

    J_gamma = length(birth_means);

    seperation_weights = {};
    seperation_covariances = {};

    seperation_weights{1} = 0.1;
    seperation_covariances{1} =diag([100,100,50,50]);

    J_beta = length(seperation_weights);
    %A_beta = zeros(size(A));
    A_beta = A;
    d_beta = 20;


    means = birth_means;
    covariances = birth_covariances;
    weights = birth_weights;

    % number of targets
    J = 0;
    
    % truncation threshold
    T = 1e-5;
    
    % merging threshold
    U = 1;

    % maximum number of Gaussian terms
    J_max = 100;

        

    for k = time_steps
        i = 0;
        predicted_means = {};
        predicted_covariances = {};
        predicted_weights = {};


        for j = 1:J_gamma
            i = i + 1;
            predicted_weights{i} = birth_weights{j};
            predicted_means{i} = birth_means{j};
            predicted_covariances{i} = birth_covariances{j};
        end
        
        for j = 1:J_beta
            for l = 1:J
                i = i + 1;
                predicted_weights{i} = weights{l}*seperation_weights{j};
                predicted_means{i} = d_beta + A_beta*means{l};
                predicted_covariances{i} = seperation_covariances{j} + A_beta * covariances{l} * A_beta';
            end
        end

        for j = 1:J
            i = i + 1;
            predicted_weights{i} = P_s * weights{j};
            predicted_means{i} = A * means{j};
            predicted_covariances{i} = G * Q * G' + A * covariances{j} * A';
        end
        
        predicted_J = i;

        predicted_etas = cell(1,predicted_J);
        predicted_innovation_covariances = cell(1,predicted_J);
        Kalman_gains = cell(1,predicted_J);
        estimated_covariances = cell(1,predicted_J);

        
        for j = 1:predicted_J
            predicted_etas{j} = C * predicted_means{j};
            predicted_innovation_covariances{j} = R + C * predicted_covariances{j} * C';
            Kalman_gains{j} = predicted_covariances{j} * C' / predicted_innovation_covariances{j};
            estimated_covariances{j} = (eye(length(predicted_covariances{j})) - Kalman_gains{j} * C) * predicted_covariances{j};
        end

        for j = 1:predicted_J
            weights{j} = (1 - P_d) * predicted_weights{j};
            means{j} = predicted_means{j};
            covariances{j} = predicted_covariances{j};
        end

        l = 0;
        xvalues = detections{2,k+1};
        yvalues = detections{3,k+1};
        for d = 1:length(xvalues)
            if d ~=0
                z = [xvalues(d);yvalues(d)];
                l = l + 1;
                for j = 1:predicted_J 
                    weights{l*predicted_J+j} = P_d * predicted_weights{j} * mvnpdf(z,predicted_etas{j},predicted_innovation_covariances{j});
                    means{l*predicted_J+j} = predicted_means{j} + Kalman_gains{j} * (z-predicted_etas{j});
                    covariances{l*predicted_J+j} = estimated_covariances{j};
                end
                weightsum = 0;
                for i = 1:predicted_J
                    weightsum = weightsum + weights{l*predicted_J+i};
                end
    
                for j = 1:predicted_J
                    kappa = beta_FA * V * 1/V;
                    weights{l*predicted_J+j} = weights{l*predicted_J+j}/(kappa+weightsum);
                end
            end
        end
        J = l*predicted_J+predicted_J;
        [weights,means,covariances,J] = pruning(weights,means,covariances,T,U,J_max);

        
        estimated_states{k+1} = extract_multiple_target_state(weights,means);
    end

end