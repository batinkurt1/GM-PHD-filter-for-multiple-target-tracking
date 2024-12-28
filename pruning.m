function [w_pruned,m_pruned,P_pruned,l] = pruning(w,m,P,T,U,J_max)
    w_pruned = {};
    m_pruned = {};
    P_pruned = {};
    
    l = 0;
    w = cell2mat(w);
    I = find(w > T);

    while ~isempty(I)
        l = l + 1;
        [~, j_idx] = max(w(I));
        j = I(j_idx);

        L = [];
        for i = I
            a = m{i};
            b = m{j};
            mahalanobis_distance = (a - b)' / P{j} * (a - b);
            if mahalanobis_distance <= U
                L = [L, i];
            end
        end
        
        w_tilde = 0;
        m_tilde = zeros(size(m{1}));
        P_tilde = zeros(size(P{1}));

        for i = L
            w_tilde = w_tilde + w(i);
        end

        for i = L
            m_tilde = m_tilde + (1/w_tilde)*(w(i)*m{i});
        end

        for i = L
            P_tilde = P_tilde + (1/w_tilde)*(w(i)*(P{i}+(m_tilde - m{i})*(m_tilde - m{i})'));
        end

        w_pruned = [w_pruned, w_tilde];
        m_pruned = [m_pruned, m_tilde];
        P_pruned = [P_pruned, P_tilde];
        I = setdiff(I, L);

        if l > J_max
            [w_pruned, sort_idx] = sort(cell2mat(w_pruned), 'descend');
            
 
            m_pruned_sorted = {};
            P_pruned_sorted = {};

            for i = sort_idx
                m_pruned_sorted{end + 1} = m_pruned{i};
                P_pruned_sorted{end + 1} = P_pruned{i};
            end
            w_pruned = num2cell(w_pruned(1:J_max));
            m_pruned = m_pruned(1:J_max);
            P_pruned = P_pruned(1:J_max);
            l = J_max;
            break;
        end
    end

end