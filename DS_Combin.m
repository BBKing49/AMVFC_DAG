function fin_u = DS_Combin(U)
    % alpha: Cell array containing all Dirichlet distribution parameters.
    % Return: Combined Dirichlet distribution parameters.
    

    classes = size(U{1}, 2); % Assuming all alpha have the same number of classes
    for v = 1:size(U,1)-1
        if v == 1
            alpha_a = DS_Combin_two(U{1}, U{2}, classes);
        else
            alpha_a = DS_Combin_two(alpha_a, U{v+1}, classes);
        end
    end
    fin_u = alpha_a-1;
end

function new_alpha = DS_Combin_two(alpha1, alpha2, classes)
        % alpha1: Dirichlet distribution parameters of view 1
        % alpha2: Dirichlet distribution parameters of view 2
        % classes: Number of classes
        % Return: Combined Dirichlet distribution parameters
        
        N = size(alpha1,1);
        alpha = {alpha1, alpha2};
        b = cell(1, 2);
        S = cell(1, 2);
        E = cell(1, 2);
        u = cell(1, 2);
        
        for v = 1:2
            S{v} = sum(alpha{v}, 2);
            E{v} = alpha{v} - 1;
            b{v} = E{v} ./ S{v};
            u{v} = classes ./ S{v};
        end
        
        % b^0 @ b^(0+1)
        for i = 1:N
            sumC = 0;
            for j  = 1:classes
                bb(i,j) = b{1}(i, j)*b{2}(i, j);
                b1u2(i,j) = b{1}(i, j).*u{2}(i);
                b2u1(i,j) = b{2}(i, j).*u{1}(i);
                for k =1:classes
                    if j~=k
                        sumC = sumC + b{1}(i, j)*b{2}(i, k);
                    end
                end
            end
            ba = (bb(i,:)+b1u2(i,:)+b2u1(i,:))/(1-sumC);
            ua = (u{1}(i,:).*u{2}(i,:))/(1-sumC);        
            new_S = classes/ua;
            new_e = new_S*ba;
%             new_e = mapminmax(new_e,0,1)';
            new_alpha(i,:) = new_e+1;
        end

end