function [q,B,J] = compute_qandB(h,E,T)
    q = zeros(length(E),1);
    for i=1:length(E)
        if h(E(i,1)) >= h(E(i,2))
            B(i,E(i,1)) = -1; B(i,E(i,2)) = 1;
            q(i) = ((1/T(i,i))*(h(E(i,1))-h(E(i,2))))^(1/1.852);
        else
            B(i,E(i,2)) = -1; B(i,E(i,1)) = 1;
            q(i) = ((1/T(i,i))*(h(E(i,2))-h(E(i,1))))^(1/1.852);
        end
    end
    J = -B'; % matrix B is necessary to retrieve the demands (not used now)
                 % but useful to keep it
end