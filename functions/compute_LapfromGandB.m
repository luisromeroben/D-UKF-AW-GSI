function [Quad,Omega,Phi] = compute_LapfromGandB(h,G,B,smooth,WA,E)

N=length(h);

Omega = zeros(N);
for p = 1:length(E)
    i = E(p,1); j = E(p,2);
    Omega(i,j) = G(i,j)^0.54*(B(i,j)*(h(i)-h(j)))^(-0.46);
    Omega(j,i) = Omega(i,j);
end

dPhi = zeros(N,1);
for i=1:N
    neig = find(WA(i,:)>0);
    for k = neig
    dPhi(i) = dPhi(i) + Omega(i,k);
    end
end

dPhi=sum(Omega);
Phi=diag(dPhi);

% Phi = diag(dPhi); %Phi = diag(sum(Omega));

if smooth == 1
    Quad = Phi - Omega;
else
    Phim2 = diag(dPhi.^(-2));
    Quad = (Phi - Omega)'*Phim2*(Phi - Omega);
end

end