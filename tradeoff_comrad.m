function [X] = tradeoff_comrad(rho,H,Y,power,X_dir)
[~,N] = size(H);
[~,L] = size(Y);
p1 = sqrt(rho);
p2 = sqrt(1-rho);
A = [p1*H.',p2*eye(N)].';
B = [p1*sqrt(power)*Y;p2*X_dir];
Q = A'*A;
G = A'*B;
%G_inv = inv(G*G');
e = eig(Q);
[V,~] = eig(Q);
l_min = min(e);
for i = 1:N
    tr_fe(i) = real(trace (G*G'*V(:,i)*V(:,i)'));
end
func = @(xx)sum(tr_fe'./(e+xx*ones(N,1)))+xx*power*L;
LB = -l_min;
UB = 1;
EPSILON = 10^-4;
l = LineSearchGoldenSection(func,LB,UB,EPSILON);
X = pinv(Q+l*eye(N))*G;
end

