function [ X ] = Orthogonal_Com_Rad(H,Y,power)
[K,N] = size(H);
[~,L] = size(Y);
u = sqrt(power/N);
[U,S,V] = svd(sqrt(L/N)*power*H'*Y);
X = u*sqrt(L)*U*eye(N,L)*V';
end

