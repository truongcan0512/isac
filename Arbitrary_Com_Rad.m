function [ X ] = Arbitrary_Com_Rad(H,Y,power,F)
[K,N] = size(H);
[~,L] = size(Y);
H_wave = H*F;
[U,S,V] = svd(sqrt(L*power)*H_wave'*Y);
X = sqrt(L)*F*U*eye(N,L)*V';
end