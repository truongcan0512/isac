function [ R ] = waveform_design_multibm_covmat_new( Pd_theta,T,a,theta,power)
% The beampattern matching design (see page 193 in textbook)
cvx_solver sedumi
cvx_begin quiet
variable alph % scale factor
variable R(T,T) hermitian semidefinite
expression u(length(theta))
for i=1:length(theta)
    u(i)=(alph*Pd_theta(i)-a(:,i)'*R*a(:,i));
end
minimize norm(u,2)
subject to
% diag(R)==ones(N,1)*power/N;
trace(R)==power;
cvx_end
end

