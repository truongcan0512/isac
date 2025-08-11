function [ R_d ] = my_waveform_design(T,A_t,A_i,power)

% Radar waveform design with CVX
eta = 1.8;                            % PAPR upper bound (e.g., max 2x average)
cvx_begin sdp
cvx_begin quiet
variable R_d(T, T) hermitian semidefinite
maximize(sum(diag(A_t'*R_d*A_t)) - sum(diag(A_i'*R_d*A_i)) + eps)  % total gain at targets
subject to
trace(R_d) == power;                             % total power constraint
for n = 1:T
    R_d(n,n) <= eta * (1/T) * trace(R_d);      % PAPR constraint per antenna
end
cvx_end
end