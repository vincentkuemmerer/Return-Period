function [theta,M,Z] = ExtremalIndex(Hs,Threshold)

W = Hs>Threshold; % Binary vector of exceedances
N = nnz(W); % Number of exceedances
S = find(W); % Index of ith exceedance
T = zeros(1,length(N)); % Preallocation of T (inter-exceedance time)

for j = 1:(N-1) % can be rewritten with diff(T), see https://arxiv.org/pdf/1605.07006.pdf
     T(j) = S(j+1) - S(j); % inter-exceedance time
end

theta = (2*(sum(T-1))^2) / ((N-1)*sum((T-1).*(T-2))); % when max(T)>2, theta is extremal index estimator
    
% Using extremal index estimator 
M = 1 / theta; % Mean cluster size M
Z = theta * N; % Total cluster count
end
