function genSynth(N, T, sig)
% N: number of time series (a multiple of 4)
% T: Length of the time series
% sig: variance of the noise process

K = [0.9 0 0 0; 1 0.9 0 0; 1 0 0.9 0; 1 0 0 0.9];
A = kron(eye(N/4), K);
A = A./norm(A).*0.5;
series = zeros(N, T);
series(:, 1) = randn(N, 1);
for t = 2:T
    series(:, t) = A*series(:, t-1) + sig*randn(N, 1);
end

save('synth.mat', 'series', 'A')
