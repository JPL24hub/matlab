function [Cwf, Cep] = MIMOCapacity(R_t, R_r, SNR)
% This function calculates the ergodic and outage capacity of a MIMO Rayleigh 
% channel considering no CSIT (equal power allocation) and perfect CSIT
% (waterfilling power allocation). In both cases perfect CSIR is assumed. The
% channel is assumed to be spatially correlated according to a Kronecker
% model but temporally uncorrelated.
% -------------------------------------------------------------------------
% Inputs:
% R_t : Transmit correlation matrix of size n_t*n_t
% R_t : Receive  correlation matrix of size n_r*n_r
% SNR : Average SNR at each receive antenna in dB
% Outputs:
% Cwf : Waterfilling capacity structure (containing ergodic and outage capacity)
% Cep : Equal power  capacity structure (containing ergodic and outage capacity)
% Example:
% calculating the capacity of a 3*4 i.i.d. Rayliegh channel at 20dB SNR
% [Cwf Cep] = MIMOCapacity(eye(e), eye(4), 20)
% -------------------------------------------------------------------------
% Omid Darvishi
% 2008-03-24
SNR = 10 ^ (SNR/10);
n_t = length(R_t);      % number of Tx antennas
n_r = length(R_r);      % number of Rx antennas
Pout = 0.1;             % outage probability
nsample = 1000;         % number of channel realizations
% nsample = 100 / Pout;
R_r_sqrt = R_r^0.5;
R_t_sqrt = R_t^0.5;
cwf = zeros(nsample, 1);
cep = zeros(nsample, 1);
Hw = (randn(n_r, n_t, nsample) + sqrt(-1) * randn(n_r, n_t, nsample)) / sqrt(2);
for i = 1:nsample
    H = R_r_sqrt * Hw(:,:,i) * R_t_sqrt;
    lambda = eig(H * H');
    [mu pl] = WFL(lambda, SNR);
    cwf(i) = sum(log2(1 + pl .* lambda));
    cep(i) = sum(log2(1 + SNR/n_t * lambda));
end
Ce = mean([cwf cep]);
C  = sort([cwf cep]); 
Co  = C(round(Pout * nsample),:);
Cwf.ergodic = Ce(1);
Cep.ergodic = Ce(2);
Cwf.outage  = Co(1);
Cep.outage  = Co(2);
function [mu, pl] = WFL(lambda, SNR)
% Calculates water filling level.
% Inputs
%     lambada: vector of the channel(H*H') eigenvalues 
%     SNR:   SNR(linear)
% Outputs
%     mu: water filling level
%     pl: vector of the power level of eigen-modes
    
[lambda idx] = sort(lambda, 'descend');
lambda = lambda(find(lambda > 0));      % ignoring non-positive eigenvalues
pl = -1;
try
    while (min(pl) < 0)
        mu = (SNR + sum(1 ./ lambda)) / length(lambda);
        pl = mu - 1 ./ lambda;
        lambda = lambda(1:end-1);
    end
catch
    disp('There exists no water filling level for the input eigenvalues. Check your data and try again')
end
pl = [pl; zeros(length(idx) - length(pl), 1)]; % assigning zero power for weak eigen-modes 
pl(idx) = pl; % rearranging the power levels 