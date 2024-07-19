function [sigOut] = reconstructSig2(omega,sigIn, K, eT12)
%reconstructSig2(omega,sigIn, k)

%RECONSTRUCTSIG reconstruct the signal based on GC source and ommega
%   Inverse of Granger causality (GC)
%
% Author: Peter Rogelj, peter.rogelj@upr.si
%
%  sigOut_(n) = (1-K)*T1 +
%                  K *(  sum_(i=1..N) ( omegaT_(i) * T1_(n-i) +
%                        sum_(i=1..N) ( omegaS_(i) * S1_(n-i) +
%                        eT12 ) 
%
% Inputs:
%  omega - weights for past samples, order N -> N samples, such that the
%          first weight is for the latest sample, last for the oldest,
%          size=[2,N] = [omegaT; omegaS]
%  sigIn - signal from which to reconstruct [T,S]
%  K     - temporal weights (the same size as sigIn)
%  eT12  - innovative process - regression error after bivariate regression
%
%  
% Output:
%  sigOUT - the result of reconstruction
if nargin<3
    k=ones(size(sigIn));
else
    if size(sigIn,2)~=size(K,2) 
        error('Invalid size of K - should be equal to the size of sigIn!')
    end
end
if nargin<4
    eT12=zeros(1,size(sigIn,2));
end

sigOut = zeros(1,size(sigIn,2));
N = size(omega,2);

for snr=1:2
    firstColumn = [zeros(1,N-1), sigIn(snr,1)];
    Xin = flipud(hankel(firstColumn, sigIn(snr,:)));
    sigOut= sigOut + omega(snr,:) * Xin;
end
sigOut = sigOut + eT12;
sigOut = K.*sigOut;




sigOut = sigOut + (1-K).*sigIn(1,:); 
end

