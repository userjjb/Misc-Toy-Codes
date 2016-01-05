%close all; clear;
 
N = 4;  % number of elements

theta_deg = 15;
T= pi*sin(theta_deg*pi/180);
 
phase_deg = [0 180];
Nphase= length(phase_deg);
phase_rad = phase_deg*pi/180;
 
% Compute Array Factor
Ind= ceil(bsxfun(@rdivide,bsxfun(@mod,[0:Nphase^N-1]',Nphase.^[N:-1:1])+1,Nphase.^[N-1:-1:0]));
AFnew = permute(sum(exp(1i*bsxfun(@plus,phase_rad(Ind),T*[0:N-1])),2),[2 1])
 
%AF1111 = exp(1i*phase_rad(1)) + exp(1i*(T + phase_rad(1))) + ...
    %exp(1i*(2*T + phase_rad(1))) + exp(1i*(3*T + phase_rad(1)));