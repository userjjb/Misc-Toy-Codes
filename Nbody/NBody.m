clc
clear all
close all

dt = 0.00001;
N = 600;
vmag = 5;
rng(0, 'simdTwister')
theta = rand(N, 1)*2*pi;
%rr = rand(N, 1)*0.4 + 0.1;
t = 1./[(N-1)/4:5*(N-1)/4]';
rr = cumsum((0.5/sum(t))*t);
x1 = rr.*cos(theta);
x2 = rr.*sin(theta);
q = rand(N, 1) + 1;
% v1 = vmag*(x2);
% v2 = vmag*(-x1);
vr = x1.^2 + x2.^2;
v1 = vmag*-x2./vr;
v2 = vmag*x1./vr;

s=[0:1000]/1000;
p=polyfitB(s,s.^(3/2),2,1e-4);

for i = 0:10000
    r1 = bsxfun(@minus,x1,x1');
    r2 = bsxfun(@minus,x2,x2');

    rmag = r1.^2 + r2.^2;
    %rmag = (rmag.*polyval(p,rmag))+1e-3;
    rmag = rmag.^2*p(1) + rmag*p(2) +p(3);

    r1 = r1./rmag;
    r2 = r2./rmag;

    r1(1:N+1:end) = 0;
    r2(1:N+1:end) = 0;

    v1 = v1 - r1*q*dt;
    v2 = v2 - r2*q*dt;

    x1 = x1 + v1*dt;
    x2 = x2 + v2*dt;
    
    if not(mod(i,100))
        scatter(x1,x2,'.r')
        %hold on
        %quiver(x1,x2,v1,v2)
        axis([-1.5 1.5 -1.5 1.5])
        axis square
        drawnow
    end
end