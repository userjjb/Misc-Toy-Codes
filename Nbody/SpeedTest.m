clc
clear all
close all

s=[0:1000]/1000;
p=polyfitB(s,s.^(3/2),2,1e-4);

dt = 0.00001;
N = 2^7* 100;
vmag = 80;
rng(0, 'simdTwister')
theta = rand(N, 1)*2*pi;
rr = rand(N, 1)*0.4 + 0.1;
x1 = rr.*cos(theta)+1.5;
x2 = rr.*sin(theta)+1.5;
q = rand(N, 1) + 1;
v1 = vmag*(x2-1.5);
v2 = vmag*(1.5-x1);
tic
r1 = bsxfun(@minus,x1,x1');
r2 = bsxfun(@minus,x2,x2');

rmag = r1.^2 + r2.^2;
rmag = (rmag.*sqrt(rmag))+1e-4;
%approximate version (faster) \/
%rmag = rmag.^2*p(1) + rmag*p(2) +p(3);

r1 = r1./rmag;
r2 = r2./rmag;

r1(1:N+1:end) = 0;
r2(1:N+1:end) = 0;

v1 = v1 - r1*q*dt;
v2 = v2 - r2*q*dt;

x1 = x1 + v1*dt;
x2 = x2 + v2*dt;
toc