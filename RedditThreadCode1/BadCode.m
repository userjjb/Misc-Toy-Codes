close all; clear;
 
N = 4;  % number of elements
 
theta_scan = 1:181;
theta_deg = 15;
theta_rad = theta_deg*pi/180;
 
phase_deg = [0 180];
phase_rad = phase_deg*pi/180;
 
% Compute Array Factor
 
AF1111 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF1112 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF1121 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF1122 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF1211 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF1212 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF1221 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF1222 = exp(1i*phase_rad(1,1)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
 
AF2111 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF2112 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF2121 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF2122 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,1))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF2211 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF2212 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,1))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
AF2221 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,1)));
AF2222 = exp(1i*phase_rad(1,2)) + exp(1i*(pi*sin(theta_rad) + phase_rad(1,2))) + ...
    exp(1i*(2*pi*sin(theta_rad) + phase_rad(1,2))) + exp(1i*(3*pi*sin(theta_rad) + phase_rad(1,2)));
 
% Combine Results
AF = [AF1111 AF1112 AF1121 AF1122 AF1211 AF1212 AF1221 AF1222 ...
    AF2111 AF2112 AF2121 AF2122 AF2211 AF2212 AF2221 AF2222];