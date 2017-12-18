close all
clear all

n = 4;
test1 = [4,3,4,2;1,1,4,2;5,3,2,1;5,1,5,2;];
test2 = [4,3,4,2;1,1,4,4;5,3,4,1;5,1,5,2;];
A = toeplitz([[2 1] repmat(0,1,n-2) [1] repmat(0,1,n^2-(n+1))]);



