L=100;
S= single(round(rand(L)));
A= single(toeplitz([[1 1] zeros(1,L-2)]));

for it=1:3000
    imagesc(S); drawnow
    S=single(min(1,max(0,int8(S)+1-mod(int8(A*S*A)+4,7))));
end