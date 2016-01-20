close all
L=500;
%S= single(round(rand(L)));
S=zeros(L);
A= single(toeplitz([[1 1] zeros(1,L-2)]));
colormap([1 1 1; 0 0 0]);

acorn=[0,1,0,0,0,0,0;0,0,0,1,0,0,0;1,1,0,0,1,1,1;];
tomas=[0,0,1,0,0;1,1,0,0,0;0,1,0,0,0;1,0,0,1,0;0,0,0,0,1;0,1,0,0,1;0,0,1,0,1;0,1,0,0,0;];

thing= tomas;
a= round(size(thing)/2);
b= size(thing)-1-a;
S(L/2-a(1):L/2+b(1),L/2-a(2):L/2+b(2))=thing;

for it=1:40000
    S=single(sign(max(0,int8(S)+1-mod(int8(A*S*A)+4,7))));
    if mod(it,30)==0
        imagesc(S); axis equal; axis([0 L 0 L]); text(0,0,num2str(it)); drawnow;
    end
end