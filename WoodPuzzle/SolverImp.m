clear all; close all
tic
shapes=...
{
    [0 1 0;
     1 1 1;
     0 1 0],
     
    [0 0 1 0 0;
     0 0 1 0 0;
     1 1 1 1 1;
     1 0 1 0 0],
     
    [0 1 0 0 0;
     0 1 0 0 1;
     1 1 1 1 1;
     0 1 0 0 0;
     0 1 0 0 0],
     
    [0 0 1 1 0;
     0 0 0 1 0;
     0 0 0 1 0;
     1 1 1 1 1;
     0 0 0 1 0],
     
    [0 0 1 0 0;
     0 0 1 0 0;
     1 1 1 1 1;
     0 0 1 0 0;
     0 1 1 0 0],
     
    [0 0 1 0;
     1 1 1 1;
     0 0 1 0;
     0 0 1 0;
     0 1 1 0],
     
    [0 1 0;
     0 1 0;
     1 1 1;
     0 1 0;
     1 1 0],
     
    [0 1 0;
     1 1 1;
     0 1 0;
     1 1 0],
     
    [0 0 1 1;
     0 0 1 0;
     1 1 1 1;
     0 0 1 0]
};

sw=[3,4,6,5,2,1,7,8,9];
for n=1:numel(shapes)
tt{n}= shapes{sw(n)};
end
shapes= tt;

L=10;

for n= 1:numel(shapes)
    s= shapes{n};
    pos= prod(L-size(s)+1);
    temp= zeros(L,L,pos*4);
    it= 1;
    I= size(s,1)-1;
    J= size(s,2)-1;
    for j=1:L-size(s,2)+1
        for i=1:L-size(s,1)+1
            temp(i:i+I,j:j+J,it+pos*0)=rot90(s,0);
            temp(j:j+J,i:i+I,it+pos*1)=rot90(s,1);
            temp(i:i+I,j:j+J,it+pos*2)=rot90(s,2);
            temp(j:j+J,i:i+I,it+pos*3)=rot90(s,3);
            it=it+1;
        end
    end
    if or(n==1,n==6); temp=temp(:,:,1:pos); end
    perms{n}= temp;
    if n~=6; perms{n}(:,:,end+1:end+size(temp,3))= reshape(flipud(reshape(temp,L,[])),L,L,[]); end
end

for n=1:numel(perms); perms{n}= int8(perms{n}); end

A= toeplitz([[-1 1] zeros(1,L-2)]);
Ax= A; Ax(:,[1,3,end-2,end])=0; Ax(1:3,2)=[0,-1,2]; Ax(end-2:end,end-1)=[2,-1,0]; Ay=Ax';

over= reshape(bsxfun(@plus,perms{1},reshape(perms{2},L,L,1,[])),L,L,1,[]);
fprintf('2 %9i ',size(over,4))

ind= find(squeeze(not(any(any(over>1)))));
fprintf('%8i %f',length(ind),size(over,4)/length(ind))
over= single(over(:,:,1,ind));

ind2= find(squeeze(not(any(any((mtimesx(Ay,over)+mtimesx(over,Ax))>3)))));
fprintf('%7i \n',length(ind2))
over= int8(over(:,:,1,ind2));

c= floor((ind(ind2)-1)/size(perms{1},3))+1;
r= ind(ind2)-(c-1)*size(perms{1},3);
used{1}=r; used{2}=c;
pat= perms{1}(:,:,r)+2*perms{2}(:,:,c);

for n= 3:numel(perms)
    over= reshape(bsxfun(@plus, perms{n}, over),L,L,1,[]);
    fprintf('%i %9i ',n,size(over,4))
    seg= size(over,4)/4;
    
    ind= find(squeeze(not(any(any(over>1)))));
    fprintf('%8i ',length(ind))
    over= single(over(:,:,1,ind));
    
    ind2= find(squeeze(not(any(any((mtimesx(Ay,over)+mtimesx(over,Ax))>3)))));
    fprintf('%7i %f \n',length(ind2),(seg*4)/length(ind2))
    over= int8(over(:,:,1,ind2));
    
    c= floor((ind(ind2)-1)/size(perms{n},3))+1;
    r= ind(ind2)-(c-1)*size(perms{n},3);
    usedP{n-2}=c; used{n}=r;
    pat= pat(:,:,c)+n*perms{n}(:,:,r);
end
toc 