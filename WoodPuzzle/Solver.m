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
    if or(n<2,n==6); temp=temp(:,:,1:pos); end
    perms{n}= temp;
end
for n=1:numel(perms); perms{n}= int8(perms{n}); end

over= reshape(bsxfun(@plus,perms{1},reshape(perms{2},L,L,1,[])),L,L,1,[]);
ind= find(squeeze(not(any(any(over>1)))));
over= over(:,:,1,ind);
c= floor((ind-1)/size(perms{1},3))+1;
r= ind-(c-1)*size(perms{1},3);
used{1}=r; used{2}=c;
pat= perms{1}(:,:,r)+2*perms{2}(:,:,c);

for n= 3:numel(perms)
    over= reshape(bsxfun(@plus, perms{n}, over),L,L,1,[]);
    seg= size(over,4)/4;
    size(over,4)
    ind= find(squeeze(not(any(any(over(:,:,1,1:end/4)>1)))));
    ind= [ind; seg+find(squeeze(not(any(any(over(:,:,1,end/4+1:end/2)>1)))))];
    ind= [ind; 2*seg+find(squeeze(not(any(any(over(:,:,1,end/2+1:3*end/4)>1)))))];
    ind= [ind; 3*seg+find(squeeze(not(any(any(over(:,:,1,3*end/4+1:end)>1)))))];
    length(ind)
    over= over(:,:,1,ind);
    c= floor((ind-1)/size(perms{n},3))+1;
    r= ind-(c-1)*size(perms{n},3);
    usedP{n-2}=c; used{n}=r;
    pat= pat(:,:,c)+n*perms{n}(:,:,r);
end
toc