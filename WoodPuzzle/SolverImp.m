%Solves a wooden puzzle of interlocking pieces, the outer perimiter and the square ring 2 blocks
%inward doesn't need ot be covered, but everything else must be.

%This takes the semi-naive approach of trying all permutations. It is only semi-naive because:
%(1)
%it discards intermediate states (where all perms of a new piece has been added) with overlapping
%pieces. Initially few states have overlapping pieces, so the total number of intermediate states
%continues to grow with each new piece added. As pieces are added the number of new states slows as
%there are far fewer legal places for the new piece. It then decays to far fewer, and then only
%solved states at the end.
%(2)
%It also looks for "hopeless" intermediate states. These are ones where a square that needs to be
%covered can never be covered becuase it is already surrounded. Currently it only looks for easy
%hopeless states where the square is surrounded directly adjacent on the NSEW sides. You could get
%fancier by looking for less obvious hopeless states where  the blockage isn't immediately
%adjacent, but past a certain blocked area size this becomes trickier since small pieces could fit.

%TODO:
%-Read in arbitrary piece shapes from file
%-Automatically select a piece addition order to stall state number growth as best as we can until
%the state growth drops off from increased piece density
%-Detect piece symmetry and tag them to be treated specially during perms generation

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

sw=[3,4,6,5,2,1,7,8,9]; %Magic! Manual ordering choosing the smallest set of perms for each piece first
for n=1:numel(shapes)
tt{n}= shapes{sw(n)};
end
shapes= tt;

L=10; %Side length of the puzzle

for n= 1:numel(shapes)
    s= shapes{n};
    pos= prod(L-size(s)+1); %Number of positions the piece can be translated to within the board
    temp= zeros(L,L,pos*4);
    it= 1;
    I= size(s,1)-1;
    J= size(s,2)-1;
    for j=1:L-size(s,2)+1
        for i=1:L-size(s,1)+1
            %Create an array of all translations and rotations
            temp(i:i+I,j:j+J,it+pos*0)=rot90(s,0);
            temp(j:j+J,i:i+I,it+pos*1)=rot90(s,1);
            temp(i:i+I,j:j+J,it+pos*2)=rot90(s,2);
            temp(j:j+J,i:i+I,it+pos*3)=rot90(s,3);
            it=it+1;
        end
    end
    if or(n==1,n==6); temp=temp(:,:,1:pos); end %The first piece does not need to be rotated, since the second piece that it is combined with already has been. Piece 6 is symmetrical
    perms{n}= temp;
    %if n~=6; perms{n}(:,:,end+1:end+size(temp,3))= reshape(flipud(reshape(temp,L,[])),L,L,[]); end %Also account for the piece being flipped over, 6 is symmetric
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
pat= perms{1}(:,:,r)+2*perms{2}(:,:,c);
%used{1}=r; used{2}=c;

mem_limit= 3e9; %This is the max size of "over" so that find() doesn't fail (on 16Gb systems), if it is we should split the generating perms{n} into pieces
for n= 3:numel(perms)
    splitB= round(linspace(1,size(over,4),1+ceil((size(perms{n},3)*size(over,4)*L^2)/mem_limit)));
        splitB= splitB(2:end);
    splitA= [1, splitB(1:end-1)+1];
    fprintf('Splitting into %i pieces \n', numel(splitA))
    
    seg=0; seg2=0; ind=[]; ind2=[]; overt=[];
    for p= 1:numel(splitA)
        overAB= over(:,:,1,splitA(p):splitB(p)); %Split "over" into memory fitting chunks that span AB
        overAB= reshape(bsxfun(@plus, perms{n}, overAB),L,L,1,[]); %Can't overwrite previous full intermediate state "over", since we will need it for later perm pieces

        indt= find(squeeze(not(any(any(overAB>1)))));
        ind= [ind; seg+indt];
        seg= seg+ size(overAB,4);
        fprintf('%i %i %9i ',n,p,seg)
        fprintf('%8i ',length(ind))
        overAB= single(overAB(:,:,1,indt));
        
        if n== numel(perms) %We don't worry about hopeless states after the last piece is placed
            ind2t= [1:size(overAB,4)]';
        else
            ind2t= find(squeeze(not(any(any((mtimesx(Ay,overAB)+mtimesx(overAB,Ax))>3)))));
        end
        ind2= [ind2; seg2+ind2t];
        seg2= seg2+ size(overAB,4);
        fprintf('%7i %f \n',length(ind2),seg/length(ind2))
        overt= cat(4, overt, int8(overAB(:,:,1,ind2t))); %Concat with previous piece solutions
    end
    c= floor((ind(ind2)-1)/size(perms{n},3))+1;
    r= ind(ind2)-(c-1)*size(perms{n},3);
    pat= pat(:,:,c)+n*perms{n}(:,:,r);
    %usedP{n-2}=c; used{n}=r;
    over= overt;
end
test= [1,1,1,1,1,1,1,1,1,1;1,0,0,0,0,0,0,0,0,1;1,0,1,1,1,1,1,1,0,1;1,0,1,0,0,0,0,1,0,1;1,0,1,0,0,0,0,1,0,1;1,0,1,0,0,0,0,1,0,1;1,0,1,0,0,0,0,1,0,1;1,0,1,1,1,1,1,1,0,1;1,0,0,0,0,0,0,0,0,1;1,1,1,1,1,1,1,1,1,1;];
okpat=pat(:,:,all(all(bsxfun(@or,test,pat))));
toc