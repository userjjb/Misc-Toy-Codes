clear all
L=6;

scramble= randi(L^2,1,12); %Sequence of random light positions to push
s=  [1 0 0 0 1;...
     0 0 0 0 0;...
     0 0 0 0 0;...
     0 0 0 0 0;...
     1 0 0 0 1];
s= zeros(L);
s= uint8(s(:)); %uint8 for speed when later using xor()

moves= zeros(L,L,L,L); %Set of responses to all possible light pushes 
for i=1:L
    for j=1:L
        a=i-1;
        b=j-1;
        k=i+1;
        l=j+1;
        if i==1
            a=L;
        elseif i==L
            k=1;
        end
        if j==1
            b=L;
        elseif j==L
            l=1;
        end
        moves(i,j,i,j)=1;
        moves(k,j,i,j)=1;
        moves(i,l,i,j)=1;
        moves(a,j,i,j)=1;
        moves(i,b,i,j)=1;
    end
end
moves= uint8(reshape(moves,L^2,[])); %This matrix defines the full system response and is actually block diagonal

for i=1:numel(scramble)
    s= xor(s,moves(:,scramble(i))); %xor() is how we switch lights on and off with a move since xor(0,1)->1 and xor(1,1)->0
end

% seen= uint32(hash*double(s)); %Un-comment if using seen/unseen filtering below
hash= 2.^[0:L^2-1]; %hash just assigns the unique decimal value to the binary array
p=s; in=0;
for i=1:20 %Set some maximum depth to search
    [r,c]= find(p);
    rsvd{i}= r;     %r are the numeric positions of each of the possible lights to flip for a given permutation
    csvd{i}= in(c); %c is the numbering of the permutations to try, in order; we will try flipping each of the lights in each of these perms. We save c(in), which is the ordering *before* we eliminated non-unique perms
    on= (c-1)*L^2+r; %'on' is the relative column position of each the lights we *want* to flip, this excludes the columns that we flip needlessly in the next step
    p= reshape(bsxfun(@xor,p,moves),L^2,[]); %for speed we don't flip only the lights that are on for a given perm column, otherwise we couldn't vectorize across all perms with the same operation, instead we try flipping all lights
    p= p(:,on);     %now we pick out only the results for the cases where there *actually* was a light on
    [h,in]= unique(uint32(hash*p)); %to identify duplicates we calc a hash for all perms and find only the unique perms
    fprintf('%i %i \n',numel(on),numel(h))
    if find(h==0) %Adjust solution hash for non-empty boards if desired
        trace=in(find(h==0));
        break
    end
    %If we search a significant fraction of the solution space, we will start to see many repeated
    %states that have been previously iterated from. To save the trouble of re-computing work
    %already done we should detect these previously seen states. NOTE: This is fairly expensive due
    %to the sizes of 'h' and the full list of previously seen states, so this isn't usually worth it
%     unseen= find(sum(bsxfun(@eq,h,seen'),1)<1);
%     seen= [seen h(unseen)];
%     fprintf('%i \n',numel(h)-numel(unseen))
    p= reshape(p(:,in),L^2,1,[]); %If using the above seen/unseen technique replace 'in' with 'in(unseen)'
    if i==20
        fprintf('Exiting without solution found')
    end
end

for i=numel(rsvd):-1:1 %Reconstruct solution permutations
    solve(i)= rsvd{i}(trace);
    trace= (csvd{i}(trace));
end
