function [r] = solver(a,b,n)

N = numel(a);
maxtries = ceil(min(N/1.2,450));
X = size(a,1);
Y = size(a,2);
popmap = a+b;
apop = sum(a(:));
bpop = sum(b(:));
ideal_a = apop/n;
ideal_b = bpop/n;
ideal_pop = ideal_a + ideal_b;

RRsave{1} = basic_alloc(a,b,n);
RRsave{2} = basic_alloc(a',b',n)';

basesize = floor(N/n);

sizes = ones(n,1)*basesize;
temp = mod(N,n);
sizes(end-temp+1:end) = sizes(end-temp+1:end)+1;

minscore  = 1e5;
dpop=zeros(n,1); anet = dpop; poperr = dpop;
for k = 1:2
    rr = RRsave{k};
    
    for i=1:n
        rri=rr==i;
        dpop(i) = sum(popmap(rri));
        anet(i) = sum(a(rri))-sum(b(rri));
        poperr(i) = dpop(i) - ideal_pop;
    end
    totscore = sum(abs(poperr)) * max(abs(anet));
    if totscore < minscore,
        minscore = totscore;
        RRbest = RRsave{k};
        used = k;
        %disp([k minscore])
    end
end
      
RR=RRbest;

if used == 1 && Y/n < 2 && X/Y > 1.4
    RR = RRsave{2};
elseif used == 2 && X/n < 2 && Y/X > 1.4
    RR = RRsave{1};
end    

if Y/n < 2.5 && X/n < 2.5 && X > 4
    tryboxes = 1;
    rows1 = floor(X/2);
    X1 = 1:rows1;
    X2 = rows1+1:X;
    pop1 = sum(sum(popmap(X1,:)));
    pop2 = sum(sum(popmap(X2,:)));
    n1 = round(n * pop1 / (pop1+pop2));
    n1 = max(min(n1,n-1),1);
    n2  = n - n1;
    
    a1 = a(X1,:);    a2 = a(X2,:);
    b1 = b(X1,:);    b2 = b(X2,:);
    r_1 = basic_alloc(a1,b1,n1);
    r_2 = basic_alloc(a2,b2,n2)+n1;
    r3 = [r_1;r_2];
    RRsave{3} = r3;
else
    tryboxes = 0;
end

r = RRbest;
K = [3.5 2.0  0  0.25  0.6  0.58];  
[RRbest,newminscore] = optRR(a,b,n,RR,minscore,maxtries,X,Y,popmap,ideal_pop,K);
if ~isempty(RRbest)
    r = RRbest;
    minscore = newminscore;
end
if maxtries > 175
    maxtries = 150;
    K = [3.5 2.0  0  0.166  0.4  0.4];  
    [RRbest,newminscore] = optRR(a,b,n,RR,minscore,maxtries,X,Y,popmap,ideal_pop,K);
    if ~isempty(RRbest)
        r = RRbest;
        minscore = newminscore;
    end
end
if tryboxes
    K = [3.5 2.0  0  0.25  0.6  0.58];  
    [RRbest,newminscore] = optRR(a,b,n,RRsave{3},minscore,maxtries,X,Y,popmap,ideal_pop,K);
    if ~isempty(RRbest)
        r = RRbest;
        minscore = newminscore;
    end
    
end


for i=1:n
    rri=rr==i;
    anet(i) = sum(a(rri))-sum(b(rri));
end
if (max(abs(anet))/abs(sum(a(:))/n -sum(b(:))/n) > 1.5)  && (minscore > 3);
    for ii=1:3
        [RRbest, newminscore] = DrinkDrPepper(a,b,n,r,minscore,popmap,ideal_pop);
        if ~isempty(RRbest)
            r = RRbest;
            minscore = newminscore;
        else
            break
        end
    end
end

K = [0.2 0. 07 0.03 0.01 0.003];
for ii =1:5
    [RRbest, newminscore] = fixpop(a,b,n,r,minscore,popmap,ideal_pop, K(ii));
    if ~isempty(RRbest)
        r = RRbest;
        minscore = newminscore;
    end
end

end

function rk = basic_alloc(a,b,n)
N = numel(a);
fixind = zeros(size(a));
fixind(1:end) = 1:length(fixind(:));
fixind(:,2:2:end)=flipud(fixind(:,2:2:end));
Ma=mean(mean(a));
Mb=mean(mean(b));
Mab=(Ma+Mb)/2;
rk = ones(N,1)*(1+n);
totalPop = sum([a(:); b(:)]);
idealDistrictPop = totalPop/n;
in=0;
for k=1:n
    while sum(a(rk<=k))+sum(b(rk<=k))<idealDistrictPop*(k)-Mab
        in=in+1;
        rk(fixind(in))=k;
    end
end
rk(fixind(in+1:N))=n;
rk = reshape(rk,size(a,1),size(a,2));
end

function [RRbest, newminscore] = optRR(a,b,n,RR,minscore,maxtries,X,Y,popmap,ideal_pop,K)

abdiff = a-b;

%repopulate the sums for the best solution
dpop=zeros(1,n); anet = dpop; poperr = dpop;
for i=1:n
    RRi = RR==i;
    dpop(i) = sum(popmap(RRi));
    anet(i) = sum(a(RRi))-sum(b(RRi));
    poperr(i) = dpop(i) - ideal_pop;
end


RRbest = [];
Rpad = 3*n * ones(X+2,Y+2);

done = 0; tries = 0;
movedcnt = zeros(size(RR));
phase = 1;

anet_bias = anet;
anet_bias(1) = anet_bias(1) + anet(2);
anet_bias(n) = anet_bias(n) + anet(n-1);
for i=2:n-1
    anet_bias(i) = anet_bias(i) + anet(i+1)/2 + anet(i-1)/2;
end
anet_bias = anet_bias - mean(anet_bias);
anet_b_gain = 1.0;
cantmove = zeros(size(RR));

while (~done && tries < maxtries)
    
    tries = tries + 1;
    if (tries > maxtries/1.75)
        movedcnt = 0*movedcnt;  phase = 3;
    end
    if (mod(tries,20)==0)
        cantmove = zeros(size(RR));
    end
    Rpad(2:end-1,2:end-1) = RR;
    Rpadup = Rpad; Rpadleft = Rpad; Rpaddown = Rpad; Rpadright = Rpad;
    Rpadup(1:end-1,:)=Rpad(2:end,:);
    Rpadleft(:,1:end-1) = Rpad(:,2:end);
    Rpaddown(2:end,:)=Rpad(1:end-1,:);
    Rpadright(:,2:end) = Rpad(:,1:end-1);

    movable = ones(size(RR));
    movable=movable & (movedcnt < 1) & (popmap > 0.0001) & (cantmove==0);
    movescore = -100*ones(size(RR));
    
    anetmap = zeros(size(RR));
    poperrmap = anetmap;
    anet_x = anet + anet_b_gain * anet_bias;
    anet_wt = ones(size(anet)); 
    
   for jj = 1:n;
        RRjj=RR==jj;
        anetmap(RRjj) = anet_x(jj);
        poperrmap(RRjj) = poperr(jj);
    end
    
    s1 = abs(anetmap) - abs(anetmap - abdiff);
    s2 = abs(poperrmap) - abs(poperrmap - popmap);
    for i=1:n
        grabable = ((Rpadup == i) | (Rpaddown == i) |(Rpadright == i) |(Rpadleft == i)) & (Rpad~=i);
        grabable = grabable(2:end-1,2:end-1);

        if phase > 1,
            gain =K(3);
            anet_b_gain = 0;
        elseif (tries<maxtries*K(4))  %4   
            gain = K(1);
        elseif (tries<maxtries*K(5))  
            gain = K(2);
            anet_b_gain = 0;
        else
            phase = 2;
            gain = K(3);
        end
        
        movescore_k1 = movescore;
        movescore(grabable) = gain * ( abs(anet_x(i)) - abs(anet_x(i) + abdiff(grabable)) + s1(grabable) ) + abs(poperr(i)) - abs(poperr(i) + popmap(grabable)) + s2(grabable) - 0.0*popmap(grabable);
%        movescore(grabable) = gain * ( anet_wt(i) * (abs(anet_x(i)) - abs(anet_x(i) + abdiff(grabable))) + s1(grabable) ) + abs(poperr(i)) - abs(poperr(i) + popmap(grabable)) + s2(grabable) - 0.0*popmap(grabable);
        movescore = max(movescore,movescore_k1);
        
    end
    
    [RRpad,CRpad]=size(Rpad);
    nomove = 1; tt =0;
    while (nomove && tt < 10)
        
        tt = tt+1;  
        movescore(movable==0) = -100;
        [maxsc,loc] = max(movescore(:));
        %disp([tries loc])   
        if maxsc > -40   %%%0
            from = RR(loc);
            [movI,movJ] = ind2sub(size(RR),loc);
            %check for contiguity and which districts could take the piece
            xaround = [movI, movI+1, movI+2, movI+2, movI+2, movI+1, movI, movI];
            yaround = [movJ, movJ, movJ, movJ+1, movJ+2, movJ+2, movJ+2, movJ+1];
            temp = (yaround-1)*RRpad + xaround;
           %temp_o = sub2ind(size(Rpad),xaround,yaround);
            surr_districts = Rpad(temp);
            surr_districts(end+1) = surr_districts(1);
            temp = diff(surr_districts~=from);
            temp = find(temp>0);
            if length(temp)<=1
                %can move without breaking contiguity
                border_districts = surr_districts(2:2:8);
                temp = border_districts(border_districts~=from & border_districts <= n);
                if max(temp) == min(temp);
                    to = max(temp);
                else
                    q = max(temp); qq = min(temp);  %dont worry about a 3rd
                    movescore1 =  ( abs(poperr(q)) - abs(popmap(loc) + poperr(q)) ) + ( abs(anet_x(q)) - abs(anet_x(q) + abdiff(loc)) );
                    movescore2 =  ( abs(poperr(qq)) - abs(popmap(loc) + poperr(qq)) ) + ( abs(anet_x(qq)) - abs(anet_x(qq) + abdiff(loc)) );
                    if movescore1> movescore2,
                        to = q;
                    else
                        to = qq;
                    end
                end
                %Move it!
                RR(loc) = to;
                movedcnt(loc) = movedcnt(loc) + 1;
%                 for i=[to from]
%                     dpop(i) =  sum(popmap(RR==i)); anet(i) = sum(a(RR==i))-sum(b(RR==i));
%                     poperr(i) = dpop(i) - ideal_pop;
%                     %dist_score(i) = abs(poperr(i)) + abs(anet(i))^2;
%                 end
                %incremental evalution of errors
                dpop(to)=dpop(to)+popmap(loc);anet(to)=anet(to)+a(loc)-b(loc);
                poperr(to)=dpop(to)-ideal_pop;
                dpop(from)=dpop(from)-popmap(loc);anet(from)=anet(from)-(a(loc)-b(loc));
                poperr(from)=dpop(from)-ideal_pop;
                
                totscore = sum(abs(poperr)) * max(abs(anet));

                if totscore < minscore,
                    minscore = totscore;
                    RRbest = RR;
                end
                nomove = 0;
            else
                cantmove(loc) = 1;   
            end
        else
            anet_b_gain = 0;
            if phase > 2,
               done = 1;
            elseif phase > 1,
                movedcnt = 0*movedcnt;  phase = 3;
            else
                phase = 2;
            end    
            nomove = 0;
        end
        
    end   % if ~nomove

end    % main loop for swaps
newminscore = minscore;
end

function [RRbest, newminscore] = fixpop(a,b,n,RR,minscore,popmap,ideal_pop, K, dist_list)

RRbest = [];
X = size(a,1);
Y = size(a,2);


Rpad = 3*n * ones(X+2,Y+2);
Rpad(2:end-1,2:end-1) = RR;
Rpadup = Rpad; Rpadleft = Rpad; Rpaddown = Rpad; Rpadright = Rpad;
Rpadup(1:end-1,:)=Rpad(2:end,:);
Rpadleft(:,1:end-1) = Rpad(:,2:end);
Rpaddown(2:end,:)=Rpad(1:end-1,:);
Rpadright(:,2:end) = Rpad(:,1:end-1);

if nargin > 8, 
    special = 1; 
else 
    special = 0;
end

dpop=zeros(1, n); anet = dpop; poperr = dpop;
for i= 1:n
    RRi = (RR==i);
    dpop(i) =  sum(popmap(RRi)); anet(i) = sum(a(RRi))-sum(b(RRi));
    poperr(i) = dpop(i) - ideal_pop;
end

if (~special)
    %disp('****')
    %disp(round(poperr*1000))
    
    temp = diff(sign(poperr));
    best = (temp ~= 0) .* min(abs(poperr(1:end-1)),abs(poperr(2:end)));
    work = find(best > (K/n) *sum(abs(poperr)));
    workval = best(work);
    [junk, temp] = sort(workval,'descend');
    districts = work(temp);
else 
    districts = min(dist_list);
end


useddist=zeros(size(poperr));

for ii=districts,
    
    i = ii;
    j = i+1;
    if max(useddist([i j]))<1
        grabable2 = ((Rpadup == j) | (Rpaddown == j) |(Rpadright == j) |(Rpadleft == j)) & (Rpad==i);
        grabable1 = ((Rpadup == i) | (Rpaddown == i) |(Rpadright == i) |(Rpadleft == i)) & (Rpad==j);
        grabable1 = grabable1(2:end-1,2:end-1);
        grabable2 = grabable2(2:end-1,2:end-1);
        cells1 = find(grabable1);
        cells2 = find(grabable2);
        pop1 = popmap(cells1);
        pop2 = popmap(cells2);
        pop1 = repmat(pop1(:),1,length(pop2));
        pop2 = repmat(pop2(:),1,size(pop1,1))';
        diffpop = pop2 - pop1;
        temp = abs(poperr([i j]));
        if poperr(i) > 0,
            list = find(diffpop>0 & diffpop>min(temp) & diffpop<max(temp));
            list2 = find(diffpop>0 & diffpop<min(temp));
        else
            list = find(diffpop<0 & -diffpop>min(temp) & -diffpop<max(temp));
            list2 = find(diffpop<0 & -diffpop<min(temp));
        end
        %improve = abs(poperr(i)-poperr(j)) - abs(poperr(i) - diffpop) - abs(poperr(2) +diffpop)
        
        list = [list; list2];
        for kk = 1:length(list),
            [swap(1), swap(2)] = ind2sub(size(diffpop),list(kk));
            Q = cells1(swap(1));
            W = cells2(swap(2));
            [iQ, jQ] = ind2sub(size(RR),Q);
            [iW, jW] = ind2sub(size(RR),W);
            temp1 = abs(iQ-iW)+abs(jQ-jW);
            if temp1 > 1
                Yes1 = Check_Contig(Q,RR,Rpad);
                Yes2 = Check_Contig(W,RR,Rpad);
                if Yes1 && Yes2
                    RRtry = RR;
                    RRtry(Q) = RR(W);
                    RRtry(W) = RR(Q);
                    for i= [RRtry(Q) RRtry(W)]              
                        dpop(i) =  sum(popmap(RRtry==i)); anet(i) = sum(a(RRtry==i))-sum(b(RRtry==i));
                        poperr(i) = dpop(i) - ideal_pop;
                    end
                    totscore = sum(abs(poperr)) * max(abs(anet));
                    if totscore < minscore,
                        minscore = totscore;
                        RR = RRtry; 
                        RRbest = RR;
                        useddist(i) = 1; useddist(j) = 1;
                        %disp([round(poperr*1000) round(minscore*10)])
                        break
                    end
                end
            end
        end
    end
end
newminscore = minscore;
end

function [Yes] = Check_Contig(loc,RR,Rpad)
            from = RR(loc);
            [movI,movJ] = ind2sub(size(RR),loc);
            %check for contiguity and which districts could take the piece
            xaround = [movI, movI+1, movI+2, movI+2, movI+2, movI+1, movI, movI];
            yaround = [movJ, movJ, movJ, movJ+1, movJ+2, movJ+2, movJ+2, movJ+1];
            temp = sub2ind(size(Rpad),xaround,yaround);
            surr_districts = Rpad(temp);
            surr_districts(end+1) = surr_districts(1);
            temp = diff(surr_districts~=from);
            temp = find(temp>0);
            if length(temp)<=1
                %can move without breaking contiguity
                Yes = 1;
            else
                Yes = 0;
            end
end

function [RRbest, newminscore] = DrinkDrPepper(a,b,n,RR,minscore,popmap,ideal_pop)

RRbest = [];
X = size(a,1);
Y = size(a,2);

Rpad = 3*n * ones(X+2,Y+2);
Rpad(2:end-1,2:end-1) = RR;
Rpadup = Rpad; Rpadleft = Rpad; Rpaddown = Rpad; Rpadright = Rpad;
Rpadup(1:end-1,:)=Rpad(2:end,:);
Rpadleft(:,1:end-1) = Rpad(:,2:end);
Rpaddown(2:end,:)=Rpad(1:end-1,:);
Rpadright(:,2:end) = Rpad(:,1:end-1);

dpop=zeros(1,n); anet =dpop; poperr = dpop;
for i= 1:n                   
    dpop(i) =  sum(popmap(RR==i)); anet(i) = sum(a(RR==i))-sum(b(RR==i));
    poperr(i) = dpop(i) - ideal_pop;
end

[~,wdist] = max(abs(anet));
xdist_list = [wdist-1 wdist+1];
xdist_list = xdist_list(xdist_list > 0 & xdist_list<= n);
if length(xdist_list) > 1   
    check1 = abs(anet(wdist) - anet(xdist_list(1)));
    check2 = abs(anet(wdist) - anet(xdist_list(2)));
    if check2 > check1
        xdist_list = [xdist_list(2) xdist_list(1)];
    end
end

abdiff = a-b;
populated = popmap > 0;    

if anet(wdist) < 1;
     abdiff = -abdiff;
end

done = 0;
for style = 1:2
    if (done == 1)
        break
    end
    for xdist = xdist_list
        
        if (done == 1)
            break
        end
        achange_tgt = abs(anet(wdist) - anet(xdist))/2;
        
        i = wdist;
        j = xdist;
        pei = poperr(i); pej = poperr(j);
        if sign(pei)==sign(pej)
            if pei < pej
                popdiff_allow_pos = max(abs(pei), abs(pej));
                popdiff_allow_neg = -1*min(abs(pei), abs(pej));
            else
                popdiff_allow_pos = min(abs(pei), abs(pej));
                popdiff_allow_neg = -1*max(abs(pei), abs(pej));
            end
            
        else
            temp = sign(pei) * max(abs(pei), abs(pej));
            popdiff_allow_pos = max( temp , 0 );
            popdiff_allow_neg = max( temp , 0 );
        end
        
        grabable1 = ((Rpadup == j) | (Rpaddown == j) |(Rpadright == j) |(Rpadleft == j)) & (Rpad==i);
        grabable2 = ((Rpadup == i) | (Rpaddown == i) |(Rpadright == i) |(Rpadleft == i)) & (Rpad==j);
        grabable1 = grabable1(2:end-1,2:end-1) & populated;
        grabable2 = grabable2(2:end-1,2:end-1) & populated;
        cells1 = find(grabable1);
        cells2 = find(grabable2);
        
        if ~isempty(cells1) && ~isempty(cells2);
            
            abdiff1 = abdiff(cells1);
            abdiff2 = abdiff(cells2);
            abdiff1 = repmat(abdiff1(:),1,length(abdiff2));
            abdiff2 = repmat(abdiff2(:),1,size(abdiff1,1))';
            achange = abdiff1 - abdiff2;
            
            achangeerr = achange - achange_tgt;
            achangeerr(achange > achange_tgt*2  | achange < 0) = 999;
            
            if style==1
                pop1 = popmap(cells1);
                pop2 = popmap(cells2);
                pop1 = repmat(pop1(:),1,length(pop2));
                pop2 = repmat(pop2(:),1,size(pop1,1))';
                diffpop = pop1 - pop2;
                achangeerr(diffpop >  popdiff_allow_pos) = 999;
                achangeerr(diffpop <  popdiff_allow_neg) = 999;
            else
                diffpop = zeros(size(achange));
            end
            
            
            
            listtemp = find (achangeerr(:) < 99);
            if ~isempty(listtemp)
                ace_list = abs(achangeerr(listtemp));
                ace_list = ace_list(1:min(10,end));
                [~,tt] = sort(ace_list);
                list = listtemp(tt);
                
                for kk = 1:length(list),
                    [swap(1), swap(2)] = ind2sub(size(diffpop),list(kk));
                    Q = cells1(swap(1));
                    W = cells2(swap(2));
                    [iQ, jQ] = ind2sub(size(RR),Q);
                    [iW, jW] = ind2sub(size(RR),W);
                    temp1 = abs(iQ-iW)+abs(jQ-jW);
                    if temp1 > 1
                        Yes1 = Check_Contig(Q,RR,Rpad);
                        Yes2 = Check_Contig(W,RR,Rpad);
                        if Yes1 && Yes2
                            if style == 1,
                                temp = RR(Q);
                                RR(Q) = RR(W);
                                RR(W) = temp;
                                for i= [RR(Q) RR(W)]
                                    dpop(i) =  sum(popmap(RR==i)); anet(i) = sum(a(RR==i))-sum(b(RR==i));
                                    poperr(i) = dpop(i) - ideal_pop;
                                end
                                totscore = sum(abs(poperr)) * max(abs(anet));
                                if totscore < minscore,
                                    minscore = totscore;
                                    RRbest = RR;
                                end
                                done = 1;
                                break
                            else
                                RRtest = RR;
                                RRtest(Q) = RR(W);
                                RRtest(W) = RR(Q);
                                [RRtest2, newminscore] = fixpop(a,b,n,RRtest,minscore,popmap,ideal_pop, 0, [RR(W) RR(Q)]);
                                if ~isempty(RRtest2) && (newminscore < minscore)
                                    RRbest = RRtest2;
                                    minscore = newminscore;
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

newminscore = minscore;
end