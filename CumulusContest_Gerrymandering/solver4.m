function r = solver4(a,b,n)

N = numel(a);
maxtries = ceil(min(N/1.10,450));
Y = size(a,2);
X = size(a,1);
popmap = b+a;
apop = sum(a(:));
bpop = sum(b(:));
ideal_a = apop/n;
ideal_b = bpop/n;
ideal_pop = ideal_a + ideal_b;

RRsave{1} = basic_alloc(a,b,n);
RRsave{2} = basic_alloc(a',b',n)';

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
K = [3.1 1.7  0  0.1 0.2 0.2];  
[RRbest,newminscore] = optRR(a,b,n,RR,minscore,maxtries,X,Y,popmap,ideal_pop,K);
if ~isempty(RRbest)
    r = RRbest;
    minscore = newminscore;
end
if maxtries > 175
    maxtries = 100;
    [RRbest,newminscore] = optRR(a,b,n,RR,minscore,maxtries,X,Y,popmap,ideal_pop,K);
    if ~isempty(RRbest)
        r = RRbest;
        minscore = newminscore;
    end
end
if tryboxes
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
if (max(abs(anet))/abs(sum(a(:))/n -sum(b(:))/n) > 1.3)  && (minscore > 2);
    for ii=1:4
        [RRbest, newminscore] = DrinkDrPepper(a,b,n,r,minscore,popmap,ideal_pop);
        if ~isempty(RRbest)
            r = RRbest;
            minscore = newminscore;
        else
            break
        end
    end
end

K = [0.2 0.01 0.005 0.005 0.003];
for ii =1:5
    [RRbest, newminscore] = fixpop(a,b,n,r,minscore,popmap,ideal_pop, K(ii));
    if ~isempty(RRbest)
        r = RRbest;
        minscore = newminscore;
    end
end

r = improve_solution(r,a,b,n);

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
ab = a+b;
for k=1:n
    while sum(ab(rk<=k))<idealDistrictPop*(k)-Mab
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

done = false;
tries = 0;
movedcnt = zeros(X,Y);
phase = 1;

anet_bias = anet;
anet_bias(1) = anet_bias(1) + anet(2);
anet_bias(n) = anet_bias(n) + anet(n-1);
for i=2:n-1
    anet_bias(i) = anet_bias(i) + anet(i+1)/2 + anet(i-1)/2;
end
anet_bias = anet_bias - mean(anet_bias);
anet_b_gain = 1.0;
cantmove = zeros(X,Y);

while (~done && tries < maxtries)
    
    tries = tries + 1;
    if (tries > maxtries/1.75)
        movedcnt = 0*movedcnt;  phase = 3;
    end
    if (mod(tries,20)==0)
        cantmove = zeros(X,Y);
    end

    Rpad(2:end-1,2:end-1) = RR;
    Rpadup = Rpad; Rpadleft = Rpad; Rpaddown = Rpad; Rpadright = Rpad;
    Rpadup(1:end-1,:)=Rpad(2:end,:);
    Rpadleft(:,1:end-1) = Rpad(:,2:end);
    Rpaddown(2:end,:)=Rpad(1:end-1,:);
    Rpadright(:,2:end) = Rpad(:,1:end-1);

    movable   = ones(X,Y);
    movescore = -100*movable;
    movable   = movable & (movedcnt < 1) & (popmap > 0.0001) & (cantmove==0);
    anetmap   = zeros(X,Y);    
    poperrmap = anetmap;
    anet_x = anet + anet_b_gain * anet_bias;
    
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
if gain == 0
    movescore(grabable) = abs(poperr(i)) - abs(poperr(i) + popmap(grabable)) + s2(grabable);
else
    movescore(grabable) = gain * ( abs(anet_x(i)) - abs(anet_x(i) + abdiff(grabable)) + s1(grabable) ) + abs(poperr(i)) - abs(poperr(i) + popmap(grabable)) + s2(grabable);
end
        movescore = max(movescore,movescore_k1);
        
    end
    
    nomove = 1; tt =0;
    while (nomove && tt < 10)
        
        tt = tt+1;  
        movescore(movable==0) = -100;
        [maxsc,loc] = max(movescore(:));
        %disp([tries loc])   
        if maxsc > -40   %%%0
            from = RR(loc);
            movJ = floor((loc-1)/X) + 1;
            movI = loc - (movJ-1)*X;           

            %check for contiguity and which districts could take the piece
            xaround = [movI, movI+1, movI+2, movI+2, movI+2, movI+1, movI, movI];
            yaround = [movJ, movJ, movJ, movJ+1, movJ+2, movJ+2, movJ+2, movJ+1];
            temp = (yaround-1)*(X+2) + xaround;
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
               done = true;
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
    temp = diff(sign(poperr));
    best = (temp ~= 0) .* min(abs(poperr(1:end-1)),abs(poperr(2:end)));
    work = find(best > (K/n) *sum(abs(poperr)));
    workval = best(work);
    [~, temp] = sort(workval,'descend');
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
        Xp   = length(pop1);
        Yp   = length(pop2);
        pop1 = pop1(:)*ones(1,Yp);
        pop2 = ones(Xp,1)*pop2(:)';
        diffpop = pop2 - pop1;
        temp = abs(poperr([i j]));
        if poperr(i) > 0,
            list1 = find(diffpop>0 & diffpop>min(temp) & diffpop<max(temp));
            list2 = find(diffpop>0 & diffpop<min(temp));
        else
            list1 = find(diffpop<0 & -diffpop>min(temp) & -diffpop<max(temp));
            list2 = find(diffpop<0 & -diffpop<min(temp));
        end
        
        list = [list1; list2];
        
        for kk = 1:length(list),
            
            swap(2) = floor((list(kk)-1)/Xp) + 1;
            swap(1) = list(kk) - (swap(2)-1)*Xp;
            Q = cells1(swap(1));
            W = cells2(swap(2));
            jQ = floor((Q-1)/X) + 1;
            iQ = Q - (jQ-1)*X;
            jW = floor((W-1)/X) + 1;
            iW = W - (jW-1)*X;

            temp1 = abs(iQ-iW)+abs(jQ-jW);
            if temp1 > 1
                Yes1 = Check_Contig(Q,RR,X,Rpad);
                Yes2 = Check_Contig(W,RR,X,Rpad);
                if Yes1 && Yes2
                    RRtry = RR;
                    RRtry(Q) = RR(W);
                    RRtry(W) = RR(Q);

                    to = RRtry(Q);
                    from = RRtry(W);
                    xdpop = dpop; xanet = anet; xpoperr = poperr;
                    xdpop(to)=xdpop(to)+popmap(Q)-popmap(W); xanet(to)=xanet(to)+(a(Q)-b(Q)-a(W)+b(W));
                    xpoperr(to)=xdpop(to)-ideal_pop;
                    xdpop(from)=xdpop(from)-popmap(Q)+popmap(W); xanet(from)=xanet(from)-(a(Q)-b(Q)-a(W)+b(W));
                    xpoperr(from)=xdpop(from)-ideal_pop;
                     
                    totscore = sum(abs(xpoperr)) * max(abs(xanet));
                    if totscore < minscore,
                        minscore = totscore;
                        RR = RRtry; 
                        dpop = xdpop; anet = xanet; poperr = xpoperr;
                        RRbest = RR;
                        useddist(i) = 1; useddist(j) = 1;
                        break
                    end
                end
            end
        end
    end
end
newminscore = minscore;
end

function [Yes] = Check_Contig(loc,RR,X,Rpad)
            from = RR(loc);
            movJ = floor((loc-1)/X) + 1;
            movI = loc - (movJ-1)*X;
            
            %check for contiguity and which districts could take the piece
            xaround = [movI, movI+1, movI+2, movI+2, movI+2, movI+1, movI, movI];
            yaround = [movJ, movJ, movJ, movJ+1, movJ+2, movJ+2, movJ+2, movJ+1];
            temp    = (yaround-1)*(X+2) + xaround;
            surr_districts = Rpad(temp);
            surr_districts(end+1) = surr_districts(1);
            temp = diff(surr_districts==from);
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

done = false;
for style = 1:2
    if (done)
        break
    end
    for xdist = xdist_list
        
        if (done)
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
            Xp      = length(abdiff1);
            Yp      = length(abdiff2);
            abdiff1 = abdiff1(:)*ones(1,Yp);
            abdiff2 = ones(Xp,1)*abdiff2(:)';                           
            achange = abdiff1 - abdiff2;
            
            achangeerr = achange - achange_tgt;
            achangeerr(achange > achange_tgt*2  | achange < 0) = 999;
            
            if style==1
                pop1 = popmap(cells1);
                pop2 = popmap(cells2);
                Xp   = length(pop1);
                Yp   = length(pop2);
                pop1 = pop1(:)*ones(1,Yp);
                pop2 = ones(Xp,1)*pop2(:)';
                diffpop = pop1 - pop2;
                achangeerr(diffpop >  popdiff_allow_pos) = 999;
                achangeerr(diffpop <  popdiff_allow_neg) = 999;
            end
            
            listtemp = find (achangeerr(:) < 99);
            if ~isempty(listtemp)
                ace_list = abs(achangeerr(listtemp));
                ace_list = ace_list(1:min(10,end));
                [~,tt] = sort(ace_list);
                list = listtemp(tt);
                
                for kk = 1:length(list),
                                        
                    swap(2) = floor((list(kk)-1)/Xp) + 1;
                    swap(1) = list(kk) - (swap(2)-1)*Xp;
                    
                    Q = cells1(swap(1));
                    W = cells2(swap(2));
                    
                    jQ = floor((Q-1)/X) + 1;
                    iQ = Q - (jQ-1)*X;
                    jW = floor((W-1)/X) + 1;
                    iW = W - (jW-1)*X;

                    temp1 = abs(iQ-iW)+abs(jQ-jW);
                    if temp1 > 1
                        Yes1 = Check_Contig(Q,RR,X,Rpad);
                        Yes2 = Check_Contig(W,RR,X,Rpad);
                        if Yes1 && Yes2
                            if style == 1,
                                temp = RR(Q);
                                RR(Q) = RR(W);
                                RR(W) = temp;

                                to = RR(Q);
                                from = RR(W);
                                xdpop = dpop; xanet = anet; xpoperr = poperr;
                                xdpop(to)=xdpop(to)+popmap(Q)-popmap(W); xanet(to)=xanet(to)+(a(Q)-b(Q)-a(W)+b(W));
                                xpoperr(to)=xdpop(to)-ideal_pop;
                                xdpop(from)=xdpop(from)-popmap(Q)+popmap(W); xanet(from)=xanet(from)-(a(Q)-b(Q)-a(W)+b(W));
                                xpoperr(from)=xdpop(from)-ideal_pop;
                     
                                totscore = sum(abs(xpoperr)) * max(abs(xanet));
                                if totscore < minscore,
                                    minscore = totscore;
                                    RRbest = RR;
                                end
                                done = true;
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


function [r,penalty] = improve_solution(r,a,b,N)

% Size of Rectanglia
[Rows,Cols] = size(a);

% Population and disparity of Rectanglia
POP  = a+b;
DIS  = a-b;
Z    = a*0;
U    = Z + 1;

% Identify Regions and Frontier
[v,L,F] = get_labels(r,1,Z,N);

if ~v
    r        = U;
    r(1:N-1) = N:-1:2;
end

% Calculate District Metrics
total_pop = sum(POP(:));
ideal_pop = total_pop/N;

POP_ERROR = zeros(N,1);
DIS_ERROR = zeros(N,1);
NUM_ELEMS = zeros(N,1);

for d = 1:N
    
    % Distric indexes
    idx_d    = r==d;
    
    % Calculate District Metrics
    POP_d    = POP(idx_d);
    DIS_d    = DIS(idx_d);
    
    POP_ERROR(d) = sum(POP_d);
    DIS_ERROR(d) = sum(DIS_d);
    NUM_ELEMS(d) = sum(idx_d(:));
    
end

penalties = [sum(abs(POP_ERROR - ideal_pop)); max(abs(DIS_ERROR))];
penalty   = penalties(1)*penalties(2);

if ~v
    return;
end

ch    = 1;
it    = 1;

while ch > 0 && it < 20
        
    ch = 0;
    
    [r_border, c_border] = find(F);
    
    % For all border points
    for m = 1:length(r_border)
                
        % Get Connected Districts
        r_m    = r_border(m);
        c_m    = c_border(m);
        idx_m  = (c_m-1)*Rows + r_m;
        
        if F(idx_m)
            
            % Get Connected Districts
            d_m    = r(idx_m);
            l_m    = L(idx_m);
            pe_m   = POP_ERROR(d_m);
            de_m   = DIS_ERROR(d_m);
            
            [~, r_con, c_con] = get_connections(r,L,Rows,Cols,r_m,c_m,d_m,0);
            
            % For all connected Districts
            for f = 1:length(r_con)
                
                r_f   = r_con(f);
                c_f   = c_con(f);
                idx_f = (c_f-1)*Rows + r_f;
                d_f   = r(idx_f);
                l_f   = L(idx_f);
                pe_f  = POP_ERROR(d_f);
                de_f  = DIS_ERROR(d_f);
                
                pe_mf     = abs(pe_m - ideal_pop) + abs(pe_f - ideal_pop);
                pe_mf_new = abs(pe_m + POP(idx_f) - ideal_pop) + abs(pe_f - POP(idx_f) - ideal_pop);
                
                de_m_new  = abs(de_m + DIS(idx_f));
                de_f_new  = abs(de_f - DIS(idx_f));
                
                % Test if is worth to join the point f to the distric m
                cond_1 = (pe_mf_new < pe_mf && de_m_new <= penalties(2) && de_f_new <= penalties(2));
                                
                if NUM_ELEMS(d_f) > 1 && (cond_1 || get_score_diff(DIS_ERROR,DIS(idx_f),d_m,d_f,pe_mf,pe_mf_new,penalties) < 0)

                    % Get neightbours
                    [~, r_con_f, c_con_f] = get_connections(r,Z,Rows,Cols,r_f,c_f,d_f,1);
                    
                    % Update District
                    r(idx_f) = d_m;
                    L(idx_f) = l_m;
                    
                    
                    % Check isolation of neightbours
                    rn              = r(min(r_con_f):max(r_con_f),min(c_con_f):max(c_con_f));
                    valid_isolation = get_labels(rn,d_f,-1*(rn~=d_f),1);
                    
                    % if isolation
                    if ~valid_isolation
                        
                        r(idx_f) = d_f;
                        L(idx_f) = l_f;
                        
                    else
                        
                        % Update Errors and Element
                        POP_ERROR(d_m) = POP_ERROR(d_m) + POP(idx_f);
                        POP_ERROR(d_f) = POP_ERROR(d_f) - POP(idx_f);
                        
                        DIS_ERROR(d_m) = DIS_ERROR(d_m) + DIS(idx_f);
                        DIS_ERROR(d_f) = DIS_ERROR(d_f) - DIS(idx_f);
                        
                        NUM_ELEMS(d_m) = NUM_ELEMS(d_m) + 1;
                        NUM_ELEMS(d_f) = NUM_ELEMS(d_f) - 1;
                        
                        penalties = [sum(abs(POP_ERROR - ideal_pop)); max(abs(DIS_ERROR))];

                        % Update Frontiers
                        [~, r_con_f, c_con_f] = get_connections(r,Z,Rows,Cols,r_f,c_f,d_f,2);
                        [~, r_con_m, c_con_m] = get_connections(r,Z,Rows,Cols,r_m,c_m,d_m,2);
                        
                        r_front = [r_con_m; r_con_f];
                        c_front = [c_con_m; c_con_f];
                        for fn = 1:length(r_front)
                            r_fn      = r_front(fn);
                            c_fn      = c_front(fn);
                            idx_fn    = (c_fn-1)*Rows + r_fn;
                            d_fn      = r(idx_fn);
                            F(idx_fn) = get_connections(r,L,Rows,Cols,r_fn,c_fn,d_fn,0);
                        end
                        
                        ch = ch + 1;
                        break;
                        
                    end
                    
                end
                
            end
        end
        
    end
    
    it = it + 1;
    
end

penalty = penalties(1)*penalties(2);

end

function Sd = get_score_diff(DIS_ERROR,DIS_f,d_m,d_f,pe_mf,pe_mf_new,penalties)

DIS_ERROR(d_m) = DIS_ERROR(d_m) + DIS_f;
DIS_ERROR(d_f) = DIS_ERROR(d_f) - DIS_f;

penalty_ini    = penalties(1)*penalties(2);
penalties(1)   = penalties(1) - pe_mf + pe_mf_new;
penalties(2)   = max(abs(DIS_ERROR));
penalty_end    = penalties(1)*penalties(2);

Sd = penalty_end - penalty_ini;

end

function [v,L,F,k] = get_labels(D,d,L,n)

[r_adj,c_adj] = find(D==d,1,'first');
[R,C]         = size(D);
F             = zeros(R,C);
k             = 1;

while ~isempty(r_adj)
    
    r     = r_adj(1);
    c     = c_adj(1);
    idx   = (c-1)*R + r;
    r_adj = r_adj(2:end);
    c_adj = c_adj(2:end);
          
    if L(idx) == 0
        
        [is_fr, r_con, c_con] = get_connections(D,L,R,C,r,c,d,1);
        
        F(r,c) = is_fr;
        L(r,c) = k;
        r_adj  = [r_adj; r_con]; %#ok<AGROW>
        c_adj  = [c_adj; c_con]; %#ok<AGROW>
                
    end
    
    if isempty(r_adj) && any(any(L==0))
        d = d + 1;
        [r_adj,c_adj] = find(D==d,1,'first');
        if ~isempty(r_adj)
            k = k + 1;
        else
            [r_adj,c_adj] = find(L==0,1,'first');
            d = D(r_adj,c_adj);
            if ~isempty(r_adj)
                k = k + 1;
            end
        end
    end
    
end

v = k == n;

end

function [f, r_adj, c_adj] = get_connections(D,L,R,C,r,c,d,type)

if r > 1 && r < R && c > 1 && c < C
    d_adj = [D(r+1,c); D(r-1,c); D(r,c-1); D(r,c+1)];
    l_adj = [L(r+1,c); L(r-1,c); L(r,c-1); L(r,c+1)];
    r_adj = [r+1; r-1; r  ; r  ];
    c_adj = [c  ; c  ; c-1; c+1];
elseif r == 1 && r < R && c == 1 && c < C
    d_adj = [D(r+1,c); D(r,c+1)];
    l_adj = [L(r+1,c); L(r,c+1)];
    r_adj = [r+1; r  ];
    c_adj = [c  ; c+1];
elseif r == 1 && r < R && c > 1 && c < C
    d_adj = [D(r+1,c); D(r,c-1); D(r,c+1)];
    l_adj = [L(r+1,c); L(r,c-1); L(r,c+1)];
    r_adj = [r+1; r  ; r  ];
    c_adj = [c  ; c-1; c+1];
elseif r == 1 && r < R && c > 1 && c == C
    d_adj = [D(r+1,c); D(r,c-1)];
    l_adj = [L(r+1,c); L(r,c-1)];
    r_adj = [r+1; r  ];
    c_adj = [c  ; c-1];
elseif r > 1 && r < R && c == 1 && c < C
    d_adj = [D(r+1,c); D(r-1,c); D(r,c+1)];
    l_adj = [L(r+1,c); L(r-1,c); L(r,c+1)];
    r_adj = [r+1; r-1; r  ];
    c_adj = [c  ; c  ; c+1];
elseif r > 1 && r < R && c > 1 && c == C
    d_adj = [D(r+1,c); D(r-1,c); D(r,c-1)];
    l_adj = [L(r+1,c); L(r-1,c); L(r,c-1)];
    r_adj = [r+1; r-1; r  ];
    c_adj = [c  ; c  ; c-1];
elseif r > 1 && r == R && c == 1 && c < C
    d_adj = [D(r-1,c); D(r,c+1)];
    l_adj = [L(r-1,c); L(r,c+1)];
    r_adj = [r-1; r  ];
    c_adj = [c  ; c+1];
elseif r > 1 && r == R && c > 1 && c < C
    d_adj = [D(r-1,c); D(r,c-1); D(r,c+1)];
    l_adj = [L(r-1,c); L(r,c-1); L(r,c+1)];
    r_adj = [r-1; r  ; r  ];
    c_adj = [c  ; c-1; c+1];
elseif r > 1 && r == R && c > 1 && c == C
    d_adj = [D(r-1,c); D(r,c-1)];
    l_adj = [L(r-1,c); L(r,c-1)];
    r_adj = [r-1; r  ];
    c_adj = [c  ; c-1];
elseif r == 1 && c == 1 && c < C
    d_adj = D(r,c+1);
    l_adj = L(r,c+1);
    r_adj = r;
    c_adj = c+1;
elseif r == 1 && c > 1 && c < C
    d_adj = [D(r,c-1); D(r,c+1)];
    l_adj = [L(r,c-1); L(r,c+1)];
    r_adj = [r  ; r  ];
    c_adj = [c-1; c+1];
elseif r == 1 && c > 1 && c == C
    d_adj = D(r,c-1);
    l_adj = L(r,c-1);
    r_adj = r;
    c_adj = c-1;
elseif r == 1 && r < R && c == 1
    d_adj = D(r+1,c);
    l_adj = L(r+1,c);
    r_adj = r+1;
    c_adj = c;
elseif r > 1  && r < R && c == 1
    d_adj = [D(r+1,c); D(r-1,c)];
    l_adj = [L(r+1,c); L(r-1,c)];
    r_adj = [r+1; r-1];
    c_adj = [c  ; c  ];
elseif r > 1  && r == R && c == 1
    d_adj = D(r-1,c);
    l_adj = L(r-1,c);
    r_adj = r-1;
    c_adj = c;
else
    d_adj = D(r,c);
    l_adj = 0;
    r_adj = r;
    c_adj = c;
end

idx_c = d == d_adj;
nc    = sum(idx_c);
lc    = length(idx_c);

switch type
    
    case 0
        
        idx_v = ~idx_c;        
        f     = nc ~= lc;
        r_adj = r_adj(idx_v);
        c_adj = c_adj(idx_v);
 
    case 1
        
        idx_v = l_adj == 0 & idx_c;
        
        c     = nc > 0 || sum(sum(D==d))==1;
        f     = nc ~= lc;
                
        if c 
            r_adj = r_adj(idx_v);
            c_adj = c_adj(idx_v);
        else
            r_adj = [];
            c_adj = [];
        end
        
    case 2
        
        idx_v = l_adj == 0;
        
        c     = nc > 0 || sum(sum(D==d))==1;
        f     = nc ~= lc;
        
        if c
            r_adj = r_adj(idx_v);
            c_adj = c_adj(idx_v);
        else
            r_adj = [];
            c_adj = [];
        end        
        
end



end