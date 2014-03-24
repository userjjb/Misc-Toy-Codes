%Calculates probability of rolling a value above a fractional limit of the
%maximum dice roll. Accepts threshold values below which the chance of
%"success" depends the value rolled.
function chance = above_limit_roll(s,n,frac_limit,varargin)
    if size(varargin,2)>0
        crit_threshold= cell2mat(varargin(1));
    else
        crit_threshold= -1;
    end
    
    max_roll= s*n;
    limit= ceil(frac_limit*max_roll);
    sum= 0;
    for k=limit:max_roll
        if k<crit_threshold
            mult= min((k+1)/100,1);
        else
            mult=1;
        end
        sum= sum+(mult*roll_binom(s,n,k));
    end
    chance= sum;
end

%Uses recursive method to calculate probability of rolling a 'k' with 
%'n' d 's' dice
function chance = roll_r(s,n,k)
    if n==1
        if k>=1 && k<=s
            chance= 1/s;
            else
            chance= 0;
        end
        return
    end
    sum= 0;
    %Heuristic to "split the difference" on the two recursion branches
    x= ceil(n/2);
    y= n-x;
    for i=1:k-n+1
            sum= sum + (roll_r(s,x,i)*roll_r(s,y,k-i));
    end
    chance= sum;
end

%uses binomial method to calculate probability of rolling a 'k' with 'n' d
%'s' dice
function chance = roll_binom(s,n,k)
    chance= 1/(s^n);
    sum=0;
    for i=0:floor((k-n)/s)
        sum= sum+ (((-1)^i) * binom_coef(n,i) * binom_coef(k-s*i-1,n-1));
    end
    chance= chance*sum;
end

%Uses multiplicative formula to calculate binomial coefficient
function coef = binom_coef(n,k)
    product=1;
    for i=1:k
        product= product*((n-(k-i))/i);
    end
    coef= product;
end