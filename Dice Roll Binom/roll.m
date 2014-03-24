%Calculates probability of rolling a 'k' with 'n' d 's' dice
function chance = roll(s,n,k)
    if n==1
        if 1<=k<=s
            chance= 1/s;
        else
            chance= 0;
        end
        return
    end
    sum= 0;
    for i=1:k-n+1
            sum= sum + (roll(s,1,i)*roll(s,n-1,k-i));
    end
    chance= sum;
end

%Calculates probability of rolling a value above a fractional limit of the
%maximum dice roll
function chance = above_limit(s,n,frac_limit)
    max_roll= s*n;
    limit= ceil(frac_limit*max_roll);
    sum= 0;
    for k=limit:max_roll
        sum= sum+roll(s,n,k);
    end
    chance= sum;
end