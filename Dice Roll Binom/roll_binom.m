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