nrow = randi([3,10])+3;
ncol = randi([3,10]);
n = randi([2,8])

a = randi(40, nrow, ncol);
b = randi(40, nrow, ncol);

ne = nrow*ncol
t = a+b
tot = sum(t(:))

% if nrow>ncol %tall and skinny
%     a = a';
%     b = b';
% end
% r = 0*a;

%function r = solver(a,b,n)
    [nrow,ncol] = size(a);
    r = a*0+n;
    assert(n<20)
    if nrow>n
        r=floor(linspace(1,n+0.99,ncol))'*ones(1,ncol);
    elseif ncol>n
        r=ones(nrow,1)*floor(linspace(1,n+0.99,nrow));
    else
        r(1:n) = 1:n;
    end
%end

