n = 17;
N = reshape(1:n^2,n,n);
order = [];
for i=1:(n-1)/2
    if mod(i,2)
        order = [order, N(n-(i-1):-1:1+i, i)', N(i, 1+(i-1):n-i)];
    else
        order = [order, N(1+(i-1):n-i, n-(i-1))', N(n-(i-1), n-(i-1):-1:1+i)];
    end
end

order = [order, (n^2+1)/2];

for i=(n-1)/2:-1:1
    if mod(i,2)
        order = [order, N(n-(i-1), 1+i:n-(i-1)), N(n-i:-1:1+(i-1), n-(i-1))'];
    else
        order = [order, N(i, n-i:-1:1+(i-1)), N(1+i:n-(i-1), i)'];
    end
end

milk = zeros(1,n^2);
milk(order)=N;
milk = reshape(milk,n,n)-(n^2+1)/2;

inner_corner = tril(ones(13))*[1,1,3,3,5,5,7,7,9,9,11,11,13]';

inner=[1];
outer=[4];
for i=1:length(inner_corner)-1
    inner = [inner, inner_corner(i):inner_corner(i+1)];
    outer = [outer, outer(end)+2:outer(end)+2+(inner_corner(i+1)-inner_corner(i))];
end