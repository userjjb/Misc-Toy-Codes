function r = solver3(a,b,n)
    temp = floor(linspace(1,n+1,numel(a)+1));
    r = reshape(temp(1:end-1),size(a));
    r(:,2:2:end)=flipud(r(:,2:2:end));
    [y,x] = size(a);
    
    if n==2
        pop = a+b;
        mask = bsxfun(@eq,r,reshape(1:n,[1,1,n]));
        districtPop = squeeze(sum(sum(a.*mask)) + sum(sum(b.*mask)));

        totalPop = sum([a(:); b(:)]);
        idealDistrictPop = totalPop/n;
        penalty1 = sum(abs(districtPop - idealDistrictPop));
        
        A = toeplitz([[1 1] zeros(1,x-2)]);
        B = toeplitz([[0 1] zeros(1,y-2)]);
        basis = ones(size(a));
        basis_avg = (B*basis + basis*A);
        r_avg = (B*r + r*A)./basis_avg;
        border = r_avg~=floor(r_avg);
        
        [~, index] = max(districtPop);
        swap = find(and(border, r==index));
        maxswap = 20;
        if numel(swap)>maxswap
            [~,i] = sort(pop(swap));
            swap = swap(i(1:maxswap));
        end
        
        if min(pop(swap)) < penalty1/2
            combos = dec2bin(0:2^numel(swap)-1)-48;
            [~, index2] = min(abs(combos*pop(swap) - penalty1/2));
            r(swap(combos(index2,:)==1)) = not(index-1)+1;
        end
    end
end