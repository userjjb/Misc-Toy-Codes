function r = solver2(a,b,n)
    temp = floor(linspace(1,n+1,numel(a)+1));
    r = reshape(temp(1:end-1),size(a));
    r(:,2:2:end)=flipud(r(:,2:2:end));
end