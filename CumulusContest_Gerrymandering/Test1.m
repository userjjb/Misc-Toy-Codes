clearvars
rng(412)
h = waitbar(0,'Please wait...');

results = zeros(640,1);
count = 0;
tic
for count=1:640
    nrow = randi([3,20]);
    ncol = randi([3,20]);
    n = randi([2,10]);
    if randi(2)==2
        n = randi([10,30]);
    end

    a = randi(40, nrow, ncol);
    b = randi(40, nrow, ncol);
    n = max(2, floor(min(n,numel(a)/5)));
    
    r = solver4(a,b,n);
    
    mask = bsxfun(@eq,r,reshape(1:n,[1,1,n]));
    districtPop = squeeze(sum(sum(a.*mask)) + sum(sum(b.*mask)));

    totalPop = sum([a(:); b(:)]);
    idealDistrictPop = totalPop/n;
    penalty1 = sum(abs(districtPop - idealDistrictPop));
    
    penalty2 = max(abs(squeeze(sum(sum(a.*mask)) - sum(sum(b.*mask)))));

    results(count) = penalty1 * penalty2;
    waitbar(count/640, h, sum(results))
end
close(h)
sum(results)/640
toc