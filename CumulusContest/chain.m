function b = chain(c,n)
    b = chain3(c,n);
end

%JJB#1 - Naive first test function
function b = chain1(c,n)
    switchback = bsxfun(@plus,n*repmat([1:n-2; n-2:-1:1]',1,(n-1)/2),1:n-1);
    Order = [n:-1:1,switchback(:)',n*[2:n-1],n^2:-1:n^2-(n-1)]; %"Ladder" ordering
    centered_start = floor((n^2-length(c))/2)+1;
    b = zeros(n);
    b(Order(centered_start:centered_start+length(c)-1)) = c;
end

%JJB#2 - Naive first test function
function b = chain2(c,n)
    switchback = bsxfun(@plus,n*repmat([1:n-2; n-2:-1:1]',1,(n-1)/2),1:n-1);
    Order = [n:-1:1,switchback(:)',n*[2:n-1],n^2:-1:n^2-(n-1)]; %"Ladder" ordering
    centered_start = floor((n^2-length(c))/2)+1;
    b = zeros(n);
    b(Order(centered_start:centered_start+length(c)-1)) = 1:length(c);
    for j=1:10
        randi(1000,1000)*randi(1000,1000);
    end
end

function b = chain3(c,n)
milky = [-128,-127,-126,-125,-124,-123,-122,-121,-120,-119,-118,-117,-116,-115,-114,-113,144;-129,98,97,96,95,94,93,92,91,90,89,88,87,86,85,-112,143;-130,99,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,84,-111,142;-131,100,-73,50,49,48,47,46,45,44,43,42,41,-60,83,-110,141;-132,101,-74,51,-32,-31,-30,-29,-28,-27,-26,-25,40,-59,82,-109,140;-133,102,-75,52,-33,18,17,16,15,14,13,-24,39,-58,81,-108,139;-134,103,-76,53,-34,19,-8,-7,-6,-5,12,-23,38,-57,80,-107,138;-135,104,-77,54,-35,20,-9,2,1,-4,11,-22,37,-56,79,-106,137;-136,105,-78,55,-36,21,-10,3,0,-3,10,-21,36,-55,78,-105,136;-137,106,-79,56,-37,22,-11,4,-1,-2,9,-20,35,-54,77,-104,135;-138,107,-80,57,-38,23,-12,5,6,7,8,-19,34,-53,76,-103,134;-139,108,-81,58,-39,24,-13,-14,-15,-16,-17,-18,33,-52,75,-102,133;-140,109,-82,59,-40,25,26,27,28,29,30,31,32,-51,74,-101,132;-141,110,-83,60,-41,-42,-43,-44,-45,-46,-47,-48,-49,-50,73,-100,131;-142,111,-84,61,62,63,64,65,66,67,68,69,70,71,72,-99,130;-143,112,-85,-86,-87,-88,-89,-90,-91,-92,-93,-94,-95,-96,-97,-98,129;-144,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128;];
inner = [1,1,2,2,3,4,5,5,6,7,8,8,9,10,11,12,13,13,14,15,16,17,18,18,19,20,21,22,23,24,25,25,26,27,28,29,30,31,32,32,33,34,35,36,37,38,39,40,41,41,42,43,44,45,46,47,48,49,50,50,51,52,53,54,55,56,57,58,59,60,61,61,62,63,64,65,66,67,68,69,70,71,72,72,73,74,75,76,77,78,79,80,81,82,83,84,85;];
outer = [4,6,7,9,10,11,12,14,15,16,17,19,20,21,22,23,24,26,27,28,29,30,31,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,62,63,64,65,66,67,68,69,70,71,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108,109,110,111,112;];

wm = 17;
L = length(c);
%Smallest square of odd side length that contains chain square +2, without exceeding board square
wc = min(round(sqrt(L)/2)*2+1+2, n);
dc = (wm-wc)/2; %Gap between milky square and chain square
cc = (wc^2+1)/2; %Center of chain square
milky = milky(1+dc:wm-dc,1+dc:wm-dc)+cc; %Select centered square of milky that fits chain square

Mr = ([1:wc]-ceil(wc/2)).^2;
M = bsxfun(@plus,Mr,Mr');

dn = (n-wc)/2;
if or(mod(L,2), L==200) %No special treatment for odd L yet
    temp = zeros(wc^2,1);
    b= zeros(n);
    
    r = floor((L-1)/2);
    temp(cc-r:cc+(L-r)-1) = 1:L;
    
    %Embed the chain square into the board square
    b(1+dn:n-dn,1+dn:n-dn) = temp(milky);
else
    Mo(milky) = M;
    b= zeros(n);
    r = floor((L-1)/2);
    temp = zeros(1,wc^2);
    
    temp(cc-r:cc+(L-r)-1) = c;
    score1 = temp*Mo';
    
    big = find(max(c)==c);
    temp(:)=0;
    temp(cc-inner(L/2-2):cc+outer(L/2-2)) = c;
    Mo = Mo(cc-inner(L/2-2):cc+outer(L/2-2));
    score2 = toeplitz([c(1),c(end:-1:2)],c)*Mo';
    [score2,ind] = min(score2);
    
    temp(:)=0;
    if score1<score2
        temp(cc-r:cc+(L-r)-1) = 1:L;
        b(1+dn:n-dn,1+dn:n-dn) = temp(milky);
    else
        temp(cc-inner(L/2-2):cc+outer(L/2-2)) = circshift([1:L]', ind-1);
        b(1+dn:n-dn,1+dn:n-dn) = temp(milky);
    end
end

end %END OF MAIN FUNCTION




