%Copyright (c) 2015, Josh Bevan
%All rights reserved.
%This code is licensed under the BSD 3-Clause License, see License.txt
%in the folder that contains this file for the full license.

%Improved precision self-enumerating pangram creator
close all
clear all
clc

N=3; %Desired number of pangram decimal accuracy
assert(N<16,'Insufficient floating point precision available')
FloatErr=eps(1);
AtoZ= char(65:90); %All uppercase latin letters

%English names for numbers 1-25, if we expect the percent incidence of a
%letter is greater than 25 than we should expand this list
NumbersText= {'ZERO' 'ONE' 'TWO' 'THREE' 'FOUR' 'FIVE' 'SIX' 'SEVEN'...
    'EIGHT' 'NINE' 'TEN' 'ELEVEN' 'TWELVE' 'THIRTEEN' 'FOURTEEN'...
    'FIFTEEN' 'SIXTEEN' 'SEVENTEEN' 'EIGHTEEN' 'NINETEEN' 'TWENTY'...
    'TWENTYONE' 'TWENTYTWO' 'TWENTYTHREE' 'TWENTYFOUR' 'TWENTYFIVE'};
%We count letters by converting text to uppercase then doing a histogram of
%the text compared against the uppercase letters
%Later on we make use of the fact that the row vector NumbersCount(N+1,:)
%contains the letter count for the number 'N'
for i=1:length(NumbersText)
    NumbersCount(i,:)= histc(upper(NumbersText{i}),AtoZ);
end

%This seed starts off the iteration with a pretty good guess
seedDEP= ['This sentence is dedicated to Chris Patuzzo\Lee Sallows and to within two decimal places'...
    'four point zero percent of the letters in this sentence are a’s'...
    'zero point one two percent are b’s, three point nine percent are c’s'...
    'zero point eight percent are d’s, nineteen point one percent are e’s'...
    'one point five percent are f’s, zero point four percent are g’s'...
    'one point five percent are h’s, six point eight percent are i’s'...
    'zero point one two percent are j’s, zero point one two percent are k’s'...
    'one point one percent are l’s, zero point three percent are m’s'...
    'twelve point one percent are n’s, eight point one percent are o’s'...
    'seven point three percent are p’s, zero point one two percent are q’s'...
    'nine point nine percent are r’s, five point six percent are s’s'...
    'nine point nine percent are t’s, zero point seven percent are u’s'...
    'one point four percent are v’s, zero point seven percent are w’s'...
    'zero point five percent are x’s, zero point three percent are y’s'...
    'and one point six percent are z’s.'];
seed= ['This sentence is dedicated to Chris Patuzzo\Lee Sallows and to within three decimal places'...
    'three point six zero percent of the letters in this sentence are a’s'...
    'zero point one zero percent are b’s, three point five zero percent are c’s'...
    'zero point seven four percent are d’s, nineteen point three zero percent are e’s'...
    'one point eight zero percent are f’s, zero point four two percent are g’s'...
    'one point eight zero percent are h’s, seven point four two percent are i’s'...
    'zero point one zero percent are j’s, zero point one zero percent are k’s'...
    'zero point seven four percent are l’s, zero point two one percent are m’s'...
    'ten point seven one percent are n’s, nine point two three percent are o’s'...
    'five point eight three percent are p’s, zero point one zero percent are q’s'...
    'nine point seven six percent are r’s, six point one five percent are s’s'...
    'nine point six five percent are t’s, zero point nine five percent are u’s'...
    'one point seven zero percent are v’s, one point seven zero percent are w’s'...
    'one point five nine percent are x’s, zero point one zero percent are y’s'...
    'and two point six five percent are z’s.'];
SeedCount= histc(upper(seed),AtoZ);
for j=[1:length(AtoZ)]+4%[1:length(AtoZ)*(N-1)-4]+5
    SeedCount= SeedCount+ NumbersCount(mod(j-1,10)+1,:);
end
%We don't explicitly store the pangram in letters/words. The fixed language
%is unchanging while the word percent values can be found from the numeric
%values. Therefore we can simply store the percent values that *generate*
%what *would* be the paragraph and operate on them.
Percent= 100*SeedCount/sum(SeedCount); Actual=Percent;

amble= ['This sentence is dedicated to Chris Patuzzo\Lee Sallows and to within' 'decimal place' 'of the letters in this sentence' 'and'];
repeated= 'point percent are s';
nLetterF= histc(upper(amble),AtoZ) + 26*histc(upper(repeated),AtoZ);
nLetterF= nLetterF+1; %a-z are each used once at least in the per letter desc
if N>1; nLetterF(19)= nLetterF(19)+1; end %For pluralizing "decimal place(s)"
nLetterF= nLetterF + NumbersCount(N+1,:); %Add N from "within N decimal place(s)"

LoopL=100;
OuterEnd=round(400000/LoopL)*8000; %How long should we search? For my machine 400-> 7 seconds
Best=0; BestA=20; nLetter=SeedCount; h=waitbar(0); tic
for Mutate=0:OuterEnd-1
    for tries=1:LoopL %We should iterate just long enough until the convergence cycle repeats
        nLetterOLD=nLetter;
        nLetter=nLetterF... %Fixed letter count plus variable counts from numeric percents
            + sum(NumbersCount( floor(Percent)+1 , : ))...
            + sum(NumbersCount( floor(10*(Percent-floor(Percent)))+1 , : ))...
            + sum(NumbersCount( floor(10*(10*Percent-floor(10*Percent)))+1 , : ))... %Always make sure the last decimal is a round() and not a floor()
            + sum(NumbersCount( round(10*(100*Percent-floor(100*Percent)))+1 , : ));
        Actual= 100*nLetter/sum(nLetter); %This is the *actual* percent incidence of the language that *would* be generated from the percent values
        
%         D(tries+Mutate*LoopL)=sum(abs(round(100*Actual)-round(100*Percent))<FloatErr);
%         if D(tries+Mutate*LoopL)>Best
%             Best= D(tries+Mutate*LoopL);
%             BestPA(1,:)= Percent;
%             BestPA(2,:)= Actual;
%         end
        %DA(tries+Mutate*LoopL)= sum(abs(nLetter-nLetterOLD));        
        if sum(abs(nLetter-nLetterOLD))<1
            BestA(1,end+1)= sum(abs(nLetter-nLetterOLD));
            BestA(2,end)= tries+Mutate*LoopL;
            BestPAa(1,:,end+1)= nLetterOLD;
            BestPAa(2,:,end)= nLetter;
            BestPAa(3,:,end)= Percent;
            BestPAa(4,:,end)= Actual;
        elseif sum(abs(nLetter-nLetterOLD))<BestA(1,end)
            BestA(1,end+1)= sum(abs(nLetter-nLetterOLD));
            BestA(2,end)= tries+Mutate*LoopL;
            BestPAa(1,:)= nLetterOLD;
            BestPAa(2,:)= nLetter;
            BestPAa(3,:)= Percent;
            BestPAa(4,:)= Actual;
        end
        %Dn(tries+Mutate*100*LoopL)= sum(abs(histc(num2str(Percent-(Percent>=10).*floor(Percent),'%.3f'),'0123456789')-sum(Percent>=10)*([1:10]<2) - histc(num2str(Actual-(Actual>=10).*floor(Actual),'%.3f'),'0123456789')-sum(Actual>=10)*([1:10]<2)));
        %nL(tries+Mutate*100*LoopL)=sum(nLetter);

        Percent=Actual;
    end
    %"Mutate" by randomly selecting a digit's word and replacing it. This
    %word may not exist, but this is unlikely and at worst wastes 1 mutate
    %cycle
    nLetter= nLetter - NumbersCount(floor(10*rand)+1,:) + NumbersCount(floor(10*rand)+1,:);
    Percent= 100*nLetter/sum(nLetter);
    waitbar(Mutate/OuterEnd,h,sprintf('%i',round((OuterEnd-Mutate)*toc/Mutate)))
end
close(h); toc

%Compare our closest success
% Percent=BestPA(1,:);
% nLetter=nLetterF...
%     + sum(NumbersCount( floor(Percent)+1 , : ))...
%     + sum(NumbersCount( floor(10*(Percent-floor(Percent)))+1 , : ))...
%     + sum(NumbersCount( floor(10*(10*Percent-floor(10*Percent)))+1 , : ))... %Always make sure the last decimal is a round() and not a floor()
%     + sum(NumbersCount( round(10*(100*Percent-floor(100*Percent)))+1 , : ));
% BestPA/100*sum(nLetter)