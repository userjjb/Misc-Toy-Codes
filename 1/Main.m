clear all
Gi=[4.478:0.0001:10];
Gii=[0:0.0001:4];
f1=fliplr([29:0.1:30]);
D = zeros(length(f1),1);

for indx1=1:length(Gi);
    disp(Gi(indx1));
    for indx2=1:length(Gii); 
        for indx3=1:length(f1);
            pf=imag(Jackle29(Gi(indx1),Gii(indx2) ,9.81 ,1005.9 ,0.007 ,555.2376239373641 ,0.000652006 ,6.2024 ,f1(indx3)));
            D(indx3)=pf;
        end
        if issorted(D);continue
        end
        [pks,locs] = findpeaks(D);
        pks=max(pks);

        if isempty(pks);continue
        end
        if f1(locs) ~= 29.4;continue
        end
        if pks > 2.118e-05;continue
        end
        if pks < 2.113e-05;continue
        end
        plot(f1, abs(D)); xlabel('frequency (Hz)'); ylabel('chi''');
        hold on 
    end
end
hold off