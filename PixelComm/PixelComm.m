%Pixel Comm utility for Pixel Aether
%s_out= PixelComm(type, s_in)
%type can be either 'a2t' or 't2a' for ascii to tile translations
%s_in should be either the filename of a 256 color bitmap of a screencap of the tile message OR
%the text you would like to translate to tiles

function s_out= PixelComm(type, s_in)
    switch type
        case 'a2t'
            A2T(s_in)
            s_out='';
        case 't2a'
            s_out= T2A(s_in);
    end
end

function A2T(s_in)
    S=dec2base(s_in,4)';
    for it=1:size(S,1)
        for jt=1:size(S,2)
            SM(it,jt)= str2num(S(it,jt));
        end
    end
    colormap([.75 .75 .75; 1 .97 .42; 0.10 .32 1; 1 0 .18])
    imagesc(SM); axis equal; axis([0.5 length(s_in)+.5 0.5 4.5])
end

function s_out= T2A(s_in)
    [a,map]=imread(s_in);
    [yf,xf]=find(not(a==58|a==59|a==42),1,'first');
    [yl,xl]=find(not(a==58|a==59|a==42),1,'last');
    aa=a(yf:yl,xf:xl);
    aaa=reshape(permute(reshape(aa,size(aa,1),28,[]),[2 1 3]),28,35,[]);
    type=reshape(uint8(mode(mode(double(aaa)))),4,[]);
    type(type==191)=1; type(type==7)=0; type(type==217)=2; type(type==71)=3;
    s_out= char((4.^[3:-1:0])*double(type));
end