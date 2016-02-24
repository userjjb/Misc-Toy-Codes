function Icrop= AutoCrop(file, thresh)
%Automatically crop pictures of documents, assumes docs are roughly aligned with frame and
%sufficient contrast is present to detect edges
%-"thresh" is a threshold value for detection, 0.3 works well from testing
I= imread(file);
Ig= double(rgb2gray(I));
[H,W]=size(Ig);

%Apply Scharr operator for edge detection
%-We are convolving each pixel with the Scharr operator, for x-direction as an example:
% |3    0   -3 |
% |10   0   -10| * A
% |3    0   -3 |
%where A is our pixel. This does two things: finds the x-gradient and smooths along the y-direction.
%-We could use larger size kerenels, but the 3x3 is the canonical example.
%-The more common operator is the Sobel-Feldman one, but Scharr has better rotational invariance.

%Gx1 and Gx2 are differential and smoothing components of the Scharr operator respectively. It is
%convenient and more efficient to construct the Scharr operator to convolve the entire image in one
%shot
Gx1=toeplitz([1 0 -1 zeros(1,W-3)],[1 zeros(1,W-5)]);
Gx2=toeplitz([3 zeros(1,H-5)],[3 10 3 zeros(1,H-3)]);
Dx=Gx2*Ig*Gx1;
Gy1=toeplitz([3 zeros(1,W-5)],[3 10 3 zeros(1,W-3)]);
Gy2=toeplitz([1 0 -1 zeros(1,H-3)],[1 zeros(1,H-5)]);
Dy=Gy2'*Ig*Gy1';

DxM=sum(Dx)/max(abs(sum(Dx))); %Normalize the gradient
DyM=sum(Dy,2)/max(abs(sum(Dy,2)));
[~,Px] = findpeaks(abs(DxM),'MinPeakHeight',thresh); %This is perhaps lazy, but it works
[~,Py] = findpeaks(abs(DyM),'MinPeakHeight',thresh);
%-We have to be careful about which peak we choose, in particular if the document is small compared
%to the total image and the background has high contrast lines we get a false positive from being
%premature, but choosing the first and last edges is easy and typically works.
%-A smarter implementation would find where the candidate left and right edges terminate and verify
%this is approximately where the top and bottom edges terminate.
Ax= Px(1);
Bx= Px(end);
Ay= Py(1);
By= Py(end);

Icrop= I(Ay:By,Ax:Bx,:);
end