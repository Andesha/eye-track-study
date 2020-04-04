function out = eISC_fixationHeatmap(points,kern,w,h,duration)
% out = eISC_fixationHeatmap(points,kern,w,h)
% ------------------------------------------------------------------------
% Creates heatmaps of gaze locations.
%
% Kernel dimensions should be odd so that the maximum of the kernel is
% in the middle pixel. This function will return an error if size in
% either direction is even. Kernels should also be square. If you do not
% want your kernel to be isotropically distributed you may still use a
% square matrix.
%
% Inputs:
% points:   Fixation locations in a fixations*coordinates matrix,
%           i.e. each row corresponds to x and y coordinates of fixation.
% kern:     External kernel or standard deviation of kernel in PIXELS.
% w:        Width of the image in pixels.
% h:        Height of the image in pixels.
%in
% Output:
% out:      Heatmap of fixations as the normalized sum of fixation kernels
%           at each fixation location.
%
% Version 0.01
% 10.4.2012 Juha Lahnakoski
% juha.lahnakoski@aalto.fi

%Defaults for the image size from a single experiment
if nargin<4 || isempty(h)
    h=1200;
end;
if nargin<3 || isempty(w)
    w=1920;
end;

if nargin>=2 && ~isempty(kern) && min(size(kern))==1
    %If kern-parameter is a plain number we use it as the radius and create
    %a new kernel
    kern=eISC_gaussKernel(kern);
    kernRadius=(size(kern,1)-1)/2;
elseif nargin>=2 && ~isempty(kern) && min(size(kern))>1
    %If kern-parameter is a matrix let's calculate the radius as the 
    kernRadius=(max(size(kern))-1)/2;
    
    %We do not want to use kernels where the size in any direction is even
    if min(mod(size(kern),2))==0%max(iseven(size(kern)))
        error('eISC_fixationHeatmap supports only kernels where the size along each dimension is odd.');
    end;
    
    %Also, we only want square kernels for now
    if min(size(kern))~=max(size(kern))
        error('eISC_fixationHeatmap supports only square kernels.');
    end;
        
else
    %If kern-parameter is not defined let's make the default kernel
    kern=eISC_gaussKernel;
    kernRadius=(size(kern,1)-1)/2;
end;

if size(kern,1)~=size(kern,2)
    %Give an error if the kernel is not square
    error('eISC currently allows only square kernels.');
end;

%Normalize the kernel just in case it wasn't already
kern=kern-min(kern(:));
kern=kern/max(kern(:));

%Create empty output
out=zeros(h,w);

%Find the non-NaN entries
idx=find(~isnan(points(:,1)).*~isnan(points(:,2)));
x=points(idx,2);
y=points(idx,1);
if nargin<5 || isempty(duration)
    duration=ones(length(x),1);
end

%Loop through fixations taking into account that some fixations will go
%over the image edges
for jj=1:length(x)
    %This calculates where the kernel ends
    xLo=round(x(jj))-kernRadius;
    xHi=round(x(jj))+kernRadius;
    yLo=round(y(jj))-kernRadius;
    yHi=round(y(jj))+kernRadius;
    %This calculates where we should cut the kernel if it goes outside the
    %image
    xKernLo=max(1,-xLo+2);
    xKernHi=min(2*kernRadius+1,2*kernRadius+1-(xHi-w));
    yKernLo=max(1,-yLo+2);
    yKernHi=min(2*kernRadius+1,2*kernRadius+1-(yHi-h));
    xLo=max(1,xLo);
    xHi=min(w,xHi);
    yLo=max(1,yLo);
    yHi=min(h,yHi);
    %And here the kernel is added to the image.
    out(yLo:yHi,xLo:xHi)=out(yLo:yHi,xLo:xHi)+duration(jj)*kern(yKernLo:yKernHi,xKernLo:xKernHi);
end;
%Normalization of the image
if sum(out(:))~=0
out=out/sum(out(:));
end;

end
