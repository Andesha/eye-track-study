function out = eISC_gaussKernel(sigma,kernRadius,viewDistance,viewWidth,viewResolution,viewScaling)
% out = eISC_gaussKernel(sigma,kernRadius)
%-------------------------------------------------------------------------------
%	Creates a 2D gaussian kernel for plotting gaze heat maps.
% All inputs are optional, but it is advicable to change them to match the experimental
% parameters (view distance, width and resolution) for optimal results.
% Default parameters are based on a single experiment.
%
% Inputs:
% sigma:			Standard deviation of the fixation kernel. By default
%                   this will be used as deviation in number of PIXELS.
%                   If view all parameters are defined this is given in
%                   visual angle (in degrees) and calculated accordingly.
% kernRadius:		Radius of the kernel matrix (i.e. output will be a
%                   square matrix with 2*kerRadius+1 columns and rows)
% viewDistance: 	Distance from the subjects' eyes to the screen in
%                   arbitrary units.
% viewWidth:		Width of the image. Should in the same units as
%                   viewDistance.
% viewResolution:	Horizontal resolution of the image
% viewScaling:		Height of a pixel on the screen divided by the width
%                   of a pixel. This is NOT the same as vertical/horizontal
%                   resolution, but rather a measure of the shape of
%                   individual image pixels in case some sort of
%                   anisotropic scaling has been applied. The results
%                   should not change massively if this is omitted and
%                   isotropic voxels are assumed if the	images have not
%                   been stretched excessively.
%
% Ouput:
%	out:			The created 2D gaussian kernel matrix
%
% Version 0.01
% 10.4.2012 Juha Lahnakoski
% juha.lahnakoski@aalto.fi

if nargin==0 || isempty(sigma)
	%These are the default settings corresponding to sigma ~1 degree
	%at 34 cm viewing distance and pixel size ~1024pix/28cm
	sigma=34*tan(pi/180)*1024/28;
end;
%Here is the full definition of the sigma, if the parameters are given
if nargin>=5 && ~isempty(viewDistance) && ~isempty(viewWidth)...
             && ~isempty(viewResolution)
         
	sigma=viewDistance*tan(sigma*pi/180)*viewResolution/viewWidth;
    
end;
%Default radius of the kernel is 3 standard deviations
if nargin<2 || isempty(kernRadius)
	kernRadius=ceil(3*sigma);
end;
%Set the view scaling
if nargin<6 ||isempty(viewScaling)
	viewScaling=1;
end;

[X,Y]=ndgrid(viewScaling*-kernRadius:viewScaling:kernRadius*viewScaling,...
                         -kernRadius:kernRadius);
out = exp(-(X.^2+Y.^2)/(2*sigma^2)) / sqrt(2*pi*sigma^2);

%Normalize the kernel
out=out/sum(out(:));
end
