function [r,cMat] = eISC_spatialSimilarity(fixMaps)
% [r,cMat] = eISC_spatialSimilarity(fixMaps)
% ------------------------------------------------------------------------
% Calculates the intersubject spatial correlations of fixation heatmaps.
% Alternative similarity or distance measures may be added later.
%
% Inputs:
% fixMaps:  width*height*subjects or height*width*subjects array of
%           fixation heatmaps, i.e. fixMaps(:,:,i) is the fixation maps of
%           the i'th subject.
%
% Output:
% r:        Mean pairwise spatial correlation coefficient across subjects
% cMat:     Upper triangle entries of the pairwise correlation matrix
%
% Version 0.01
% 10.4.2012 Juha Lahnakoski
% juha.lahnakoski@aalto.fi

%Calculate the correlations and select the upperm triangle entries without
%the diagonal (i.e. triu(...,1))
passThrough = reshape(fixMaps,[],size(fixMaps,3));
cMat=corr(passThrough);
cMat=cMat(find(triu(ones(size(cMat)),1)));

%Calculate the mean using the Fisher Z-transform first (atanh) and then
%transforming back.
r=tanh(nanmean(atanh(cMat)));

end
