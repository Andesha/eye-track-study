%% eISC_runAnalysis.m: Script for calculating fixation similarities over
%  time windows. Requires other files in the eISCtools-folder.
%-------------------------------------------------------------------------
% This file contains some basic examples of how to use the few functions
% that have currently been written. If there are any questions or
% suggestions for additions, please email me.
%
% Results will be saved to eISC_results structure as follows. Please note
% that the fixation heatmaps are not saved because of excessive memory load
% in long experiments.
%
% eISC_results.eISC:            Mean similarity over all subjects (in time)
% eISC_results.eISCmat:         Upper triangle entries of similarity matrix
%                               over all subjects (pair X time).
%
% eISC_results.groupComp.groupMeans:  Mean similarities within groups.
%                                     (groups X timepoints)
% eISC_results.groupComp.betweenGroups: Mean similarities between groups.
%                                       One value per time point for two
%                                       groups, otherwise size is 
%                                       [groups X groups X timepoints]
% eISC_results.groupComp.p:     Uncorrected p-values of group comparison
%                               One value per time point for two
%                               groups, otherwise size is 
%                               [groups X groups X timepoints]
% eISC_results.groupComp.h:     Results of hypothesis testing at the
%                               specified p-value threshold and multiple
%                               comparisons correction (0=no difference,
%                               1=difference). Size as above.
% More to be added...
%
% Version 0.01
% 11.4.2012 Juha Lahnakoski
% juha.lahnakoski@aalto.fi
%-------------------------------------------------------------------------
%% Set various parameters here
%-------------------------------------------------------------------------
clear

addpath('/triton/becs/scratch/braindata/shared/BML_utilities/eISCtools/');

%This makes the log file offset correspond to time units in the eye gaze
%data (read_logfile.m returns offset in seconds, eyedata is in
%milliseconds).
eyeVsLogTime=1000;

%Time window length in milliseconds. Can be a single number to make all
%windows the same length or a vector with the size corresponding to
%winOnsets
winLength=1000;

%Here is the vector of window starting points in milliseconds from the
%beginning of the experiment. If this is empty the first window starts at
%zero and each successive window will jump one window length forward
winOnsets=[];
%winOnsets=1:40;

%Offsets of the beginning of the data of each subject in milliseconds.
%Note that the data during this offset period will be discarded and it will
%affect the timing of the analysis windows. You should set this to reflect
%the start of the experiment. Make this zero (or equal for all) if the data
%gathering started at the same time for each subject.
%Note: this is overriden in the loop reading the files, so if you want to
%set these manually you should edit the next section
offsets=[];

%Heatmap parameters for fixation plotting. If these are empty (i.e. set to
%[]) the defaults will be used. Note that kernelSigma is in pixels if the
%viewing distance, width and resolution are not defined, otherwise it is
%degrees of the visual field. The kernelRadius-parameter refers to the actual
%size of the matrix containing the kernel (sides are 2*kernelRadius+1).
%Default value is usually fine.
%View scaling refers to the scaling of the image (if for some reason the
%scaling is different for x and y-axes.
%For more info, see eISC_gaussKernel.m

kernelSigma=[];
kernelRadius=[];
viewDistance=[];
viewWidth=[];
viewResolution=[];
viewScaling=[];

%% Data loading
%Here is the actual data (subData cell-array with each cell corresponding
%to one subject and first column of the contained matrix corresponding to
%the time code, and second and third column to the x and y-coordinates).
%Now this reads the inputs from file format of the old gaze camera at AMI
%centre, and experiment starting time from log file of Presentation
%software, but if your data is in a different format you can just load it
%into the array manually.
%There are many ways to read this data depending on the file format etc.
%Eventually we would like to add functions for reading the data here once
%we know what types of files we need to read.

filenames={
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1002/Moral_T1002_1/movie1.asc'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1003/Moral_T1003_2/movie1.asc'
    };


lognames={
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1002_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1003_2_movie1.log'
};

for k=1:length(filenames)
    
    %This reads the start of video event assumed to be the start of the
    %experiment
    offsets(k)=read_logfile_MORAL(lognames{k});
    
    %This reads the eye tracking data
    %load(sprintf('%s',filenames{k})) %Choose this when you have saved the splitted files in mat format
    eyeData=read_eyedata(filenames{k});
    
    %These are the default locations of time code, horizontal and vertical
    %positions in the data field, respectively.
    tempData=eyeData.Fixations(:,[1 4 5]);
    
    %Throw away the first samples and remove the offset from the remaining
    %samples
    %tempData(find(tempData(:,1)<offsets(k)*eyeVsLogTime),:)=[];
    tempData(:,1)=tempData(:,1)-offsets(k)*eyeVsLogTime;
    
    %Save the data to cell array because there may be a different number of
    %data points for each subject
    subData{k}=tempData;
end;

%%------------------------------------------------------------------------
%% Screen size parameters

%Screen size information. Here we assume that the screen is in the middle
%of the calibration window but does not fill the whole window. Any samples
%outside of the screen are omitted. If you want to use all the samples you
%may set the screen with and height equal to the calibration window size.
calibrationWidth=eyeData.CalibrationAreaSize(1);%1024
calibrationHeight=eyeData.CalibrationAreaSize(2);%768

screenWidth=760;
screenHeight=578;

screenLims=[floor((calibrationWidth-screenWidth)/2) ...
    ceil((calibrationWidth+screenWidth)/2);
    floor((calibrationHeight-screenHeight)/2) ...
    ceil((calibrationHeight+screenHeight)/2)];

%%------------------------------------------------------------------------
%% Group comparison setup
% These set the properties of the group comparison. There would be various
% ways to do comparisons between groups. The only one implemented for now
% is comparing the similarity values within and between groups.

%This sets if the group comparison is performed (1) or not (0)
groupComparison=1;

%This sets the group IDs for each subject. This should have as many entries
%as there are subjects, otherwise we get an error. Please use simple
%numeric IDs here and keep track which group is which.
%groups=[1 1 1 1 2 2 2 2];
groups=[1 2];

%Paired test (1) or not (2). If you have the same subjects in the same
%order you probably want a paired test. Otherwise leave this to zero.
pairedSamples=0;

%p-value threshold
pThr=0.05;

%Correction approach. This corrects for the number of time points with the
%following options:
%0: No correction
%1: FDR correction, pID
%2: FDR correction, pN
%3: Bonferroni corrrection
mcCor=1;





% -------------------------------------------------------------------------
%% Run the analysis
% There should be no need to change things below this unless new analyses
% are added. Setting the parameters above should be enough.
% -------------------------------------------------------------------------

%Save the number of subjects to make the code a bit shorter
nSubs=length(subData);
tLength=zeros(nSubs,1);
for sub=1:nSubs
    tLength(sub)=max(subData{sub}(:,1));
end;
%Stop the comparisons when the data for the shortest dataset ends
tEnd=min(tLength);
%-------------------------------------------------------------------------
%This is for equal length windows
if isempty(winOnsets)
    
    %Pre-allocate the memory for results
    eISC_results.eISC=zeros(round(tEnd/winLength),1);
    eISC_results.eISCmat=zeros(round(tEnd/winLength),nSubs*(nSubs-1)/2);
    
    fixMaps=zeros([diff(screenLims') nSubs]);
    
    fprintf('Processing time window: ');
    for t=1:tEnd/winLength
        fprintf('%i,',t);
        if mod(t,100)==0
            fprintf('\n');
        end;
        for sub=1:nSubs
            
            tPoints=find((subData{sub}(:,1)>(t-1)*winLength).*(subData{sub}(:,1)<(t*winLength)));
            fixMaps(:,:,sub)=eISC_fixationHeatmap(subData{sub}(tPoints,2:3)-repmat(screenLims(:,1),[1 length(tPoints)])',...
                eISC_gaussKernel(kernelSigma,...
                kernelRadius,...
                viewDistance,...
                viewWidth,...
                viewResolution,...
                viewScaling),...
                screenWidth,screenHeight);
            
        end;
        [eISC_results.eISC(t),eISC_results.eISCmat(t,:)]=...
            eISC_spatialSimilarity(fixMaps);
        
        if groupComparison==1
            if length(unique(groups))<=2
            [eISC_results.groupComp.p(t),...
             eISC_results.groupComp.groupMeans(:,t),...
             eISC_results.groupComp.betweenGroups(:,:,t)] =...
                    eISC_groupComparison(eISC_results.eISCmat(t,:),groups,pairedSamples);
            else
            [eISC_results.groupComp.p(:,:,t),...
             eISC_results.groupComp.groupMeans(:,t),...
             eISC_results.groupComp.betweenGroups(:,:,t)] =...
                    eISC_groupComparison(eISC_results.eISCmat(t,:),groups,pairedSamples);
            end;
            
        end;
    end;
fprintf('\n');

%-------------------------------------------------------------------------
%Next for the case where we have preset the window onsets and lengths
else
    if length(winLength)>1 && length(winOnsets) ~= length(winLength)
        error('Window onset and length vectors should be of equal length');
    end;
    if length(winLength)==1
        winLength=ones(size(winOnsets))*winLength;
    end;
    %Pre-allocate the memory for results
    eISC_results.eISC=zeros(length(winOnsets),1);
    eISC_results.eISCmat=zeros(length(winOnsets),nSubs*(nSubs-1)/2);
    
    if groupComparison==1
        if length(unique(groups))<=2
            eISC_results.groupComp.p=zeros(length(winOnsets),1);
            eISC_results.groupComp.groupMeans=...
                zeros(length(unique(groups)),length(winOnsets));
            eISC_results.groupComp.betweenGroups=zeros(length(winOnsets),1);
        else
            eISC_results.groupComp.p=...
                zeros(length(unique(groups)),length(unique(groups)),...
                length(winOnsets));
            eISC_results.groupComp.groupMeans=...
                zeros(length(unique(groups)),length(winOnsets));
            eISC_results.groupComp.betweenGroups=...
                zeros(length(unique(groups)),length(unique(groups)),...
                length(winOnsets));
        end;
    end;

    
    fixMaps=zeros([diff(screenLims') nSubs]);
    
    fprintf('Processing time window: ');
    for t=1:length(winOnsets)
        fprintf('%i,',t);
        if mod(t,100)==0
            fprintf('\n');
        end;
        for sub=1:nSubs
            
            tPoints=find((subData{sub}(:,1)>winOnsets(t)).*(subData{sub}(:,1)<winOnsets(t)+winLength(t)));
            fixMaps(:,:,sub)=eISC_fixationHeatmap(subData{sub}(tPoints,2:3)-repmat(screenLims(:,1),[1 length(tPoints)])',...
                eISC_gaussKernel(kernelSigma,...
                kernelRadius,...
                viewDistance,...
                viewWidth,...
                viewResolution,...
                viewScaling),...
                screenWidth,screenHeight);
        end;
        
        figure;for k=1:4 subplot(2,2,k);imagesc(fixMaps(:,:,k)); end
        
        
        [eISC_results.eISC(t),eISC_results.eISCmat(t,:)]=...
            eISC_spatialSimilarity(fixMaps);
        if groupComparison==1
            if length(unique(groups))<=2
            [eISC_results.groupComp.p(t),...
             eISC_results.groupComp.groupMeans(:,t),...
             eISC_results.groupComp.betweenGroups(t)] =...
                    eISC_groupComparison(eISC_results.eISCmat(t,:),groups,pairedSamples);
            else
            [eISC_results.groupComp.p(:,:,t),...
             eISC_results.groupComp.groupMeans(:,t),...
             eISC_results.groupComp.betweenGroups(:,:,t)] =...
                    eISC_groupComparison(eISC_results.eISCmat(t,:),groups,pairedSamples);
            end;
            
        end;
    end;
fprintf('\n');    
end;

%% Correcting for multiple comparisons

switch mcCor
    case 0
        eISC_results.groupComp.h=(eISC_results.groupComp.p<pThr)+(eISC_results.groupComp.p>(1-pThr));
        break;
    case 1
        p1thr=FDR(eISC_results.groupComp.p(find(eISC_results.groupComp.p)),pThr);
        eISC_results.groupComp.h=(eISC_results.groupComp.p<p1thr)+(eISC_results.groupComp.p>(1-p1thr));
        break;
    case 2
        [~,p2thr]=FDR(eISC_results.groupComp.p(find(eISC_results.groupComp.p)),pThr);
        eISC_results.groupComp.h=(eISC_results.groupComp.p<p2thr)+(eISC_results.groupComp.p>(1-p2thr));
        break;
    case 3
        pLength=size(eISC_results.groupComp.p);
        pLength=pLength(end);
        eISC_results.groupComp.h=(eISC_results.groupComp.p<(p/pLength))+(eISC_results.groupComp.p>(1-(p/pLength)));
        break;
        
end;
