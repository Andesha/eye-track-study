%Time window length in milliseconds. Can be a single number to make all
%windows the same length or a vector with the size corresponding to
%winOnsets
winLength=500;

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
folder = input('Enter Folder Containing Data: ', 's');
file = getfiles(folder);
filenames = {file(:).n};

% Import data
for k=1:length(filenames)
    %This reads the eye tracking data
    %eyeData = read_data(); %% NOT FINISHED/READY. 
    %eyeData=read_eyedata(filenames{k}); % prolly have to make our own
    disp(['Loading: ' filenames{k}]);
    eyeData = importdata([folder '/' filenames{k}]);
    
    % Text transformation of timestamps because matlab is DUMB >:(
    x = eyeData.textdata(:,3);
    x = x(2:end);
    y = str2double(x);
    colIndex = min(size(eyeData.data) + 1);
    eyeData.data(:,colIndex) = y;
    
    %These are the default locations of time code, horizontal and vertical
    %positions in the data field, respectively.
    tempData=eyeData.data(:,[colIndex 2 3]); % 1 3 4 originally
        
    %Throw away the first samples and remove the offset from the remaining
    %samples
    % FILL IN OUR OWN METHOD
    tempData = tempData(50:end,:);
    
    % Timetamp correction:
    minStamp = tempData(1,1);
    tempData(:,1) = tempData(:,1) - minStamp;
    tempData(:,1) = tempData(:,1) / 1000;
        
    %Save the data to cell array because there may be a different number of
    %data points for each subject
    subData{k}=tempData;
end

disp('Done loading files...');

clear tempData;
clear eyeData;

%%------------------------------------------------------------------------
%% Screen size parameters
%Screen size information. Here we assume that the screen is in the middle
%of the calibration window but does not fill the whole window. Any samples
%outside of the screen are omitted. If you want to use all the samples you
%may set the screen with and height equal to the calibration window size.

%% EMILY INFO Jan9. 
% display pixels: 1440 x 900
% display size(lxh): 34.29 x 25.91 (cm). 
% max screen pixels: 1920 x 1200
% distance from screen: 63 cm. 
% visual angle: 23 degrees.  
%% FINISH 

calibrationWidth = 1920; %eyeData.CalibrationAreaSize(1);%1024
calibrationHeight = 1200; %eyeData.CalibrationAreaSize(2);%768

screenWidth = 1440; %760;
screenHeight = 900; %578;

screenLims=[floor((calibrationWidth-screenWidth)/2) ...
    ceil((calibrationWidth+screenWidth)/2);
    floor((calibrationHeight-screenHeight)/2) ...
    ceil((calibrationHeight+screenHeight)/2)];

% tEnd is based on the maximum number of blocks that can fit into the
% shortest file, i.e. you can't have group correlations between files that
% go on for different amounts of time.
nSubs=length(subData);
tLength=zeros(nSubs,1);
for sub=1:nSubs
    tLength(sub)=max(subData{sub}(:,1));
end
%Stop the comparisons when the data for the shortest dataset ends
tEnd=min(tLength);
    
%Pre-allocate the memory for results
fixMaps=zeros([diff(screenLims') nSubs]);

disp('Done parameter loading and preloading...');
for t=1:tEnd/winLength
    removed = 0;
    removedAt = [];
    for sub=1:nSubs
        % Grab the row numbers which are valid in the current block.
        % Note that the number of points valid per block is variable via
        % this code.
        tPoints=find((subData{sub}(:,1)>(t-1)*winLength).*(subData{sub}(:,1)<(t*winLength)));
        
        % Check if this sub has NaN entries in the block - if so skip him.
        if invalidPoints(subData{sub},tPoints)
            removed = removed + 1;
            removedAt = [removedAt sub];
            continue;
        end
        % Using tPoints, grab the x,y pairs and correct them via repmat
        % approach.
        % Gauss kernel is fed in as the second argument.
        singleMap = eISC_fixationHeatmap(subData{sub}(tPoints,2:3),...
            eISC_gaussKernel(kernelSigma,...
                kernelRadius,...
                viewDistance,...
                viewWidth,...
                viewResolution,...
                viewScaling),...
            screenWidth,screenHeight);
        fixMaps(:,:,sub) = singleMap';
    end
    
    % Fancy stats! ~ would contain the table of pairwise correlations if
    % you wanted it though.
    outString = tylers_similarity(fixMaps,removedAt,t);
    % [r,~] = eISC_spatialSimilarity(fixMaps);
    % Eventually this output needs to be made fancy
    % disp([num2str(t*winLength) ',' num2str(r) ',' num2str(removed)]);
    disp([num2str(t*winLength) ',' num2str(removed) ',' outString]);
end
disp('Finished!');