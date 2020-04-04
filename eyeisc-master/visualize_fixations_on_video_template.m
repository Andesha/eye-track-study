addpath('/scratch/braindata/shared/BML_utilities/eISCtools')
addpath('/scratch/braindata/shared/BML_utilities/')
%% Set various parameters here
%-------------------------------------------------------------------------
%clear

%This makes the log file offset correspond to time units in the eye gaze
%data (read_logfile.m returns offset in seconds, eyedata is in
%milliseconds).
eyeVsLogTime=1000;

%Time window length in milliseconds. Can be a single number to make all
%windows the same length or a vector with the size corresponding to
%winOnsets
winLength=2000;

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
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0101/Moral_T0101_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0101/Moral_T0101_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0102/Moral_T0102_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0102/Moral_T0102_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0103/Moral_T0103_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0202/Moral_T0202_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0202/Moral_T0202_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0302/Moral_T0302_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0302/Moral_T0302_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0303/Moral_T0303_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0303/Moral_T0303_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0401/Moral_T0401_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0401/Moral_T0401_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0402/Moral_T0402_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0402/Moral_T0402_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0501/Moral_T0501_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0502/Moral_T0502_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0502/Moral_T0502_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0503/Moral_T0503_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0601/Moral_T0601_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0601/Moral_T0601_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0601/Moral_T0601_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0602/Moral_T0602_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0602/Moral_T0602_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0603/Moral_T0603_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0603/Moral_T0603_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0701/Moral_T0701_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0701/Moral_T0701_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0702/Moral_T0702_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0702/Moral_T0702_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0702/Moral_T0702_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0703/Moral_T0703_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0703/Moral_T0703_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0703/Moral_T0703_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0801/Moral_T0801_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0801/Moral_T0801_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0801/Moral_T0801_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0801/Moral_T0801_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0802/Moral_T0802_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0802/Moral_T0802_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0803/Moral_T0803_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0803/Moral_T0803_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0803/Moral_T0803_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0803/Moral_T0803_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0901/Moral_T0901_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0901/Moral_T0901_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0901/Moral_T0901_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0902/Moral_T0902_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0902/Moral_T0902_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0903/Moral_T0903_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0903/Moral_T0903_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T0903/Moral_T0903_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1001/Moral_T1001_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1001/Moral_T1001_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1001/Moral_T1001_2/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1001/Moral_T1001_2/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1002/Moral_T1002_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1002/Moral_T1002_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1003/Moral_T1003_1/movie1.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1003/Moral_T1003_1/movie2.mat'
'/scratch/braindata/bacham1/eyetracking/accepted/Moral_T1003/Moral_T1003_2/movie1.mat'
    };


lognames={
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0101_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0101_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0102_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0102_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0103_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0202_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0202_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0302_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0302_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0303_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0303_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0401_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0401_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0402_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0402_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0501_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0502_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0502_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0503_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0601_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0601_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0601_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0602_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0602_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0603_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0603_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0701_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0701_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0702_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0702_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0702_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0703_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0703_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0703_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0801_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0801_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0801_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0801_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0802_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0802_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0803_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0803_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0803_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0803_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0901_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0901_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0901_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0902_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0902_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0903_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0903_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T0903_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1001_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1001_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1001_2_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1001_2_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1002_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1002_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1003_1_movie1.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1003_1_movie2.log'
'/scratch/braindata/bacham1/data_analysis/experimental_results/logfiles/Moral_T1003_2_movie1.log'
};

for k=1:length(filenames)
    
    %This reads the start of video event assumed to be the start of the
    %experiment
    offsets(k)=read_logfile_MORAL(lognames{k});
    
    %This reads the eye tracking data
    load(sprintf('%s',filenames{k})) %Choose this when you have saved the splitted files in mat format
    %eyeData=read_eyedata_eyelink(filenames{k});
    
    %These are the default locations of time code, horizontal and vertical
    %positions in the data field, respectively.
    
    tempData=eyeData.Fixations;%data(:,[1 7 8]);
    startDiffs(k)=eyeData.inputs(1,1)-tempData(1,1);
    tempData(:,1:2)=tempData(:,1:2)-eyeData.inputs(1,1);
    %Throw away the first samples and remove the offset from the remaining
    %samples
    tempData(find(tempData(:,1)<offsets(k)*eyeVsLogTime),:)=[];
    tempData(:,1:2)=tempData(:,1:2)-tempData(1,1);
    
    %Save the data to cell array because there may be a different number of
    %data points for each subject
    subData{k}=tempData;
    meanY(k)=mean(tempData(:,5));
    meanX(k)=mean(tempData(:,5));
end;
grandMeanLoc=mean(cat(2,meanX(:),meanY(:)),2);
% grandMeanLoc(1)=1920/2;
% % Calibrate the mean fixation locations
for k=1:length(subData)
    %subData{k}(:,4)=subData{k}(:,4)+grandMeanLoc(1)-meanX(1,k);
    subData{k}(:,5)=subData{k}(:,5)+grandMeanLoc(2)-meanY(k);
    meanLoc(:,k)=mean(subData{k}(:,4:5));
end;
%% Screen size parameters

%Screen size information. Here we assume that the screen is in the middle
%of the calibration window but does not fill the whole window. Any samples
%outside of the screen are omitted. If you want to use all the samples you
%may set the screen with and height equal to the calibration window size.
calibrationWidth=eyeData.CalibrationAreaSize(1);%1024
calibrationHeight=eyeData.CalibrationAreaSize(2);%768

screenWidth=720;%calibrationWidth;
screenHeight=576;%calibrationHeight;

screenLims=[floor((calibrationWidth-screenWidth)/2) ...
    ceil((calibrationWidth+screenWidth)/2);
    floor((calibrationHeight-screenHeight)/2) ...
    ceil((calibrationHeight+screenHeight)/2)];

%% Load ROIs and compare to fixation locations
cd /scratch/braindata/jlahnako/dilemma/anno/
vidObj=VideoReader('/scratch/braindata/ryyppov1/movie_stimuli/dilemmapart1.avi');
d=dir('side_1*.bmp');
%% The ROI each subject viewed at each time point is saved into "fixes" variable
% Vectors "anna", "kate" and "both" contain annotations of timepoints when
% the corresponding person(s) were visible
vidOut=VideoWriter('/scratch/braindata/jlahnako/dilemma/gaze_demo.avi');
vidOut.FrameRate=vidObj.FrameRate;
vidOut.open;
hpl=[];
fixes=[];
cols=jet(length(subData));
toffset=000;
figure;
temp=vidObj.read(1);
imh=image(zeros(size(temp(:,:,1))));
ht=text(size(temp,2)/2,50,sprintf('%i%% Anna, %i%% Kate',0,0));
set(ht,'color',[1 1 1],'HorizontalAlignment','center');
anna=zeros(vidObj.Duration*vidObj.FrameRate,1);
kate=zeros(vidObj.Duration*vidObj.FrameRate,1);
both=zeros(vidObj.Duration*vidObj.FrameRate,1);
for t=1:vidObj.Duration*vidObj.FrameRate
    if mod(t,vidObj.FrameRate)==0
        fprintf('%i,',t);
        if mod(t,vidObj.FrameRate*10)==0
            fprintf('\n');
        end;
    end;
    tt=t*1000/25+toffset/25;
%Video only
%     im=vidObj.read(t+toffset);
%     imagesc(im);

%Video + ROI
    try
        im1=double(imread(sprintf('/scratch/braindata/jlahnako/dilemma/anno/side_1_%06i.bmp',t+toffset)));
    catch e
        im1=double(zeros(size(temp(:,:,1))));
    end;
    try
        im2=double(imread(sprintf('/scratch/braindata/jlahnako/dilemma/anno/side_2_%06i.bmp',t+toffset)));
    catch e
        im2=double(zeros(size(temp(:,:,1))));
    end;
    im1=im1(:,:,1);
    im2=im2(:,:,1);
    im=im1+2*im2;
    imVid=vidObj.read(t+toffset);
    imVid=reshape(imVid,[],3);
    imVid(logical(im1),1)=64+imVid(logical(im1),1);
    imVid(logical(im2),2)=64+imVid(logical(im2),2);
    imVid=reshape(imVid,size(im1,1),size(im1,2),[]);
    set(imh,'cdata',imVid);
        
    hold on;
    set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1]);
    set(gcf,'position',[600 500 size(im,2) size(im,1)]);
    box off;
    for k=1:length(subData)
        tpoints=find((subData{k}(:,1)<=tt).*(subData{k}(:,2)>=tt));
        if any(tpoints)
            
            fixes(t,k)=im(min(max(1,round(subData{k}(tpoints,5)-screenLims(2,1))),screenHeight),min(max(1,round(subData{k}(tpoints,4)-screenLims(1,1))),screenWidth));
            if length(hpl)>=k && ishandle(hpl(k)) && any(hpl(k))
                set(hpl(k),'Xdata',subData{k}(tpoints,4)-screenLims(1,1),'Ydata',subData{k}(tpoints,5)-screenLims(2,1),'color',cols(k,:),'MarkerSize',max(1,round(subData{k}(tpoints(1),3)/30)));
            else
                hpl(k)=plot(subData{k}(tpoints,4)-screenLims(1,1),subData{k}(tpoints,5)-screenLims(2,1),'o--','color',cols(k,:),'MarkerSize',max(1,round(subData{k}(tpoints(1),3)/30)));
            end
            
        else
            t1=find(subData{k}(:,1)<=tt,1,'last');
            t2=find(subData{k}(:,1)>=tt,1,'first');
            if any(t1) && any(t2)
                
                if length(hpl)>=k && ishandle(hpl(k)) && any(hpl(k))
                    set(hpl(k),'Xdata',subData{k}([t1 t2],4)-screenLims(1,1),'Ydata',subData{k}([t1 t2],5)-screenLims(2,1),'color',cols(k,:)/2);
                else
                    hpl(k)=plot(subData{k}([t1 t2],4)-screenLims(1,1),subData{k}([t1 t2],5)-screenLims(2,1),'o--','color',cols(k,:)/2);
                end
            end;
            fixes(t,k)=0;
        end;
        
    end;
    
    set(ht,'String',sprintf('%i%% Anna, %i%% Kate',round(100*(sum(fixes(t,:)==1)+sum(fixes(t,:)==3))/size(fixes,2)),round(100*(sum(fixes(t,:)==2)+sum(fixes(t,:)==3))/size(fixes,2))));
    both(t)=length(unique(im(:)))>2;
    anna(t)=any(im1(:));
    kate(t)=any(im2(:));
    vidOut.writeVideo(getframe);
    hold off
    
end;
vidOut.close;
