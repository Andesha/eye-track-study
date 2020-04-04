addpath('/triton/becs/scratch/braindata/shared/BML_utilities/eISCtools/');
%load /scratch/braindata/jlahnako/Tikku-paper/frameMapping.mat
%load /scratch/braindata/jlahnako/delete_after_reboot_tmp/eye/data_eye2.mat
%% Read the video file into VideoReader object
% Note: This takes quite a while, so don't run this every time if you don't
% need to
v=VideoReader('/scratch/braindata/bacham1/movie_stimuli/msk_movie_final.avi');

%%

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

filenames={'/scratch/braindata/bacham1/eyetracking/converted/MORAL38.asc'
           '/scratch/braindata/bacham1/eyetracking/converted/MORAL39.asc'
           '/scratch/braindata/bacham1/eyetracking/converted/MORAL40.asc'};

lognames={'/archive/braindata/2014/bacham1/logfiles/Moral_T0802-moral_movie1adopted_annafirst.log'
          '/archive/braindata/2014/bacham1/logfiles/Moral_T0802-moral_movie1adopted_annafirst.log'
          '/archive/braindata/2014/bacham1/logfiles/Moral_T0802-moral_movie1adopted_annafirst.log'};

for k=1:length(filenames)
    
    %This reads the start of video event assumed to be the start of the
    %experiment
    offsets(k)=read_logfile(lognames{k});
    
    %This reads the eye tracking data
    eyeData=read_eyedata_eyelink(filenames{k});
    
    starts(k)=eyeData.inputs(find(diff(eyeData.inputs(:,1))>4000,1)+1,1);
    tempData=eyeData.Fixations;
    tempData(:,1)=tempData(:,1)-starts(k);
    
    %Throw away the first samples and remove the offset from the remaining
    %samples
    tempData(find(tempData(:,1)<offsets(k)*eyeVsLogTime),:)=[];
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
calibrationWidth=1920;%eyeData.CalibrationAreaSize(1);%1024
calibrationHeight=1200;%eyeData.CalibrationAreaSize(2);%768

screenWidth=720;
screenHeight=576;

screenLims=[floor((calibrationWidth-screenWidth)/2) ...
    ceil((calibrationWidth+screenWidth)/2);
    floor((calibrationHeight-screenHeight)/2) ...
    ceil((calibrationHeight+screenHeight)/2)];
%%
for sub=1:length(filenames)
    subData{sub}(:,4)=subData{sub}(:,4)-screenLims(1,1);
    subData{sub}(:,5)=subData{sub}(:,5)-screenLims(2,1);
    subData{sub}((subData{sub}(:,4)<0),4)=nan;
    subData{sub}((subData{sub}(:,5)<0),5)=nan;
    subData{sub}((subData{sub}(:,4)>screenWidth),4)=nan;
    subData{sub}((subData{sub}(:,5)>screenHeight),5)=nan;
end;
    
%%
tt=25*25;%10*25;%(3*60+34)*25;
frameDist=1000/25; %in milliseconds
len=20;
cols=hot(256);
vidFrames=uint8(zeros(576,720,3,len));
%eyeData=[];

%%
vidWri=VideoWriter('/scratch/braindata/jlahnako/MareikeFixVideo2.avi');
vidWri.FrameRate=25;
vidWri.open;
for t=1:len*25
    %tic
    %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%06i / %06i',t,len*25);
    fprintf('%06i / %06i, fixations: ',t,len*25);
    tOff=t+tt;
    %im=imread(sprintf('/fs/tmpdata/tikkuVideoTemp/img/frame_%07i.png',tOff));
    im=v.read(tOff);
    imR=permute(im,[2 1 3]);%imresize(permute(im,[2 1 3]),[1024,768]);
    %imOverlay=imR;
    for sub=1:length(filenames)%[1 3 5]
        
        fixes=find((subData{sub}(:,1)<(tOff*frameDist)).*((subData{sub}(:,1)+subData{sub}(:,3))>(tOff*frameDist)));
        fprintf('sub%i, %i ',sub,length(fixes));
        if any(fixes)
            imFix=eISC_fixationHeatmap(subData{sub}(fixes,[5 4]),51/2,screenWidth,screenHeight)';
            if any(imFix(:))
                %fixes=find((subData{sub}(:,1)>tOff*frameDist).*(subData{sub}(:,1)<((tOff+1)*frameDist)));
                %imFix=plotFixations(a{sub}(fixes,[7 8]),51,1024,768);
                imFixInd=ceil(255*(imFix/max(imFix(:))))+1;
                imFixCol=cols(imFixInd,:);
                imR=uint8((repmat(1-imFix/max(imFix(:)),[1 1 3]).*double(imR))+...
                    255*(repmat(imFix/max(imFix(:)),[1 1 3]).*reshape(imFixCol,screenWidth,screenHeight,[])));
            end;
        end;
    end;
    fprintf('\n');
    vidFrame=imresize(permute(imR,[2 1 3]),[576,720]);
    writeVideo(vidWri,vidFrame);
    %toc
end
vidWri.close;
fprintf('...done\n');

%%


