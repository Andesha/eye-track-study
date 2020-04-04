function eyedata=read_eyedata_eyelink(eyefile)
%eyedata=read_eyedata(eyefile)
%  Read eye tracking data
%   eyefile : path to the eye track file
%   eyedata : struct with following fields with example values
%       CalibrationAreaSize: [1024 768]
%                SampleRate: 1000
%                 Threshold: Info for threshold (2 numbers, what ever they
%                            are)
%                   RawData: ...
%                 Fixations: Start time, end time, length, coordinates,
%                            pupil size
%                  Saccades: Start time, end time, duration, start coord,
%                            end coord, saccade amplitude??, some other
%                            number
%                    Blinks: Start time, end time, length
% V 0.00001 10.12.2013 Juha Lahnakoski, juha.lahnakoski@aalto.fi
% Partially based on code from Jouko Lampinen's read_eyedata function

eyedata=[];

f=fopen(eyefile,'r');
if f==-1,
    error(['Cannot open the file ' eyefile]);
end

% % find the eof and read the last line to estimate the number of lines
% % in the file, to preallocate the data matrix
% fseek(f,-2,'eof');
% while(1)
%     byte=fread(f,1,'uchar');
%     if byte==10,
%         break
%     end
%     fseek(f,-2,'cof');
% end
% lastLine=readln(f);
%
% fseek(f,0,'bof');

dataNo=0;
fixNo=0;
sacNo=0;
blnkNo=0;
inputNo=0;
% read the header
while(1)
    
    dataLine=fgetl(f);
    % some Windows versions return empty line from \n\r pair
    if isempty(dataLine),
        dataLine=fgetl(f);
    end
    
    if dataLine==-1
        break;
        
        %If the dataLine does not start with a number (time stamp) read in the
        %header info or fixation/saccade/blink event
    elseif ~isempty(dataLine) && isnan(str2double(dataLine(1))) %&& ~isnumeric(str2double(dataLine(1)))
        
        spcInd=strfind(dataLine,' ');
        tabInd=[find(dataLine==9) length(dataLine)];
        
        if length(spcInd)>1
            header1=setstr(dataLine(1:spcInd(1)-1));
            header2=setstr(dataLine(spcInd(1)+1:spcInd(2)-1));
            switch header1
                case 'EFIX',
                    fixNo=fixNo+1;
                    %                     fprintf('Fixation #%i\n',fixNo);
                    eyedata.Fixations(fixNo,:)=str2num(dataLine(spcInd(2):end));
                case 'ESACC',
                    if ~isempty(str2num(dataLine(spcInd(2):end)))
                        sacNo=sacNo+1;
                        %                     fprintf('Saccade #%i\n',sacNo);
                        eyedata.Saccades(sacNo,:)=str2num(dataLine(spcInd(2):end));
                    end;
                case 'EBLINK'
                    blnkNo=blnkNo+1;
                    %                     fprintf('Blink #%i\n',blnkNo);
                    eyedata.Blinks(blnkNo,:)=str2num(dataLine(spcInd(2):end));
            end
            switch header2
                case 'GAZE_COORDS',
                    eyedata.CalibrationAreaOffset=[str2double(dataLine(spcInd(2):spcInd(3)-1)) str2double(dataLine(spcInd(3):spcInd(4)-1))];
                    eyedata.CalibrationAreaSize=[str2double(dataLine(spcInd(4):spcInd(5)-1)) str2double(dataLine(spcInd(5):end))]+1;
                case 'THRESHOLDS',
                    eyedata.Threshold=[str2double(dataLine(spcInd(3):spcInd(4)-1)) str2double(dataLine(spcInd(4):end))];
            end
        elseif length(tabInd)>1
            k=1;
            if strcmp(dataLine(1:tabInd(1)-1),'INPUT')
                inputNo=inputNo+1;
                eyedata.inputs(inputNo,:)=str2num(dataLine(tabInd(1)+1:end));
            end;
            while 1
                if strcmp(dataLine(tabInd(k)+1:tabInd(k+1)-1),'RATE')
                    eyedata.SampleRate=str2double(dataLine(tabInd(k+1)+1:tabInd(k+2)-1));
                    break;
                else
                    k=k+1;
                    if k==length(tabInd)
                        break;
                    end;
                end;
            end;
        end
        
        % If the line starts with a time stamp it should be raw data
    else
        dataNo=dataNo+1;
        %This is just for following the progress
        if mod(dataNo,1000)==0
            fprintf('%i,',dataNo);
        end;
        if mod(dataNo,10000)==0
            fprintf('\n');
        end;
        %%%----------------------------------------------------------------------%%
        %%%       Uncomment the following three lines (if no already uncommented)%%
        %%%       to save it to the rawData field of eyedata...                  %%
        %%%----------------------------------------------------------------------%%
        %          if ~isempty(str2num(dataLine(1:end-3)))
        %              eyedata.rawData(dataNo,:)=str2num(dataLine(1:end-3));
        %          end;
        
    end
end
eyedata.dataNo=dataNo;
% print the header
% eyedata

% nSamples=ceil((lastLine(1)-firstLine(1)+1)/1000*eyedata.SampleRate);
% eyedata.data=zeros(nSamples,length(firstLine));
% eyedata.data(1,:)=firstLine;
%
% disp(sprintf('Estimated number of lines %d. Reading...',nSamples));
%
% lineIndex=2;
% while(1)
%     data=readline(f);
%     if isnan(data(1)),
%         fclose(f);
%         break;
%     else
%         eyedata.data(lineIndex,:)=data;
%         lineIndex=lineIndex+1;
%         %         % uncomment the following to display progress
%         %         if rem(lineIndex,10000)==0,
%         %             disp(n2s([lineIndex data]));
%         %         end
%     end
% end
%
% % done reading, remove the extra rows allocated
% eyedata.data(lineIndex:end,:)=[];
%
%     function numericalData=readline(f)
%         dataLine=fgetl(f);
%         % some Windows versions return empty line from \n\r pair
%         if isempty(dataLine),
%             dataLine=fgetl(f);
%         end
%         if length(dataLine)==1 && dataLine==-1,
%             numericalData=NaN;
%         else
%             numericalData=str2double(dataLine);
%         end
%     end
%
%
end
