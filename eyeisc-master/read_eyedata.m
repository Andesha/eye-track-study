function eyedata=read_eyedata(eyefile)
%eyedata=read_eyedata(eyefile)
%  Read eye tracking data
%   eyefile : path to the eye track file
%   eyedata : struct with following fields with example values
%               FileVersion: '	2'
%                Fileformat: '	2558'
%                   Subject: '	'
%                      Date: '	April 09, 2010'
%               Description: '	'
%            NumberOfPoints: '	n.a.'
%     CalibrationAreaOffset: [0 0]
%       CalibrationAreaSize: [1024 768]
%                SampleRate: 60
%                fieldNames: {1x10 cell}
%                      data: [253623x10 double]
% The first six fields are returned as unconverted text,
% the next three fields are converted to numerical values,
%  .fieldNames is cell array containing the fields names as listed in the file
%  .data is unprocessed data matrix 

% V 0.1 26.4.2010 jouko.lampinen@tkk.fi 

f=fopen(eyefile,'r');
if f==-1, 
    error(['Cannot open the file ' eyefile]);
end

% find the eof and read the last line to estimate the number of lines
% in the file, to preallocate the data matrix
fseek(f,-2,'eof');
while(1)
    byte=fread(f,1,'uchar');
    if byte==10, 
        break
    end
    fseek(f,-2,'cof');
end
lastLine=readline(f);

fseek(f,0,'bof');

% read the header 
while(1)
    dataLine=fgetl(f);
    % some Windows versions return empty line from \n\r pair
    if isempty(dataLine),
        dataLine=fgetl(f);
    end
    if dataLine(1)=='#'
        j=find(dataLine==':');
        if ~isempty(j),
            header=setstr(dataLine(2:j(1)-1));
            switch header
                case 'FileVersion',
                    eyedata.FileVersion=dataLine(j(1)+1:end);
                case 'Fileformat',
                    eyedata.Fileformat=dataLine(j(1)+1:end);
                case 'Subject',
                    eyedata.Subject=dataLine(j(1)+1:end);
                case 'Date',
                    eyedata.Date=dataLine(j(1)+1:end);
                case 'Description',
                    eyedata.Description=dataLine(j(1)+1:end);
                case '# of Pts Recorded',
                    eyedata.NumberOfPoints=dataLine(j(1)+1:end);
                case 'Offset Of Calibration Area',
                    eyedata.CalibrationAreaOffset=str2num(dataLine(j(1)+1:end));
                case 'Size Of Calibration Area',
                    eyedata.CalibrationAreaSize=str2num(dataLine(j(1)+1:end));
                case 'Sample Rate',
                    eyedata.SampleRate=str2num(dataLine(j(1)+1:end));
                otherwise
                    disp(['Unknown header field:' header]);
            end
        else 
            % no : in the header            
            if length(dataLine)>10,
                % exract the field names
                dataLine=dataLine(2:end);
                j=[0 find(dataLine==9) length(dataLine)+1];
                for k=1:length(j)-1,
                    label{k}=dataLine(j(k)+1:j(k+1)-1);
                end
                eyedata.fieldNames=label;
            end
        end
    else
        firstLine=str2num(dataLine);
        break;
    end
end

% print the header
% eyedata

nSamples=ceil((lastLine(1)-firstLine(1)+1)/1000*eyedata.SampleRate);
eyedata.data=zeros(nSamples,length(firstLine));
eyedata.data(1,:)=firstLine;

disp(sprintf('Estimated number of lines %d. Reading...',nSamples));

lineIndex=2;
while(1)
    data=readline(f);
    if isnan(data(1)),
        fclose(f);
        break;
    else
        eyedata.data(lineIndex,:)=data;
        lineIndex=lineIndex+1;
%         % uncomment the following to display progress
%         if rem(lineIndex,10000)==0,
%             disp(n2s([lineIndex data]));
%         end
    end
end

% done reading, remove the extra rows allocated
eyedata.data(lineIndex:end,:)=[];

function numericalData=readline(f)
    dataLine=fgetl(f);
    % some Windows versions return empty line from \n\r pair
    if isempty(dataLine),
        dataLine=fgetl(f);
    end
    if length(dataLine)==1 && dataLine==-1,
        numericalData=NaN;
    else
        numericalData=str2num(dataLine);
    end
end


end
