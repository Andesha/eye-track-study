function [globalTime, trialTime]=read_logfile(logfile)
% [globalTime, trialTime]=read_logfile(logfile)
% Read presentation logfile and return the first video start time
%   logfile : path to the presentation log file
%   globalTime : start time of the first video from beginning of the session (edit: from the first pulse)
%   trialTime : trial time 

% V 0.1 26.4.2010  jouko.lampinen@tkk.fi
% Slight edits 14.8.2014 juha.lahnakoski@aalto.fi

timeUnits=1/10000;

f=fopen(logfile,'r');
if f==-1, 
    error(['Cannot open the file ' logfile]);
end
firstPulse=true;
while(1)
    line=fgetl(f);
    if length(line)==1 && line==-1,
        disp(['Video start event not found in ' logfile]);
        fclose(f);
        globalTime=0; trialTime=0;
    end
    % field separator is tab, ascii 9
    tabIndices=find(line==9);

    if length(tabIndices>=6) 
        eventType=setstr(line(tabIndices(2)+1:tabIndices(3)-1));
	if strcmp(eventType,'Pulse') && firstPulse
		firstPulse=false;
		startTime=str2num(line(tabIndices(4)+1:tabIndices(5)-1));
	end
        if strcmp(eventType,'Video')
            code=setstr(line(tabIndices(3)+1:tabIndices(4)-1));
            %if strcmp(code,'video')
                globalTime=(str2num(line(tabIndices(4)+1:tabIndices(5)-1))-startTime)*timeUnits;
                trialTime=str2num(line(tabIndices(5)+1:tabIndices(6)-1))*timeUnits;
                fclose(f);
                return;
            %end
        end
    end
end

end
