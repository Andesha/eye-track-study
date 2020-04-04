function [globalTime, trialTime]=read_logfile_MORAL(logfile)
% [globalTime, trialTime]=read_logfile(logfile)
% Read presentation logfile and return the first video start time
%   logfile : path to the presentation log file
%   globalTime : start time of the first video from beginning of the session
%   trialTime : trial time 

% V 0.1 26.4.2010  jouko.lampinen@tkk.fi

timeUnits=1/10000;

fid=fopen(logfile,'r');
if fid==-1, 
    error(['Cannot open the file ' logfile]);
end
isfirst=true;
while(1)
    line=fgetl(fid);
    if length(line)==1 && line==-1,
        disp(['Video start event not found in ' logfile]);
        fclose(fid);
        globalTime=0; trialTime=0;
    end
    % field separator is tab, ascii 9
    tabIndices=find(line==9);
    
    if length(tabIndices>=6) 
        eventType=setstr(line(tabIndices(2)+1:tabIndices(3)-1));
        if isfirst && strcmp(eventType,'Pulse')
            strtTime=str2num(line(tabIndices(4)+1:tabIndices(5)-1))*timeUnits;
            isfirst=false;
        end;
        if strcmp(eventType,'Video')
            globalTime=str2num(line(tabIndices(4)+1:tabIndices(5)-1))*timeUnits-strtTime;
            trialTime=str2num(line(tabIndices(5)+1:tabIndices(6)-1))*timeUnits-strtTime;
            fclose(fid);
            return;
        end
    end
end

end
