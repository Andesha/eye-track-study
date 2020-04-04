% Read eye tracking data into program for excel. 
% cd(folder);
% for i = 1:length(filenames)
% % Load and display progress.  
% curname = filenames{i};
% disp(curname)
% f = xlsread(filenames{i});
% 
% 
% 
% end 


% tsvfiles . 
cd(folder)
f = tdfread(filenames{1}); % statistics toolbox required.

