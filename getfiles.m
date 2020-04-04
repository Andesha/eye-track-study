%% Gets a list of the file names in a folder. 
%% INPUT: the folder you want. 
%% OUTPUT: the list the names of the files you want to use. 

function FILELIST = getfiles(folder)
    list = dir(folder);   % this will be a structure. 
    % first two values in name field are '.', don't include these. 
    ls = length(list); 
    [FILELIST(1:ls-2).n] = deal(zeros); 
    
    for i = 3:length(list)
        FILELIST(i-2).n = list(i).name; 
    end 
end 
