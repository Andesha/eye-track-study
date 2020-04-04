function out = invalidPoints(data,index)

    checkX = data(:,2);
    checkY = data(:,3);
    
    checkX = checkX(index);
    checkY = checkY(index);
    
    fullSize = length(checkX);
    criteria = fullSize * 0.75;
    
    out = ((sum(isnan(checkX)) + sum(isnan(checkY))) > criteria);
end