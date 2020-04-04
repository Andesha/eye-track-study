function outString = tyler_similarity(fixMaps, removedAt, t)
    % Initialization
    outString = '';
    fixMapsSize = size(fixMaps);
    summate = zeros(fixMapsSize(1),fixMapsSize(2));
    % For each participant...
    for i=1:fixMapsSize(3)
        % If it was one we skipped due to missing data, just write N/A.
        if (ismember(i,removedAt))
            outString = '';
        else
            % Pre allocate just to be safe
  
            % Add everything together
            summate = summate + fixMaps(:,:,i);
            % Scale the summation to become the mean
            outString = '';
        end  
    end
    meanMatrix = summate./(fixMapsSize(3) - 1);
    save(['young_avg_' num2str(t) '.mat'],'meanMatrix')
end
