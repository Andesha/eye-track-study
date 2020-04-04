function outString = tyler_similarity(fixMaps, removedAt, t)
    % Initialization
    outString = '';
    fixMapsSize = size(fixMaps);

    % For each participant...
    for i=1:fixMapsSize(3)
        % If it was one we skipped due to missing data, just write N/A.
        if (ismember(i,removedAt))
            outString = [outString 'N/A,'];
        else
            holder = load(['young_IVT_MovingMedian_avg/young_avg_' num2str(t) '.mat']);
            meanMatrix = holder.meanMatrix;
            % Reshape 2D data to 1D samples and meld together
            % Done with magic matlab syntax
            linearMeanMatrix = meanMatrix(:);
            linearMeanMatrix = linearMeanMatrix';
            % linearMeanMatrix = linearMeanMatrix'; % Just flipping again - trust me
            linearIndivMatrix = fixMaps(:,:,i);
            linearIndivMatrix = linearIndivMatrix(:);
            linearIndivMatrix = linearIndivMatrix';
            statsMatrix = cat(2,linearIndivMatrix',linearMeanMatrix');
            % Do the stats!
            cMat=corr(statsMatrix);
            cMat=cMat(find(triu(ones(size(cMat)),1)));
            r=tanh(nanmean(atanh(cMat)));
            % Save to output string
            outString = [outString num2str(r) ','];
        end  
    end
end
