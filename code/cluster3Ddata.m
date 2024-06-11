function cluster_ID = cluster3Ddata (xyz, tid, radius)
    %% spatial cluster 3D data with DBSCAN
    % 
    % INPUT:
    %   xyz:
    %       N by 3 matrix, as x,y,z coordinates for each localization
    % 
    %   tid: 
    %       N by 1 matrix, as track ID for each localization
    %
    %   radius:        
    %       search radius for DBSCAN
    %
    %
    % OUTPUT:
    %   cluster_ID:
    %       N by 1 matrix, as cluster ID for each localization
    %
    % <Ziqiang.Huang@embl.de
    % 05.03.2024

    traceID = unique(tid);
    
    nDim = 3;
    if (size(xyz, 2)==2)   % 2D data
        nDim = 2;
    elseif (all(xyz(:,3)==0))
        nDim = 2;
        xyz(:, 3) = [];
    else % 3D data, not exclude 1D case
    end


    center = nan(numel(traceID), nDim);
    for i = 1 : numel(traceID)
        center(i, :) = mean(xyz(tid==traceID(i), :));
    end

    cluster_ID = dbscan(center, radius, 1);
    
end