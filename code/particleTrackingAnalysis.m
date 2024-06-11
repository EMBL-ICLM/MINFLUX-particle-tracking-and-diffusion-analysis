function particle_data = particleTrackingAnalysis(matFilePath, lengthRemoval, zScale)

    %% a MATLAB script for processing, clustering, and analysis of MINFLUX tracking data.
    % 
    % INPUT:
    %   matFilePath:
    %       the system path to the .mat format MINFLUX raw data
    % 
    %   lengthRemoval: 
    %       the number of data point to be discarded at the
    %       beginning and the end of each trace (track)
    %
    %   zScale:        
    %       z axis correction factor
    %
    %
    % OUTPUT:
    %   particle_data:
    %       MATLAB struct containing the analysis result
    %       - each row corresponding to a cluster from DBSCAN
    %       - the columns (fields) corresponding to:
    %        01, cluster        : cluster ID
    %        02, nLocalization  : number of data points
    %        03, trace          : trace ID(s) in this cluster
    %        04, axis range     : bounding box of the cluster
    %        05, coordinates    : localization coordinates (of each trace)
    %        06, time           : time stamp (of each trace)
    %        07, displacement   : calculated displacement in between each adjacent localization 
    %        08, velocity       : calculated instantenous velocity in between each adjacent localization 
    %        09, sphericity     : a measure of volume to area fraction, as an quantification of how much 3D cloud is shpere-like 
    %        10, single particle: whether cluster could from a single silicon core
    %        11, sphereFitting  : whether sphere shell fitting is succesful
    %        12, fittingError   : normalized error of the least-square sphere fitting
    %        13, center         : fitted sphere center
    %        14, diameter       : fitted sphere diameter
    %        15, segments       : track segments, that broken around 2.5* min(dt)
    %        16, MSD_raw        : computed MSD values for each segment
    %        17, tau_raw        : all possible time intervals in each segment
    %        18, D_raw          : computed diffusion coefficient values for each segment
    %        19, MSD_reg        : regularized MSD values for each segment
    %        20, tau_reg        : regularized time intervals in each segment
    %        21, D_reg          : computed regularized diffusion coefficient values for each segment
    %        22, MSD_mean       : weighted average MSD values for all data in current cluster
    %        23, tau_mean       : all possible time intervals for all data in current cluster
    %        24, D_mean         : computed weighted average diffusion coefficient values for all data in current cluster
    %
    %
    % <Ziqiang.Huang@embl.de>
    % 05.03.2024
    
    % need at least two data points in track after filtering, to perform
    % diffusion analysis
    minLoc = 2*lengthRemoval + 2;

    % load MINFLUX MATLAB format data
    data = load(matFilePath);
    vld = data.vld;
    
    % get localization coordinates, raw data
    if ~isfield(data, 'loc')    % account for Aberrior format change
        xyz = squeeze(data.itr.loc(vld, end, :));
    else
        xyz = squeeze(data.loc(vld, end, :));
    end

    % remove marginal traces
    margin = 0.01; % filter out data within 1% of the border
    x = xyz(:, 1); y = xyz(:, 2);
    xmin = min(x); xmax = max(x); xrange = xmax - xmin;
    ymin = min(y); ymax = max(y); yrange = ymax - ymin;
    xmargin = margin * xrange; ymargin = margin * yrange;
    nearBorder = x<xmin+xmargin | x>xmax-xmargin | y<ymin+ymargin | y>ymax-ymargin; % if data located within 1% to the border
    
    % get track (trace) ID
    tid = data.tid(vld);
    trace_ID = unique(tid(~nearBorder))'; % trace_ID is track ID, and the same from this point on

    % remove the beginning and end of each trace
    % make use of vld attribute, as if the removed data points are invalid
    idx_1 = arrayfun(@(x) find(tid==x, 1, 'first'), trace_ID); 
    TF_origin = vld;
    TF_origin(idx_1) = 0;
    TF_origin(end)=0;
    TF = ones(2*lengthRemoval, length(vld));
    for K = -lengthRemoval : lengthRemoval-1
        TF(K+lengthRemoval+1, :) = circshift(TF_origin, K);
    end
    TF_final = sum(TF, 1)==2*lengthRemoval;
    vld = vld & TF_final;
    
    % update trace ID array with valid localizations
    tid = data.tid(vld);
    trace_ID = unique(tid)';

    % get time stamp of data
    t = data.tim'; % no filter of time value, to compute delay
    maxDt = 2.5 * min(diff(t));   %!!! ALLOW ONLY 2 STICKNESS
    t = t(vld);

    % get valid localization coordinates, raw
    % get XYZ coordinates data, and scale Z if required
    xyz = squeeze(data.loc(vld, end, :));
    ndims = 3;
    if (all(xyz(:,3)==0))   % 2D data
        ndims = 2;
        xyz(:,3) = [];
    else
        xyz(:,3) = xyz(:,3) * zScale;
    end
    

    % cluster xyz data with DBSCAN
    clusterRadius = 2e-7; % cluster radius set to be the expected diameter of the silicon particle
    cluster_ID = cluster3Ddata (xyz, tid, clusterRadius);
    [~, ~, ic] = unique(tid);
    cid = arrayfun(@(x) cluster_ID(x), ic);
    nCluster = numel(unique(cid));
    
    % prepare result for each cluster
    particle_data = struct(); idx = 0;
    for cluster_idx = 1 : nCluster
        % processing cluster: i
        currentCluster = cid == cluster_idx; % locate data belonging to the current trace ID
        nLocalization = sum(currentCluster);
        if (nLocalization < minLoc || nLocalization < 4)
            continue;
        end
        idx = idx + 1;

        xyz_cluster = xyz(currentCluster, :);
        t_cluster = t(currentCluster);
        
        particle_data(idx).cluster = cluster_idx;
        particle_data(idx).nLocalization = nLocalization;
        particle_data(idx).trace = trace_ID(cluster_ID==cluster_idx, 1);

        nLocPerTrace = arrayfun(@(x) sum(tid==x), particle_data(idx).trace);

        % compute coordinates, time, displacement, and velocity, 
        % trace-wise, per cluster
        particle_data(idx).axis_range = [min(xyz_cluster); max(xyz_cluster)];
        particle_data(idx).coordinates = mat2cell(xyz_cluster, nLocPerTrace);
        particle_data(idx).time = mat2cell(t_cluster, nLocPerTrace);
        time_diff = cellfun(@(x) diff(x), particle_data(idx).time, 'UniformOutput', false);
        particle_data(idx).displacement = cellfun(@(x) vecnorm( diff(x), 2, 2), particle_data(idx).coordinates, 'UniformOutput', false);
        particle_data(idx).velocity = cellfun(@(x,y) x./y, particle_data(idx).displacement, time_diff, 'UniformOutput', false);
        
        particle_data(idx).sphericity = computeSphericity(xyz_cluster);
        particle_data(idx).single_particle = false;
        % fit sphere, need at least 4 points in 3D
        particle_data(idx).sphereFitting = "NA";
        if (ndims == 3)
            lastwarn('');
            [center, radii, chi2] = computeSphereFit(xyz_cluster);
            [warnMsg, ~] = lastwarn;
            if ~isempty(warnMsg)    %fprintf('sphere fit not possible for track %0d.\n\n', i);
                particle_data(idx).sphereFitting = "Fail";
            else
                particle_data(idx).sphereFitting = "Success";
            end
            particle_data(idx).fittingError = sqrt(chi2/particle_data(idx).nLocalization);
            particle_data(idx).center = center';
            particle_data(idx).diameter = radii(1)*2;
            if (particle_data(idx).diameter <= 2e-7)
                particle_data(idx).single_particle = true;
            end
        else
            particle_data(idx).sphereFitting = "2D";
            particle_data(idx).fittingError = nan;
            particle_data(idx).center = nan(1, 3);
            particle_data(idx).diameter = nan;
        end
        
        % compute MSD and diffusion coefficient, on track segments
        txyz = [t_cluster xyz_cluster];
        largeIntervel_idx = [0 find(diff(t_cluster)>maxDt)' size(txyz, 1)];
        didx = diff(largeIntervel_idx);
        out = mat2cell(txyz, didx);
        segments = out(didx>3, 1);
        particle_data(idx).segments = segments;
        nSeg = numel(segments);
        particle_data(idx).MSD_raw = cell(nSeg, 1); 
        particle_data(idx).tau_raw = cell(nSeg, 1);
        particle_data(idx).D_raw = cell(nSeg, 1);
        weight = cellfun(@(x) size(x, 1), segments);
        if (max(weight)>1e4)
            idx = idx - 1;
            continue;
        end

        diffusion_data = computeMSD (segments, false);
        particle_data(idx).MSD_raw  = diffusion_data.MSD_raw;
        particle_data(idx).tau_raw  = diffusion_data.tau_raw;
        particle_data(idx).D_raw    = diffusion_data.D_raw;
        particle_data(idx).MSD_reg  = diffusion_data.MSD_reg;
        particle_data(idx).tau_reg  = diffusion_data.tau_reg;
        particle_data(idx).D_reg    = diffusion_data.D_reg;
        particle_data(idx).MSD_mean = diffusion_data.MSD_mean;
        particle_data(idx).tau_mean = diffusion_data.tau_mean;
        particle_data(idx).D_mean   = diffusion_data.D_mean;
        
    end
    % sort data so that better sphere fitting results are pushed to front
    [~,index] = sortrows([particle_data.sphericity].'); particle_data = particle_data(index(end:-1:1)); clear index

    % save result
    [filepath, name, ~] = fileparts(matFilePath);
    matDataFile = strcat(name, '_result.mat');
    save(fullfile(filepath, matDataFile), 'particle_data', '-v7.3');

    
end