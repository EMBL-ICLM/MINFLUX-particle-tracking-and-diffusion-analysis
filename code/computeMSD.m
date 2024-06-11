function diffusion_data = computeMSD (segments, display)
    %% vectorized method to compute MSD from non-equidistant time intervals
    % 
    % INPUT:
    %   segments:
    %       cell array containing collection of segments, supposedly
    %       belonging to the same cluster.
    %       Each cell element holds tracking data of a track segment.
    %       The tracking data is N by 4 matrix, which are the time stamp,
    %       and x, y, z coordinates of each localization in a given segment.
    %
    %   display:
    %       whether to plot the result
    %
    %
    % OUTPUT:
    %   diffusion_data:
    %       MATLAB struct containing the diffusion analysis result
    %       from different computation modes. It stores raw, regularized,
    %       and weighted mean MSD, tau, and diffusion coefficient.
    %
    %       - MSD_raw:
    %           mean squared displacement at each corresponding time interval
    %
    %   tau:
    %       all possible time intervals
    %
    %   D:
    %       diffusion coefficient, as computed for each MSD-tau pair
    %
    % 
    % the vectorized MSD computation code is adopted from:
    %   Maxime Deforet, May 21 2013. Kehl, a fast (no loop) method to compute MSD 
    %   <https://www.mathworks.com/matlabcentral/fileexchange/41858-kehl-a-fast-no-loop-method-to-compute-msd>
    % the weighted average code is adopted from:
    %   Jean-Yves Tinevez, msdanalyzer
    %   <https://tinevez.github.io/msdanalyzer/>
    %
    %
    % <Ziqiang.Huang@embl.de>
    % 05.03.2024

    if nargin < 2
        display = false;
    end
    
    diffusion_data = struct();
    
    % get number of segments, remove empty segments if there's any
    N = numel(segments);
    nan_entries = cellfun('isempty', segments);
    MSD_raw = cell(N, 1); tau_raw = cell(N, 1); D_raw = cell(N, 1);
    
    % get mode time interval among all segments
    % as the base of the regularized time intervals
    dt_all = cellfun(@(x) diff(x(:, 1)), segments, 'UniformOutput', false);
    dt = vertcat(dt_all{:});
    interval = mode(round(dt*1e4)/1e4); % time interval precision 0.1ms
    MSD_reg = cell(N, 1); tau_reg = cell(N, 1); D_reg = cell(N, 1);

    % computed weighted average values, based on previous raw values,
    % code adapted from msdanalyzer
    % get weights, as number of data points
    segments_all = vertcat(segments);
    ndim = 3;
    %if (all(segments_all(:,3)==0))
    %   ndim = 2;
    %end
    num_all = cellfun(@length, segments_all);
    
    % compute MSD, tau, and D for each segments, both raw and regularized
    for i = 1 : N
        if (nan_entries(i))
            continue;
        end
        xyzt = [segments{i}(:, 2:4), segments{i}(:, 1)];
        [MSD_raw{i}, tau_raw{i}, D_raw{i}] = getMSD (xyzt, ndim);

        xyzt_reg = regularizeTimeInterval(xyzt, interval);
        [MSD_reg{i}, tau_reg{i}, D_reg{i}] = getMSD (xyzt_reg, ndim);
    end
    diffusion_data.MSD_raw = MSD_raw;
    diffusion_data.tau_raw = tau_raw;
    diffusion_data.D_raw = D_raw;
    diffusion_data.MSD_reg = MSD_reg;
    diffusion_data.tau_reg = tau_reg;
    diffusion_data.D_reg = D_reg;
    
    % compute weight averaged MSD, tau, and D
    % combine all raw MSD, and tau
    MSD_all = vertcat(MSD_raw);
    tau_all = vertcat(tau_raw);
    tau_unique = unique(vertcat(tau_all{:}));
    % number of total time intervals
    n_delays = numel(tau_unique);
    % compute weighted mean MSD
    sum_weight          = zeros(n_delays, 1);
    sum_weighted_mean   = zeros(n_delays, 1);
    for i = 1 : numel(MSD_all)
        m = MSD_all{i};
        t = tau_all{i};
        n = num_all(i);
        % Find common indices
        [~, index_in_all_delays, ~] = intersect(tau_unique, t);
        % Accumulate
        sum_weight(index_in_all_delays)           = sum_weight(index_in_all_delays)         + n;
        sum_weighted_mean(index_in_all_delays)    = sum_weighted_mean(index_in_all_delays)  + m .* n;
    end
    diffusion_data.MSD_mean = sum_weighted_mean ./ sum_weight;
    diffusion_data.tau_mean = tau_unique;
    diffusion_data.D_mean   = computeDiffCoeff(diffusion_data.MSD_mean, diffusion_data.tau_mean, ndim);


    if (display) 
        figure, plot(tau_raw, 1e12* MSD_raw);
        hold on;
        plot(tau_reg, 1e12* MSD_reg, 'r');
        hold on;
        plot(tau_mean, 1e12* MSD_mean, 'LineWidth', 5, 'Color', [.2 .2 .2 .9]);
    end

end
    
    % compute MSD with vectorized approach
    function [MSD, tau, D] = getMSD (xyzt, ndim)
        xyzt(:,end) = xyzt(:,end)-min(xyzt(:,end)); % re-zero time data
        N =  size(xyzt,1);
        [ i, j ] = find(triu(ones(N), 1));
        D = nan(1, N^2);
        D( i + N*(j-1) ) = (sum(abs( xyzt(i,1:end-1) - xyzt(j,1:end-1) ).^2, 2));
        dt = nan(1, N^2);
        dt( i + N*(j-1) ) = -( xyzt(i,end) - xyzt(j,end) );
        dt = round(dt*1e4)/1e4;     % round to 0.1 millisecond precision
        [DT,idx]=sort(dt(:));
        DD=D(idx);
        First_idx=find(DT-circshift(DT,1)~=0);
        Last_idx=find(DT-circshift(DT,-1)~=0);
        C=cumsum([0,DD]);
        MSD=(C(Last_idx+1)-C(First_idx))'./(Last_idx-First_idx+1); 
        tau=DT(First_idx);
        MSD(isnan(MSD))=[];
        tau(isnan(tau))=[];
        D = computeDiffCoeff (MSD, tau, ndim);
    end

    % compute Diffusion Coefficient from MSD-tau pair 
    function D = computeDiffCoeff (MSD, tau, ndim)
        D = MSD./tau / (2*ndim);
    end

    % regularize time interval
    function xyzt_reg = regularizeTimeInterval(xyzt, interval)
        N = size(xyzt, 1);
        T = xyzt(:, end);
        T_reg = (T(1) : interval : T(N))';
        xyzt_reg = interp1(T, xyzt, T_reg);
    end

