function factor = computeZscaleFactor (data)
    %% estimate z scale factor of MINFLUX data, based on spatial
    %  dispersion of data
    % 
    % INPUT:
    %   data:
    %       MINFLUX raw data file
    %
    %
    % OUTPUT:
    %   factor:
    %       a value between 0.5 and 1, as the z axis correction factor:
    %       z_corrected = z_original * factor
    %
    %
    % <Ziqiang.Huang@embl.de>
    % 05.03.2024
    
    factor = nan;
    vld = data.vld;
    tid = data.tid(vld);
    id = unique(tid);
    freq = histcounts(tid, [id; inf])';
    num_id = numel(id);

    if isfield(data, 'loc')
        loc = squeeze(data.loc(vld, end, :));
    else
        loc = squeeze(data.itr.loc(vld, end, :));
    end
    x = loc(:,1); y = loc(:,2); z = loc(:,3);

    % estimate z scale for each trace
    zScale_tid = zeros (num_id, 1);
    for i = 1 : num_id
        zScale_tid(i, :) = freq(i) .* estimateZscale ( ...
            x(tid==id(i)), ...
            y(tid==id(i)), ...
            z(tid==id(i)) );
    end

    % remove outlier
    outlier = isoutlier(zScale_tid);
    zScale_tid = zScale_tid(~outlier, :);

    % compute weighted average of z scale factor
    % the weight is number of localizations within a trace
    freq = freq(~outlier, :);
    if (size(zScale_tid, 1) ~= 1)
        zScale_tid = sum(zScale_tid,'omitnan') / sum(freq(~isnan(zScale_tid(:,1))));
    end
    
    % discard result that outside reasonable range
    if (zScale_tid<0.5 || zScale_tid>1.0)
        disp(' Failed to compute z scale factor from loc data!');
    else
        factor = zScale_tid;
    end

end
    
    % cmopute z scale based on spatial dispersion of the localization data
    % on each axis.
    % We use inter-quantile range to reflect the data dispersion
    function scale = estimateZscale (x, y, z)
        dispersion_x = iqr(x, 1);
        dispersion_x = dispersion_x(1);
        dispersion_y = iqr(y, 1);
        dispersion_y = dispersion_y(1);
        dispersion_z = iqr(z, 1);
        dispersion_z = dispersion_z(1);
        scale = geomean([dispersion_x, dispersion_y]) / dispersion_z;
    end