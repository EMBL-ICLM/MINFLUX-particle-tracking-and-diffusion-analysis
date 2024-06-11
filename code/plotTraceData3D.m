function plotTraceData3D (particle_data, index)
%% main function execution starts here
    % input argument control, select 1st trace if not specified
    if nargin < 2
        index = 1;
    else
        index = index(1); % control the index to be single value
    end
    
    % initiate global variables
    globalVar = struct();

    % initiate cluster specific variables
    cluster_selected = struct();    % selected cluster structure
    cluster_ID_selected = [];       % cluster ID of the selected cluster, single cluster allowed
    cluster_idx_selected = [];      % logical array corresponding to loc in selected cluster
    globalVar.cluster_array_selected = [];

    globalVar.nTraces_selected = [];                   % number of traces in selected cluster
    trace_ID_selected = [];                  % tracy ID list in selected cluster
    globalVar.trace_array_selected = [];
    globalVar.traceColor_selected = [];

    trace_ID_segment = [];          % trace ID for each segments in selected cluster
    nSeg = [];
    fig_result_name = '';
    %trace_ID_segment_array = [];    % trace ID for each MSD/tau/D entry
    %plotMSDColor_selected = [];
    % initate plot related variables


        
    %% get global variables from data
    %[globalVar.nLoc_all, globalVar.nDim, ~, globalVar.cluster_ID_all, xyz_all, xyz_selected, globalVar.cluster_array_all, globalVar.clusterColor] = initGlobalVar (particle_data);
    initGlobalVar (particle_data);
    initSelectedCluster (index);
    
    %% prepare overview and result figures
    
    % figure 901 as the overview figure
    if ~ishandle(901)
        fig_overview = figure(901);
        fig_overview.Name = 'Localization Data Overview';
        fig_overview.NumberTitle = 'off';
        fig_overview.Position = [100 400 1000 500];
    else
        fig_overview = findobj( 'Type', 'Figure', 'Number', 901);
    end

    plotOverview(); % return the axis and scatter plot obj of the overview plot
    colorOverview();

    %  figure 902 as the measurement result figure for selected trace(s)
    if ~ishandle(902)
        fig_result = figure(902);
        fig_result.Position = [300 100 1200 600];
        fig_result.Name = fig_result_name;%strcat('Measurement Result: cluster ', num2str(index), ' , sphericity: ', num2str(cluster_selected.sphericity));
        fig_result.NumberTitle = 'off';
    else
        fig_result = findobj( 'Type', 'Figure', 'Number', 902);
    end

    plotCluster(); %#ok<*SETNU> 
%% main function execuation ends here


%% member functions

    % load global variables from particle_data
    function initGlobalVar (data)

        globalVar.nDim = getDimension(data);                            % dimension of data
        globalVar.nCluster = numel(data);                               % number of clusters in data
        globalVar.cluster_ID_all = int32(vertcat(data(:).cluster));     % cluster ID list in data
        
        globalVar.loc_trace = vertcat(data(:).coordinates);             % localization coordinates as cell array, each cell belongs to a single trace
        globalVar.nTrace = numel(globalVar.loc_trace);                  % number of traces in data
        globalVar.trace_ID_all = vertcat(data(:).trace);                % trace ID list in data

        globalVar.loc_all = 1e9 * vertcat(globalVar.loc_trace{:});      % all localization coordinates of data
        globalVar.nLoc_all = size(globalVar.loc_all, 1);                % total number of localization coordinates
        globalVar.nLocPerCluster = vertcat(data.nLocalization);         % number of localization in each cluster
        globalVar.nLocPerTrace = cellfun(@(x) size(x, 1), globalVar.loc_trace);	% number of localization in each trace

        globalVar.cluster_array = repelem(globalVar.cluster_ID_all, globalVar.nLocPerCluster, 1);	% cluster ID array, denote each loc's cluster ID 
        globalVar.trace_array = repelem(globalVar.trace_ID_all, globalVar.nLocPerTrace, 1);         % tracy ID array, denote each loc's trace ID

        globalVar.clusterColor = generatePlotColors(globalVar.nCluster);	% cluster plot color scheme, use glasbey on dark to discern 
        globalVar.traceColor = generatePlotColors(globalVar.nTrace);        % trace plot color scheme, use glasbey on dark to discern
        globalVar.colorMode = 1;            % color mode for overview plot: 0, color by cluster; 1, color by seleted cluster
        
        globalVar.gui_overview = struct();  % store overview plot gui data
        globalVar.gui_result = struct();    % store result plot gui data
  
    end

    % load selected cluster ID and array
    function initSelectedCluster (index)
        cluster_selected = particle_data(index);
        cluster_ID_selected = cluster_selected.cluster;
        cluster_idx_selected = globalVar.cluster_array == cluster_ID_selected;
        globalVar.cluster_array_selected = int8(cluster_idx_selected);
        trace_ID_selected = cluster_selected.trace;  % trace IDs in cluster
        globalVar.nTraces_selected = numel(trace_ID_selected);
        %globalVar.nLocPerTrace = cellfun(@(x) size(x, 1), cluster_selected.coordinates);
        globalVar.trace_array_selected = globalVar.trace_array(cluster_idx_selected, :); % loc in cluster, marked with trace ID
        globalVar.traceColor_selected = generatePlotColors (globalVar.nTraces_selected);
        nSeg = numel(cluster_selected.segments);
        fig_result_name = strcat('Measurement Result: cluster ID:', num2str(cluster_ID_selected), ...
            ' , sphericity: ', num2str(cluster_selected.sphericity), ...
            ' , diameter: ', num2str(cluster_selected.diameter*1e9, '%3.0f'), ...
            ' nm , single particle: ');
        if (cluster_selected.single_particle)
            fig_result_name = strcat(fig_result_name, 'yes');
        else
            fig_result_name = strcat(fig_result_name, 'no');
        end
        fig_result_name = strcat(fig_result_name, ', segments: ', num2str(nSeg));
        %plotMSDColor_selected = [repmat([0 .5 .1 ], nSeg, 1); [.7 .1 .7]; [.1 .1 .1]];
        %locateSegments();
        %ssmakeMSDPlotColors();
    end

    % plot overview scatter
    function plotOverview ()
        set(0, 'currentfigure', fig_overview);
        %clf(fig_overview, 'reset');
        globalVar.gui_overview = struct();
        globalVar.gui_overview.ax = fig_overview.CurrentAxes;
        set(globalVar.gui_overview.ax, 'NextPlot', 'replacechildren');

        if (globalVar.nDim == 3)
            globalVar.gui_overview.sc = scatter3(globalVar.loc_all(:,1), globalVar.loc_all(:,2), globalVar.loc_all(:,3), 1, globalVar.cluster_array);
        elseif (globalVar.nDim == 2)
            globalVar.gui_overview.sc = scatter(globalVar.loc_all(:,1), globalVar.loc_all(:,2), 1, globalVar.cluster_array);
        else
            error('data dimension wrong!');
        end
        % add cluster ID row to the data tip
        globalVar.cluster_array_datatip = arrayfun(@(x) double(repelem(globalVar.cluster_ID_all(x), globalVar.nLocPerCluster(x))), 1:globalVar.nCluster, 'UniformOutput', false);
        globalVar.cluster_array_datatip = horzcat(globalVar.cluster_array_datatip{:})';
        globalVar.gui_overview.sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Cluster', globalVar.cluster_array_datatip);
        %sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace', globalVar.trace_array);
        title('colored by cluster ID');
        xlabel('X (nm)', 'FontSize', 20);
        ylabel('Y (nm)', 'FontSize', 20);
        zlabel('Z (nm)', 'FontSize', 20);
        axis equal;
        % adding UI components
        %ui_Load = 
        uicontrol(  'Parent',fig_overview, ...
                    'Style','pushbutton',...
                    'Units','normalized',...
                    'Position',[0.1 0.02 0.08 0.03],...
                    'Fontsize',10,...
                    'string','load',...
                    'Visible','on',...
                    'Callback',{@loadResultData});
        %ui_Color =
        uicontrol(  'Parent',fig_overview, ...
                    'Style','pushbutton',...
                    'Units','normalized',...
                    'Position',[0.2 0.02 0.08 0.03],...
                    'Fontsize',10,...
                    'string','color',...
                    'Visible','on',...
                    'Callback',{@changeOverviewColor});
        % create a pushbtton on overview figure for data selection
        %ui_Selection =
        uicontrol(  'Parent',fig_overview, ...
                    'Style','pushbutton',...
                    'Units','normalized',...
                    'Position',[0.3 0.02 0.08 0.03],...
                    'Fontsize',10,...
                    'string','select',...
                    'Visible','on',...
                    'Callback',{@BoxButtonPushed});
        globalVar.gui_overview.ax = gca;
    end

    % plot selected trace measurement results
    function plotCluster()
        %% Update the measurement result figure window
        % refresh figure
        set(0, 'currentfigure', fig_result);
        fig_result.Name = fig_result_name; %['Measurement Result: cluster ', num2str(cluster_ID_selected)];
        clf(fig_result);
        
        % create struct to store result figure window GUI data
        globalVar.gui_result = struct();
        
        %% left panel: scatter plot 
        globalVar.gui_result.ax1 = subplot('Position', [0.05, 0.15, 0.25, 0.75]);
        plotClusterScatter (globalVar.gui_result.ax1);

        %% right upper panel: from left to right: displacement, velocity, time interval
        % displacement over time
        globalVar.gui_result.ax2 = subplot(2,4,2);
        plotDisplacement (globalVar.gui_result.ax2);
        % instantaneous velocity plot with RMS
        globalVar.gui_result.ax3 = subplot(2,4,3);
        plotVelocity (globalVar.gui_result.ax3);
        % instantaneous time interval with RMS
        globalVar.gui_result.ax4 = subplot(2,4,4);
        plotInterval (globalVar.gui_result.ax4);
        %% right lower panel: from left to right: MSD plot, log(D) histogram, log(D) cumulative frequency
        % MSD curve
        globalVar.gui_result.ax5 = subplot(2,4,6);
        plotMSD (globalVar.gui_result.ax5);

        % log(D) histogram
        globalVar.gui_result.ax6 = subplot(2,4,7);
        plotDiffusionHisto (globalVar.gui_result.ax6);%, D_mean, D_reg);

        % log(D) cumulative frequency

        globalVar.gui_result.ax7 = subplot(2,4,8);
        plotDiffusionCumFreq (globalVar.gui_result.ax7);%, D_mean, D_reg);


        % add legend to result figure
        legendTxt = arrayfun(@(x) ['trace: ', num2str(trace_ID_selected(x))], (1:globalVar.nTraces_selected), 'UniformOutput', false);
        legend(globalVar.gui_result.ax2, legendTxt, 'AutoUpdate', 'off', 'Location', 'northwest');

        % add UI components: time slider
        globalVar.gui_result.vis_pos = 1;
        globalVar.gui_result.vis_rng = 10;

        pos_1 = globalVar.gui_result.ax1.Position;
        %pos_4 = globalVar.gui_result.ax4.Position;
        pos_5 = globalVar.gui_result.ax5.Position;
        pos_6 = globalVar.gui_result.ax6.Position;
        pos_7 = globalVar.gui_result.ax7.Position;
        
        % get the longest track, apply time slider around it
        time = cluster_selected.time;
        N = cellfun(@(x) numel(x), time);
        [~, idx] = max(N);
        range = numel(time{idx});

        % create a time slider to browse the trace scatter data
        %ui_Slider =
        uicontrol(  'Parent',fig_result, ...
                    'Style','slider',...
                    'Min', 1,...
                    'Max', range,...
                    'SliderStep', [1/range 10/range],...
                    'Value', 1,...
                    'Units','normalized',...
                    'Position',[pos_1(1)-0.01 pos_1(2)-0.1 pos_1(3)+0.02 0.02],...
                    'Fontsize',10,...
                    'string','time',...
                    'Visible','on',...
                    'Callback',{@framePosUpdate});
        % create a edit field for visible range of time slider 
        uicontrol(  'Parent',fig_result, ...
                    'Style','text',...
                    'Units','normalized',...
                    'Position',[pos_1(1)+pos_1(3)-0.13 pos_1(2)-0.14 0.11 0.03],...
                    'Fontsize',10,...
                    'string','visible frame range:',...
                    'Visible','on');
        %ui_range = 
        uicontrol(  'Parent',fig_result, ...
                    'Style','edit',...
                    'Units','normalized',...
                    'Position',[pos_1(1)+pos_1(3)-0.02 pos_1(2)-0.14 0.03 0.03],...
                    'Fontsize',10,...
                    'string','10',...
                    'Visible','on',...
                    'Callback',{@frameRangeUpdate});

        % create a pushbtton underneath time interval plot for time range selection
        %ui_TimeSelection =
%         uicontrol(  'Parent',fig_result, ...
%                     'Style','pushbutton',...
%                     'Units','normalized',...
%                     'Position',[pos_4(1) pos_4(2)-0.08 0.07 0.02],....
%                     'Fontsize',10,...
%                     'string','select time window',...
%                     'Visible','on',...
%                     'Callback',{@selectTimeWindow, time, xyz, ax1, ax2, ax3, ax4});

        % create a pushbtton underneath time interval plot for time range selection
        %ui_AveragePlot =
        uicontrol(  'Parent',fig_result, ...
                    'Style','checkbox',...
                    'Units','normalized',...
                    'Position',[pos_5(1)+0.03 pos_5(2)-0.1 0.07 0.03],....
                    'Fontsize',10,...
                    'string','show raw',...
                    'Visible','on',...
                    'Value', 1,...
                    'Callback',{@showRawMSD});
        uicontrol(  'Parent',fig_result, ...
                    'Style','checkbox',...
                    'Units','normalized',...
                    'Position',[pos_6(1)+0.03 pos_6(2)-0.1 0.1 0.03],....
                    'Fontsize',10,...
                    'string','show equalized',...
                    'Visible','on',...
                    'Callback',{@showRegMSD});
        uicontrol(  'Parent',fig_result, ...
                    'Style','checkbox',...
                    'Units','normalized',...
                    'Position',[pos_7(1)+0.03 pos_7(2)-0.1 0.15 0.03],....
                    'Fontsize',10,...
                    'string','show group average',...
                    'Visible','on',...
                    'Callback',{@showMeanMSD});
    end

    % plot trace scatter plot
    function plotClusterScatter (axes)
        cla(axes);
        set(axes, 'ColorOrder', globalVar.traceColor_selected, 'NextPlot', 'add');
        globalVar.gui_result.sc1_1 = cell(globalVar.nTraces_selected, 1); globalVar.gui_result.sf1_1 = [];
        xyz = cluster_selected.coordinates;       
        for i = 1 : globalVar.nTraces_selected
            xyz_trace = 1e9*(xyz{i}-cluster_selected.center);
            if (globalVar.nDim == 3)
                
                globalVar.gui_result.sc1_1{i} = plot3(axes, xyz_trace(:,1), xyz_trace(:,2), xyz_trace(:,3), 'LineWidth', 0.2, 'Marker', '.', 'MarkerSize', 10);% ,'DisplayName', strcat('trace-', num2str(trace_ID{j})));
  
            elseif (globalVar.nDim == 2)
                globalVar.gui_result.sc1_1{i} = scatter(axes, xyz_trace(:,1), xyz_trace(:,2), '.');% ,'DisplayName', strcat('trace-', num2str(trace_ID{j})));
            end
        end
        if ( ~strcmp(cluster_selected.sphereFitting, "Fail") && cluster_selected.sphericity>0.8)          % create sphere surface
            globalVar.gui_result.sf1_1 = plotSphereSurface (axes, [0 0 0], 1e9* cluster_selected.diameter / 2, 'g', 'on');
        end
        
        % highlighted trace in magenta
        globalVar.gui_result.highlight_trace = plot(axes, xyz{1}(:,1), xyz{1}(:,2));
        globalVar.gui_result.highlight_head = plot(axes, xyz{1}(:,1), xyz{1}(:,2));
        globalVar.gui_result.highlight_trace.Visible = 'off';
        globalVar.gui_result.highlight_head.Visible = 'off';
        
        % set title and x,y,z label of axes
        grid on;
        title('localization scatter plot', 'FontSize', 20);
        xlabel('X (nm)', 'FontSize', 20);
        ylabel('Y (nm)', 'FontSize', 20);
        zlabel('Z (nm)', 'FontSize', 20);
        axis equal;
    end
    
    % plot displacement
    function plotDisplacement (axes)
        cla(axes);
        set(axes, 'ColorOrder', globalVar.traceColor_selected, 'NextPlot', 'add');
        globalVar.gui_result.p2 = cell(globalVar.nTraces_selected, 1);
        time = cluster_selected.time;
        for j = 1 : globalVar.nTraces_selected
            time{j} = time{j} - time{j}(1);
            globalVar.gui_result.p2{j} = plot(axes, time{j}(1:end-1), 1e6*cumsum(cluster_selected.displacement{j}));
        end
        globalVar.gui_result.bar2 = bar(axes, 1, 0);
        globalVar.gui_result.bar2.Visible = 'off';
        title( 'displacement over time');
        xlabel('time (second)');
        ylabel('displacement (um)');
        axis tight;
    end

    % plot instantaneous velocity
    function plotVelocity (axes)
        cla(axes);
        set(axes, 'ColorOrder', globalVar.traceColor_selected, 'NextPlot', 'add');
        globalVar.gui_result.p3 = cell(globalVar.nTraces_selected, 1);
        time = cluster_selected.time;
        for i = 1 : globalVar.nTraces_selected
            time{i} = time{i} - time{i}(1);
            velocity_rms = RMS (cluster_selected.velocity{i}, 1e2);
            globalVar.gui_result.p3{i} = plot(axes, time{i}(1:end-1), 1e6*velocity_rms);
        end
        globalVar.gui_result.bar3 = bar(axes, 1, 0);
        globalVar.gui_result.bar3.Visible = 'off';
        title( 'instantaneous velocity ');
        xlabel('time (second)') ;
        ylabel('velocity (um/s)') ;
        axis tight;
        yl = ylim;
        ylim([0 yl(2)*2]);
    end

    % plot instantaneous velocity, against time
    function plotInterval (axes)
        cla(axes);
        set(axes, 'ColorOrder', globalVar.traceColor_selected, 'NextPlot', 'add');
        globalVar.gui_result.p4 = cell(globalVar.nTraces_selected, 1); %p2 = cell(nPlot, 1);
        time = cluster_selected.time;
        for i = 1 : globalVar.nTraces_selected
            time{i} = time{i} - time{i}(1);
            globalVar.gui_result.p4{i} = plot(axes, time{i}(1:end-1), 1e3*diff(time{i}), '.');
        end
        globalVar.gui_result.bar4 = bar(axes, 0, 0);
        globalVar.gui_result.bar4.Visible = 'off';
        title('time interval');
        xlabel('time (second)') ;
        ylabel('time interval (ms)') ;
        axis tight;
    end
    
    % plot MSD curve, against tau (time interval)
    function plotMSD (axes) % input MSD in meter, and tau in second
        cla(axes);
        %nPlot = numel(cluster_selected.MSD_raw);
        set(axes, 'ColorOrder', makeMSDPlotColors(), 'NextPlot', 'add');
        set(axes, 'NextPlot', 'add');
        globalVar.gui_result.p5_1 = cell(nSeg, 1); 
        %p1 = [];
        globalVar.gui_result.p5_2 = [];
        %p3 = [];
        globalVar.gui_result.p5_3 = cell(nSeg, 1);
        
        for i = 1 : nSeg
            %p1{i} = scatter(axes, vertcat(cluster_selected.tau{:}), 1e12*vertcat(cluster_selected.MSD{:}), 1, [0 .35 .35]);
            globalVar.gui_result.p5_1{i} = plot(axes, 1e3*cluster_selected.tau_raw{i}, 1e12*cluster_selected.MSD_raw{i}, '.'); % - . --
            globalVar.gui_result.p5_3{i} = plot(axes, 1e3*cluster_selected.tau_reg{i}, 1e12*cluster_selected.MSD_reg{i}, 'LineWidth', 2); % plot MSD in unit um^2
            set(globalVar.gui_result.p5_3{i}, 'Visible', 'off');
        end
        %p1.CData = plotMSDColor_selected;
        %set(p1, 'Visible', 'on');
%         for i = 1 : nPlot
%             if (isempty(cluster_selected.MSD{i}))
%                 continue;
%             end
%             p1{i} = plot(axes, cluster_selected.tau{i}, 1e12*cluster_selected.MSD{i}, '.'); % plot MSD in unit um^2
%             p1{i}.Color(4) = 1/nthroot(nPlot, 4);
%         end
        if (~isempty(cluster_selected.MSD_mean))
            globalVar.gui_result.p5_2 = plot(axes, 1e3*cluster_selected.tau_mean, 1e12*cluster_selected.MSD_mean, 'LineWidth', 5, 'Color', [.2 .2 .2 .9]); % plot MSD in unit um^2
            set(globalVar.gui_result.p5_2, 'Visible', 'off');
        end
%         if (~isempty(cluster_selected.MSD_reg))
%             globalVar.gui_result.p5_3 = plot(axes, cluster_selected.tau_reg, 1e12*cluster_selected.MSD_reg, 'Color', [.1 .1 .1]); % plot MSD in unit um^2
%             set(globalVar.gui_result.p5_3, 'Visible', 'off');
%         end
        title('MSD per time interval');
        xlabel('tau (ms)');
        ylabel('MSD (um^2)');
        axis tight;

        %colormap(axes, plotMSDColor_selected);
        


    end

    % plot diffusion histogram
    function plotDiffusionHisto (axes) % input D in meter
        cla(axes);
        %set(axes, 'ColorOrder', plotMSDColor_selected, 'NextPlot', 'add');
        set(axes, 'NextPlot', 'add');
        %nPlot = numel(cluster_selected.DiffusionCoeff);
        %h1 = cell(nPlot, 1); 
        globalVar.gui_result.h6_1 = [];	globalVar.gui_result.h6_2 = [];	globalVar.gui_result.h6_3 = [];
        
        if (~isempty(cluster_selected.D_raw))
            globalVar.gui_result.h6_1 = histogram(axes, log10(1e12*vertcat(cluster_selected.D_raw{:})));%, 'FaceColor', [0 .35 .35], 'EdgeAlpha', 0);
        end

        if (~isempty(cluster_selected.D_reg))
            globalVar.gui_result.h6_3 = histogram(axes, log10(1e12*vertcat(cluster_selected.D_reg{:})));%, 'FaceColor', [0 .35 .35], 'EdgeAlpha', 0);
        end
        set(globalVar.gui_result.h6_3, 'Visible', 'off');

        title('Diffusion Coefficient histogram (log)');
        xlabel('log(D)');
        ylabel('counts');
        axis tight;
    end


    % plot diffusion cumulative frequency
    function plotDiffusionCumFreq (ax) % input D in meter
        cla(ax);
        %set(ax, 'ColorOrder', plotMSDColor_selected, 'NextPlot', 'add');
        set(ax, 'NextPlot', 'add');
        set(fig_result, 'CurrentAxes', ax); % save time comparing to axes()
        %axes(ax);

        %nPlot = numel(cluster_selected.DiffusionCoeff);
        %h1 = cell(nPlot, 1); 
        globalVar.gui_result.h7_1 = [];  globalVar.gui_result.h7_2 = []; globalVar.gui_result.h7_3 = [];
        
        if (~isempty(cluster_selected.D_raw))
            globalVar.gui_result.h7_1 = cdfplot(log10(1e12* vertcat(cluster_selected.D_raw{:})));
        end
        %plot(ax, h1);
        if (~isempty(cluster_selected.D_reg))
            globalVar.gui_result.h7_3 = cdfplot(log10(1e12* vertcat(cluster_selected.D_reg{:})));
        end
        set(globalVar.gui_result.h7_3, 'Visible', 'off');

        title('D cumulative frequency (log)');
        xlabel('log(D)') ;
        ylabel('cumulative frequency') ;
        axis tight;
    end

    % change overview color
     function colorOverview ()
        switch globalVar.colorMode
            case 0
                plotColor = globalVar.cluster_array_selected;
                titleTxt = 'colored by selected cluster';
                d = jet(256);
                colormap(globalVar.gui_overview.ax, d(75:205, :));
            case 1
                plotColor = globalVar.cluster_array;
                titleTxt = 'colored by cluster ID';
                colormap(globalVar.gui_overview.ax, globalVar.clusterColor);
            case 2
                plotColor = globalVar.trace_array;
                titleTxt = 'colored by trace ID';
                colormap(globalVar.gui_overview.ax, globalVar.traceColor);
        end
        globalVar.gui_overview.sc.CData = plotColor;

%         if (globalVar.colorMode == 0)
%             d = jet(256);
%             colormap(ax_overview, d(75:205, :));
%             %colormap(ax_overview, globalVar.traceColor);
%         elseif (globalVar.colorMode == 1)
%             colormap(ax_overview, globalVar.clusterColor);
%         else
%             colormap(ax_overview, globalVar.traceColor);
%         end
        title(globalVar.gui_overview.ax, titleTxt);
    end


%% event callback functions
    % load new particle analysis result data
    function loadResultData (~, ~)
        ise = evalin( 'base', 'exist(''file'',''var'') == 1' );
        if (ise)
            [last_dir, ~, ~] = fileparts(evalin('base','file'));
        else
            last_dir = '/home/';
        end
        [file, path] = uigetfile('*_result.mat', 'load analysis result file', last_dir); %, '_result.mat');
        if isequal(file, 0)
           return;
        end
        assignin('base', 'file', fullfile(path, file));
        result = load(fullfile(path, file), '-mat');
        particle_data = result.particle_data;
        assignin('base', 'particle_data', particle_data);
        loadData(particle_data);
    end

    % color change button on overview figure
    function changeOverviewColor (~, ~)
        globalVar.colorMode = rem(globalVar.colorMode + 1, 3);
        colorOverview();
    end

    % trace selection button on overview figure
    function BoxButtonPushed(~, ~)
        roi = drawrectangle(globalVar.gui_overview.ax, 'StripeColor', 'none', 'LineWidth', 0.5);
        tf = inROI(roi, globalVar.loc_all(:,1), globalVar.loc_all(:,2));
        ID = mode(globalVar.cluster_array(tf)); % get the mode cluster if many selected together
        %index_new = arrayfun(@(x) find(globalVar.cluster_ID_all == x), ID);
        index_new = find(globalVar.cluster_ID_all == ID);
        roi.delete();
        if (isempty(index_new)) 
            return;
        else
            index = index_new;
        end

        initSelectedCluster (index);

        colorOverview();

        plotCluster();
    end
    
    % call back function to update visible frame position
    function framePosUpdate(src, ~)%, time, xyz, ax1, ax2, ax3, ax4)
        globalVar.gui_result.vis_pos = round(src.Value);
        frameUpdate(globalVar.gui_result.vis_pos, globalVar.gui_result.vis_rng)%, time, xyz, ax1, ax2, ax3, ax4);
    end
    % call back function to update visible frame range
    function frameRangeUpdate(src, ~)%, time, xyz, ax1, ax2, ax3, ax4)
        rng = round(str2double(src.String));
        if isnan(rng)
            return;
        else
            globalVar.gui_result.vis_rng = rng;
        end
        frameUpdate(globalVar.gui_result.vis_pos, globalVar.gui_result.vis_rng)%, time, xyz, ax1, ax2, ax3, ax4);
    end
    % call back function to slider frame update
    function frameUpdate(pos, rng)%, time, xyz, ax1, ax2, ax3, ax4)
        
        time = cluster_selected.time;
        xyz = cluster_selected.coordinates;
       
        %range = numel(time{idx});

        % locate trace with longest time range
        N = cellfun(@(x) numel(x), time);
        [~, idx] = max(N);
        %xyz_trace = xyz{idx};
        xyz_trace = 1e9*(xyz{idx}-cluster_selected.center);
        
        N = size(xyz_trace, 1);
        t = time{idx} - min(time{idx});
        box_range = range(t)/N * rng;
        box_pos = t(pos) + box_range/2;
%        t = time{idx} - min(time{idx});
%         box_range = range(t)/N * rng;
%         box_pos = t(pos) + box_range/2;
        range_max = min(N, pos + rng);
        % update scatter plot curve
        points = fnplt( cscvn( xyz_trace(pos:range_max, :)' ) );
        delete(globalVar.gui_result.highlight_trace);
        delete(globalVar.gui_result.highlight_head);
        if (globalVar.nDim == 3) 
            globalVar.gui_result.highlight_trace = plot3(globalVar.gui_result.ax1, points(1,:), points(2,:), points(3,:), '-m');
            globalVar.gui_result.highlight_head = plot3(globalVar.gui_result.ax1, points(1, end), points(2, end), points(3, end), '-pentagram');
        elseif (globalVar.nDim == 2)
            globalVar.gui_result.highlight_trace = plot(globalVar.gui_result.ax1, points(1,:), points(2,:), '-m');
            globalVar.gui_result.highlight_head = plot3(globalVar.gui_result.ax1, points(1, end), points(2, end), '-pentagram');
        end
        globalVar.gui_result.highlight_trace.Visible = 'on';
        globalVar.gui_result.highlight_trace.LineWidth = 2;
        globalVar.gui_result.highlight_head.MarkerSize = 20;
        globalVar.gui_result.highlight_head.MarkerFaceColor = 'y';

        % update displacement, velocity, and time interval plot bar
        delete(globalVar.gui_result.bar2);
        globalVar.gui_result.bar2 = bar(globalVar.gui_result.ax2, ...
            box_pos, globalVar.gui_result.ax2.YLim(2), box_range, ...
            'BaseValue', globalVar.gui_result.ax2.YLim(1), 'FaceColor', 'm', 'FaceAlpha', .6, 'EdgeColor', 'none');
        delete(globalVar.gui_result.bar3);
        globalVar.gui_result.bar3 = bar(globalVar.gui_result.ax3, ...
            box_pos, globalVar.gui_result.ax3.YLim(2), box_range, ...
            'BaseValue', globalVar.gui_result.ax3.YLim(1), 'FaceColor', 'm', 'FaceAlpha', .6, 'EdgeColor', 'none');
        delete(globalVar.gui_result.bar4);
        globalVar.gui_result.bar4 = bar(globalVar.gui_result.ax4, ...
            box_pos, globalVar.gui_result.ax4.YLim(2), box_range, ...
            'BaseValue', globalVar.gui_result.ax4.YLim(1), 'FaceColor', 'm', 'FaceAlpha', .6, 'EdgeColor', 'none');     
    end
    

    function showRawMSD (src, ~)
        status = 'on';
        if (src.Value==0)
            status = 'off';
        end
        %set(globalVar.gui_result.p5_1, 'Visible', status);
        %set(globalVar.gui_result.h6_1, 'Visible', status);
        %set(globalVar.gui_result.h7_1, 'Visible', status);
        for i = 1 : numel(globalVar.gui_result.p5_1)
            set(globalVar.gui_result.p5_1{i}, 'Visible', status);
        end
        set(globalVar.gui_result.h6_1, 'Visible', status);
        set(globalVar.gui_result.h7_1, 'Visible', status);
        
    end
    function showMeanMSD (src, ~)
        status = 'on';
        if (src.Value==0)
            status = 'off';
        end
        set(globalVar.gui_result.p5_2, 'Visible', status);
        %set(globalVar.gui_result.h6_2, 'Visible', status);
        %set(globalVar.gui_result.h7_2, 'Visible', status);
    end
    function showRegMSD (src, ~)
        status = 'on';
        if (src.Value==0)
            status = 'off';
        end
        for i = 1 :numel(globalVar.gui_result.p5_3)
            set(globalVar.gui_result.p5_3{i}, 'Visible', status);
        end
        %set(globalVar.gui_result.p5_3, 'Visible', status);
        set(globalVar.gui_result.h6_3, 'Visible', status);
        set(globalVar.gui_result.h7_3, 'Visible', status);
    end


%% support functions
    % load new data
    function loadData (particle_data)
        initGlobalVar (particle_data);
        index = 1;
        initSelectedCluster (index);
        %% prepare overview and result figures
        % figure 901 as the overview figure
        if ~ishandle(901)
            fig_overview = figure(901);
            fig_overview.Position = [100 700 1000 500];
        end
        fig_overview = figure(901);
        fig_overview.Name = 'Localzation Data Overview';
        fig_overview.NumberTitle = 'off';
        globalVar.colorMode = 1;
        plotOverview();
        colorOverview();
        %  figure 902 as the measurement result figure for selected trace(s)
        if ~ishandle(902)
            fig_result = figure(902);
            fig_result.Position = [300 300 2000 900];
        end
        fig_result = figure(902);
        fig_result.NumberTitle = 'off';
        plotCluster(); %#ok<*SETNU> 
    end

    % get data dimension
    function nDim = getDimension (data)
        cor = data(1).coordinates{:};
        nDim = size(cor, 2);
    end

    % get selected trace IDs
    function trace_index_selected = getSelectedTraceIndex (traceIDs)
        trace_index_selected = zeros(globalVar.nLoc_all, 1);
        for i = 1 : numel(traceIDs)
            trace_index_selected(globalVar.cluster_array_all==traceIDs(i)) = 1;
        end
    end

    % compute RMS of 1D data
    function rms_data = RMS (data, windowSize)
        rms_data = sqrt( movmean(data.^2, windowSize) );
    end
    
    % get segments' trace belongings
    function locateSegments ()
        trace_t1 = cellfun(@(x) x(1), cluster_selected.time);
        trace_t2 = cellfun(@(x) x(end), cluster_selected.time);
        segment_t = cellfun(@(x) [x(1,1), x(end,1)], cluster_selected.segments, 'UniformOutput', false);
        trace_ID_segment = cellfun(@(x) find(x(1)>=trace_t1 & x(2)<=trace_t2), segment_t, 'UniformOutput', false);
        trace_ID_segment = vertcat(trace_ID_segment{:});
    end

    % generate glasbey LUT plot color
    function colors = generatePlotColors (nPlot)
        bg = [1 1 1]; % white background
        n_grid = 30;  % number of grid divisions along each axis in RGB space
        x = linspace(0, 1, n_grid);
        [R,G,B] = ndgrid(x,x,x);
        rgb = [R(:) G(:) B(:)];
        if (nPlot > size(rgb,1)/3)
            error('You can''t readily distinguish that many colors');
        end
        % Convert to Lab color space, which more closely represents human perception

        C = makecform('srgb2lab');
        lab = applycform(rgb, C);
        bglab = applycform(bg, C);


        mindist2 = inf(size(rgb,1),1);
        for i = 1:size(bglab,1)-1
            dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
            dist2 = sum(dX.^2,2);  % square distance
            mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
        end
        
        % Iteratively pick the color that maximizes the distance to the nearest
        % already-picked color
        colors = zeros(nPlot,3);
        lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
        for i = 1:nPlot
            dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
            dist2 = sum(dX.^2,2);  % square distance
            mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
            [~, id] = max(mindist2);  % find the entry farthest from all previously-chosen colors
            colors(i,:) = rgb(id,:);  % save for output
            lastlab = lab(id,:);  % prepare for next iteration
        end

    end
    
    % prepare multiple plot color order
    function plotColors = makeMSDPlotColors ()
        colors = generatePlotColors(nSeg);
        %alpha =  repmat([0.4; 1.0], nSeg, 1);
        plotColors = repelem(colors, 2, 1);
        %plotColors(:, 4) = alpha;
%         nMSD = cellfun(@(x) size(x, 1), cluster_selected.MSD);
%         trace_ID_segment_array = repelem(trace_ID_segment, nMSD, 1);
%         nSeg_trace = arrayfun(@(x) sum(trace_ID_segment_array==x), 1:globalVar.nTraces_selected)';
%         plotMSDColor_selected = repelem(globalVar.traceColor_selected, nSeg_trace, 1);
%         plotMSDColor_selected = [plotMSDColor_selected; [0 .55 .55]; [.1 .1 .1]];
    end

    function sf = plotSphereSurface (axes, center, radii, color, visibility)
        r = radii * ones(50, 50); % radius is 5
        [th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
        [x, y, z] = sph2cart(th, phi, r);
        x = x + center(1);  % center in x-direction
        y = y + center(2);  % center in y-direction
        z = z + center(3);  % center in z-direction
        sf = surface(axes, x, y, z, 'FaceColor', color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        sf.Visible = visibility;
    end

end
