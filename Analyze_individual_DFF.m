clc; clear; close all;
addpath('X:\OptoLab_v4.1\function\misc')

%PARAMETERS
%----------
%READ FILES (demodulated Calcium signals)
rFiles =['\**\*Ca.mat'];



opt.savePath = '';
opt.saveDo   = true;
opt.saveFuns = {... input for all: @(figure handle,savename)
    @(h,name)print(h,name,'-r300','-dpng');...
    ...@(h,name)print(h,name,'-r300','-dpdf');...
    };

%ANALYSIS OPTIONS
opt.offset = true;
%function arguments
opt.detrend.win = 2500; %[s] 
opt.detrend.degree = 2; %polynomal degree
opt.movmean.win = 5;  %movmean time [s] 

%stages in hypnogram
Stages = {... 
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    };
%stage transitions
Transitions = {... 
    'Wake','NREM';...
    'NREM','Wake';...
      'REM','Wake';...
    'NREM','REM';...
       };
margin  = 20; %[s], plots transtions +- margin 
exclNaN = true; %exclude transitions with margin beyond data range


%PLOT PROPERTIES
props.figure = {'visible','on'}; %off = invisible / on = visible
%axis
props.axi.hyp = {'box','off'}; %figure 1, hypnogram
%props.axi.cal = {'box','off'}; %figure 1, calcium signal
props.axi.cal = {'box','off','ylim',[0,1]}; %figure 1, calcium signal
props.axi.sta = {'box','off'}; %figure 1, mean across stages
props.axi.tra = {'box','off'}; %figure 2, transitions
props.equal.ylimTRA = false; %to set all transition axis equal ylim
%plot
props.plo.hyp = {};
props.plo.cal = {};
props.plo.sta = {}; %bar plot!
props.plo.tra = {'linewidth',2};
errorTrans  = 'SEM';

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
if ~opt.saveDo
    fprintf('[\bNOTE: data will not be saved (variable opt.saveDo)]\b\n')
end

%FILE-LIST
if iscell(rFiles)
    ind = cellfun(@(x)exist(x,'file')~=2,rFiles);
    if any(ind)
        fprintf('[\bFiles removed, %i of %i (do NOT exist)!]\b\n',...
            sum(ind),numel(ind))
        fprintf(' - %s\n',rFiles{ind})
        rFiles(ind) = [];
    end
    if isempty(rFiles)
        fprintf(2,'No File Left!\n')
        return
    end
elseif ischar(rFiles)
    tmp = dir(rFiles);
    tmp([tmp(:).isdir]) = [];
    tmp = fullfile({tmp.folder},{tmp.name});
    if isempty(tmp)
        fprintf(2,'No File Found For: %s\n',rFiles)
        return
    end
    rFiles = selectionList(tmp(:),'Ca*_DFF.mat');
    if isempty(rFiles)
        fprintf(2,'No File Selected\n')
        return
    end
else
    error('Class ''%s'' not supported for variable rFiles',class(rFiles))
end
noFIL = numel(rFiles);

%INIT
if ~iscell(opt.saveFuns)
    opt.saveFuns = {opt.saveFuns};
end
%number of ...
noSTA = size(Stages,1);
noTRA = size(Transitions,1);

%PLOT PROPERTIES
tmp = cell2mat(cellfun(@(x)x(:),Stages(:,1),'UniformOutput',false));
ticks = min(tmp):max(tmp);
label(size(ticks)) = {''};
for k = 1:noSTA
    [stageNUM,stageLAB] = Stages{k,:};
    label{ismember(ticks,stageNUM)} = stageLAB;
end
props.axi.hyp = ['ytick',ticks,'yticklabel',{label},props.axi.hyp];
%check (props must not be empty)
tmp = fieldnames(props.axi);
for k = 1:numel(tmp)
    if isempty(props.axi.(tmp{k}))
        props.axi.(tmp{k}) = {'visible','on'};
    end
end
tmp = fieldnames(props.plo);
for k = 1:numel(tmp)
    if isempty(props.plo.(tmp{k}))
        props.plo.(tmp{k}) = {'visible','on'};
    end
end

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    clear res
    [rPath,rFile,rExt] = fileparts(rFiles{fil});
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    %files
    fnames.cal = fullfile(rPath,[rFile,'.mat']);
    fnames.hyp = fullfile(rPath,'Hypnogram.mat');
    if ~exist(fnames.hyp,'file')
        [~,rFile,rExt] = fileparts(fnames.hyp);
        fprintf('%s [\b%s NOT found]\b\n',indent,[rFile,rExt])
        continue
    end
    tmp = {...
        sprintf('Created by %s.m, %s',scriptName,date);'';...
        'Resulting offset is already cut from results';...
        };
    res.info.info = char(tmp);
    
    %% LOAD CA-DATA
    fprintf('%s Load Data\n',indent)
    data   = load(fnames.cal);
    fs     = data.SampRate;
    signal = data.demodulated_signal(:);
    noSAM  = numel(signal);
    fprintf('%s   Ca-signal, fs = %g Hz\n',indent,fs)
    %id for title
    id = strrep(rFile,'_','\_'); %init
    if isfield(data,'info') && isfield(data.info,'file')
        [~,id] = fileparts(data.info.file); %from recorded filename
        id = strrep(id,'_','\_');
    end
    %get offset
    if isfield(data,'offset')
        offset = data.offset;
    else
        offset = [];
    end
    if opt.offset
        off = fun_CaOffset(signal,'fs',fs,'offset',offset,...
            'title',regexprep(fnames.cal,{'\','_'},{'\\\\','\\_'}));
        if ~isequal(off,offset)
            offset = off;
            save(fnames.cal,'offset','-append')
            fprintf('%s     Offset: %i (appended to calcium-file)\n',...
                indent,offset)
        else
            fprintf('%s     Offset: %i (old value kept)\n',...
                indent,offset)
        end
    else
        if isempty(offset)
            offset = 0;
            fprintf('%s     [\bOffset not set, uses 0]\b\n',indent)
        else
            fprintf('%s     Offset: %i (read from calcium-file)\n',...
                indent,offset)
        end
    end
    
    %% LOAD HYPNOGRAM
    data = load(fnames.hyp);
    hypnogram = data.Hypnogram(:);
    if isfield(data,'fs')
        fsHYP = data.fs;
        fprintf('%s   Hypnogram, fs = %g Hz\n',indent,fsHYP)
    else %estimate
        fsHYP = round(fs/noSAM*numel(hypnogram));
        fprintf('%s   Hypnogram, fs = %g Hz [\b(estimated!)]\b\n',...
            indent,fsHYP)     
    end
    % resample hypnogram
    if fs>fsHYP %upsample
        fac = fs/fsHYP;
        if fac==round(fac)
            hypnogram = repmat(hypnogram',fac,1);
            hypnogram = hypnogram(:);
            fprintf('%s     Upsampled by factor %i\n',indent,fac)
        else
            tmp = round(resample(hypnogram,fs,fsHYP,0));
            tmp = tmp(:);
            %time shift correction 
            shift = round((1/fsHYP-1/fs)*fs); % >= 0 !!!
            if shift>0
                tmp = [repmat(tmp(1),shift,1);tmp(1:end-shift)];
            end
            hypnogram = tmp(:);
            fprintf('%s     Upsampled to %i Hz\n',indent,fs)
        end
    elseif fs<fsHYP %downsample
        
        tmp = resample(hypnogram,fs,fsHYP,0);
        tmp = tmp(:);
        %time shift correction 
        shift = round((1/fs-1/fsHYP)*fs);
        if shift>0
            tmp = [tmp(shift+1:end);repmat(tmp(end),shift,1)];
        end
        if false %test plot
            figure
            t1 = (1:numel(hypnogram))/fsHYP;
            t2 = (1:numel(tmp))/fs;
            plot(t1,hypnogram,'b','linewidth',2); hold on
            plot(t2,tmp,':r','linewidth',2);
            mima = [min(hypnogram),max(hypnogram)];
            set(gca,'ylim',mima+0.1*[-1,1]*diff(mima),'ytick',...
                mima(1):mima(2),'xlim',[0,max([t1(end),t2(end)])])
            legend('hypnogram','downsampled')
            zoom xon
            return
        end
        hypnogram = tmp(:);
        fprintf('%s     Downsampled to %i Hz\n',indent,fs)
    end
    %same size signal & hypnogram 
    if numel(hypnogram)<noSAM
        hypnogram(end+1:noSAM) = NaN;
    elseif numel(hypnogram)>noSAM
        hypnogram(noSAM+1:end) = [];
    end
    
    %CUT OFFSET
    res.info.offset = offset;
    offSIG = []; offHYP = []; %init
    if offset>0
        ind = 1:offset;
        offSIG = signal(ind);
        offHYP = hypnogram(ind);
        signal(ind)    = [];
        hypnogram(ind) = [];
        noSAM = numel(signal);
    end
    
    %% DFF
    %trend
    fprintf('%s Calculate DFF\n',indent)
    %detrend signal
    win = opt.detrend.win*fs;
    trend = signal - detrend(signal,opt.detrend.degree,... 
        1:win:noSAM,... 
        'Continuous',false);
    trend = smooth(trend,win/2);
    signalD = signal-trend;
    %smooth signal
    bins = round(opt.movmean.win*fs);
    if mod(bins,2)==0
        bins = bins+1;
    end
    signalM = movmean(signalD,bins); %smooth signal
    mi = prctile(signalD,1);
    ma = prctile(signalD,99);
    %dff (like dff)
    res.fs  = fs;
    res.dff = (signalM-mi)/(ma-mi);
    res.hypnogram = hypnogram;
    
    %% MEAN ACROSS STAGES
    res.stages.label = cell(1,noSTA);
    res.stages.mean  = NaN(1,noSTA);
    res.stages.std   = NaN(1,noSTA);
    res.stages.N     = NaN(1,noSTA);
    indCUT = true(size(hypnogram));
    if opt.cutAfterStage4
        indCUT(find(hypnogram==4,1,'last')+1:end) = false;
    end
    for sta = 1:noSTA
        [stage,label] = Stages{sta,:};
        ind = ismember(hypnogram,stage) & indCUT;
        res.stages.label{sta} = label;
        res.stages.mean(sta)  = mean(res.dff(ind));
        res.stages.std(sta)   = std(res.dff(ind));
        res.stages.N(sta)     = sum(ind);
    end
    
    %% MEAN ACROSS TRANSITIONS
    ind0 = -round(margin*fs):round(margin*fs);
    t = ind0'/fs;
    res.trans.t = t;
    res.trans.label = cell(1,noSTA);
    res.trans.mean  = NaN(numel(t),noSTA);
    res.trans.std   = NaN(numel(t),noSTA);
    res.trans.N     = NaN(1,noSTA);
    for tra = 1:noTRA
        [stage1,stage2] = Transitions{tra,:};
        nums1 = Stages{strcmpi(Stages(:,2),stage1),1};
        nums2 = Stages{strcmpi(Stages(:,2),stage2),1};
        indTRA = find(...
            ismember(hypnogram(1:end-1),nums1) & ...
            ismember(hypnogram(2:end)  ,nums2));
        %matrix of transitions
        N = numel(indTRA);
        M = NaN(numel(ind0),N);
        for k = 1:N
            ind = ind0+indTRA(k); 
            ind1 = find(ind>=1,1,'first');    
            ind2 = find(ind<=noSAM,1,'last');
            M(ind1:ind2,k) = res.dff(ind(ind1):ind(ind2));
        end
        %exclude transitions with margin beyond data range
        if exclNaN
            ind = any(isnan(M),1);
            M(:,ind) = [];
            N = size(M,2);
        end
        %res
        res.trans.label{tra}  = sprintf('%s - %s',stage1,stage2);
        res.trans.mean(:,tra) = nanmean(M,2);
        res.trans.std(:,tra)  = nanstd(M,[],2);
        res.trans.N(tra)      = N;
        
      
        
        
    end
    
    
    %% PLOT INIT
    fprintf('%s Plots\n',indent)
    clear hf;
    
    %% FIGURE 
    clear hf;
    hf(1) = figure(props.figure{:}); figLab{1} = 'procedure';
    ha = NaN(4,1);
    t  = (1:noSAM)'/fs/3600;
    to = (-offset+1:0)/fs/3600;
    xl = [min([to(:);0]),t(end)];
    %subplot, signal. filter, hugline
    ha(1) = subplot(311); hold on
    legSTR = {}; %init
    if offset>0
        plot(to,offSIG,'k');
        legSTR(end+1) = {'Offset'};
    end
    plot(t,signal,'b');
    legSTR(end+1) = {'Raw Signal'};
    plot(t,trend,'c');
    legSTR(end+1) = {'Trend'};
    if max(signal)<10^-5 %probably noise and not real signal
        title({sprintf('Signal %s, probably noise!',id),...
            sprintf('Offset = %i',offset)},'color','r')
    else
        title({sprintf('Signal %s',id),sprintf('Offset = %i',offset)})
    end
    legend(legSTR);
    tmp = [offSIG;signal]; ind = ceil(noSAM/5*4):noSAM;
    yl = [min(tmp),max(tmp)]; yl(2) = max([2*tmp(ind)-yl(1);yl(2)]);
    set(ha(1),'ylim',yl+0.05*[-1,1]*diff(yl),'box','off')
    %subplot, detrended
    ha(2) = subplot(312); hold on
    plot(t,signalD,'b');
    plot(t,signalM,'c');
    legend('Detrended','Det. & Smoothed')
    tmp = signalD; ind = ceil(noSAM/5*4):noSAM;
    yl = [min(tmp),max(tmp)]; yl(2) = max([2*tmp(ind)-yl(1);yl(2)]);
    set(ha(2),'ylim',yl+0.05*[-1,1]*diff(yl),'box','off')
    title('Detrended and Smoothed (movmean)');
    %subplot, hypnogram & signal
    ha(3) = subplot(313);
    plot(t,hypnogram,'color',ones(1,3)*0.55)
    yl = [min(hypnogram),max(hypnogram)];
    if numel(unique(yl))==1
        yl = [-.5,.5]+unique(yl);
    end
    set(ha(3),props.axi.hyp{:},'YAxisLocation','right',...
        'ylim',yl+0.1*[-1,1]*diff(yl),'xticklabel',[])
    ha(4) = axes('position',get(ha(3),'position'),'xticklabel',[]);
    plot(t,res.dff,'b');
    yl = [min(res.dff),max(res.dff)];
    set(ha(4),'ylim',yl+[-0.05,0.2]*diff(yl),props.axi.cal{:},'color','none')
    title('\DeltaF/F, hypnogram');
    %settings
    linkaxes(ha,'x')
    set(ha,'xlim',xl)
    xlabel('Time [h]')
    
    
    %% SAVE
    drawnow
    if opt.saveDo
        %save path
        if ischar(opt.savePath) && exist(opt.savePath,'dir')
            sPath = opt.savePath;
            uniqueFiles = true;
        else
            sPath = rPath;
            uniqueFiles = false;
        end
        %check
        noFIG = numel(hf);
        sFile = sprintf('%s_DFF',rFile);
        if numel(unique(lower(figLab)))~=noFIG
            error('Figure labels must be unique!')
        end
        %savenames
        snames = cell(noFIG+1,1);
        if ~uniqueFiles
            snames{1} = fullfile(sPath,sprintf('%s.mat',sFile));
            for k = 1:noFIG 
                snames{k+1} = fullfile(sPath,...
                    sprintf('%s_%s',sFile,figLab{k}));
            end
        else
            cnt = 0; tmp = true; %init
            while tmp
                cnt = cnt+1;
                snames{1} = fullfile(sPath,sprintf('%s_%i.mat',...
                    sFile,cnt));
                for k = 1:noFIG 
                    snames{k+1} = fullfile(sPath,sprintf('%s_%s_%i',...
                        sFile,figLab{k},cnt));
                end
                tmp = exist(snames{1},'file')==2;
            end
        end
        %save data
        sname = snames{1};
        save(sname,'-struct','res')
        [~,sFile,sExt] = fileparts(sname);
        fprintf('%s Saved Data  : %s\n',indent,[sFile,sExt])
        %save figures
        for fig = 1:noFIG
            sname = snames{fig+1};
            for k = 1:numel(opt.saveFuns)
                fun = opt.saveFuns{k};
                fun(hf(fig),sname)
            end
            [~,sFile] = fileparts(sname);
            fprintf('%s Saved Figure: %s\n',indent,sFile)
        end
        if fil<noFIL
            close(hf)
        end
    end
end