
%PARAMETERS
%----------
%DRIVE LETTER 
driveLetter = 'X:'; %inclusive double point !!!
addpath(genpath(fullfile(driveLetter,'OptoLab_v4.1\function')))

%DATA LIST (excel file with data info)
dataList = fullfile(driveLetter,...
    ['.xlsx']);
%column labels 
Columns = { ... {variable name, label in xls-file}
    'path'   ,'Path';...
    'animal' ,'Mouse ID';...
    'area'   ,'Target Area';...
    'sensor' ,'Sensor';...
    'trial'  ,'Trial';...
    };
%filename of calcium file
fileCA = 'Ca_DFF.mat';


fs = 128; %[Hz]

comp.error = 'SEM';

comp.axes = {... comparison within axis
   ... {variable, {element list}}
   'area',{'BF','LH'};...
   };
comp.cmp = {... pre-selection of data, plots one axis per row
   ... {variable, {element list}}
   'sensor',{'OXLight'};...
   'sensor',{'OXLight+Suvorexant'};...
   };

sel.animal = {};
sel.trial  = {}; 
sel.area   = {};
sel.sensor = {};

%SAVE OPTIONS
opt.save.figs = true;
opt.save.dat  = true; 
sPath = fullfile(driveLetter,... 
    ['']);
opt.save.funsFig = {... figures, several possible
    @(h,name,str)print(h,name,'-dpng','-r300');...
    @(h,name,str)print(h,name,'-depsc','-r300');...
    };

%PROPERTES (must not be empty, use 'visible','on' if nothing else)
props.fig = {'renderer','painter'}; 
props.axi.stages = {'box','off'};
props.axi.trans1 = {'visible','on'}; 
props.axi.trans2 = {'visible','on'}; 

%% MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));

if ~iscell(opt.save.funsFig) 
    opt.save.funsFig = {opt.save.funsFig};
end
if numel(driveLetter)~=2 || driveLetter(2)~=':'
    error(sprintf(['driveLetter must have 2 characters inclusive ',...
        'double point.\nE.g. driveLetter = ''C:''']))
end

%PREPARE DATA LIST
[~,~,DATA] = xlsread(dataList);
ind = cellfun(@ischar,DATA); 
DATA(ind) = cellfun(@strtrim,DATA(ind),'UniformOutput',false);
Columns   = cellfun(@strtrim,Columns,'UniformOutput',false);
ind = cellfun(@isnumeric,DATA);
DATA(ind) = cellfun(@num2str,DATA(ind),'UniformOutput',false);
%split label, Data
labelsXLS = DATA(1,:);
DATA(1,:) = [];
%get column number
for k = 1:size(Columns)
    [field,label] = Columns{k,:};
    col = find(strcmpi(labelsXLS,label));
    n   = numel(col);
    switch n
        case 0
            error('Label ''%s'' not found in DATA',label)
        case 1
            columns.(field) = col;
        otherwise
            error('Label ''%s'' must be unique (found N = %i)',label,n)
    end
end
%remove non char, non paths
ind = cellfun(@(x) ischar(x) && numel(x)>=2 && x(2)==':',...
    DATA(:,columns.path));
DATA(~ind,:) = [];
%change to correct drive letter
DATA(:,columns.path) = cellfun(@(x)[driveLetter,x(3:end)],...
    DATA(:,columns.path),'uniformoutput',false);
%remove in-existing paths
ind = cellfun(@(x)exist(x,'dir'),DATA(:,columns.path));
DATA(~ind,:) = [];
%sort data
[~,ind] = natsort(DATA(:,columns.path));
DATA = DATA(ind,:);
if isempty(DATA)
    fprintf(2,'No Valid Data Found!\n')
    return
end

%SELECT DATA
%pre-selection
tmp = fieldnames(sel);
for k = 1:numel(tmp)
    label = tmp{k};
    dat   = sel.(label);
    if isempty(dat)
        continue
    end
    if ~all(cellfun(@ischar,dat))
        error('Data selection sel.%s must be cell of chars!',label)
    end
    
    DATA(~ismember(lower(DATA(:,columns.(label))),lower(dat)),:) = [];
end
if isempty(DATA)
    fprintf(2,'No Data Left After Pre-Selection!\n')
    return
end
%select by paths
clear tmp;
tmp.sort = false; tmp.unique = false;
[~,ind] = selectionList(DATA(:,columns.path),[],tmp);
DATA = DATA(ind,:); clear tmp;
if isempty(DATA)
    fprintf(2,'No Data Path Selected!\n')
    return
end

%ANIMAL/AREAS/SENSORS/TRIALS
animals = natsort(unique(DATA(:,columns.animal)));
areas   = natsort(unique(DATA(:,columns.area)));
sensors = natsort(unique(DATA(:,columns.sensor)));
trials  = natsort(unique(DATA(:,columns.trial)));
%number of ...
noANI = numel(animals);
noSEN = numel(sensors);
noARE = numel(areas);
noTRI = numel(trials);
%print out
fprintf('Available Data of Selected Paths:\n')
n = numel(num2str(max([noANI,noSEN,noARE,noTRI])));
fprintf('  Animals (N = %*i): %s\n',n,noANI,strjoin(animals,', '))
fprintf('  Areas   (N = %*i): %s\n',n,noARE,strjoin(areas,', '))
fprintf('  Sensors (N = %*i): %s\n',n,noSEN,strjoin(sensors,', '))
fprintf('  Trials  (N = %*i): %s\n',n,noTRI,strjoin(trials,', '))

%COMPARISONS & CHECK
[noAXI,tmp] = size(comp.axes);
if tmp~=2
    error('Variable comp.axes must have two columns')
end
[noCMP,tmp] = size(comp.cmp);
if tmp~=2
    error('Variable comp.cmp must have two columns')
end
%fill empty cells
for k = 1:noAXI
    [label,data] = comp.axes{k,:};
    if ~iscell(data)
        error('Entries in comp.axes(:,2) must be cells')
    end
    if isempty(data)
        comp.axes{k,2} = natsort(unique(DATA(:,columns.(label))));
    end
end
for k = 1:noCMP
    [label,data] = comp.cmp{k,:};
    if ~iscell(data)
        error('Entries in comp.cmp(:,2) must be cells')
    end
    if isempty(data)
        comp.axes{k,2} = natsort(unique(DATA(:,columns.(label))));
    end
end
%colors compare
colsCMP = lines(noCMP);

%% READ DATA (per axis and comparison)
%LOOPES
fprintf('\nLOAD DATA\n')
columnAXI = columns.(comp.axes{1});
columnCMP = columns.(comp.cmp{1});
DataSTA = cell(noAXI,noCMP);
DataTRA = cell(noAXI,noCMP);
clear stages 
%axis loop
for axi = 1:noAXI
    [label,vars] = comp.axes{axi,:};
    col = columns.(lower(label));
    ind = ismember(lower(DATA(:,col)), lower(vars));
    DAT = DATA(ind,:);
    fprintf(' Axis: %s\n',strjoin(unique(DAT(:,col))))
    
    %compare loop
    for cmp = 1:noCMP
        [label,vars] = comp.cmp{cmp,:};
        col = columns.(lower(label));
        ind = ismember(lower(DAT(:,col)), lower(vars));
        D = DAT(ind,:);
        noFIL = size(D,1);
        fprintf('   Comp: %s (N = %i)\n',strjoin(unique(D(:,col))),noFIL)
        
        %file loop
        datSTA = []; datTRA = [];
        for fil = 1:noFIL
            rPath = D{fil,columns.path};
            
            %LOAD DATA
            rname = fullfile(rPath,fileCA);
            if exist(rname,'file')~=2
                fprintf(2,'   File NOT found: %s\n',[rPath,fileCA])
                
                continue
            end
            data = load(rname);
            
            %INIT DATA
            if ~exist('stages','var')
                stages       = data.stages.label;
                transitions  = data.trans.label;
                noSTA = numel(stages);
                noTRA = numel(transitions);
                t = data.trans.t(:)'; %row vector
            elseif ~isequal(stages,data.stages.label) || ...
                    ~isequal(transitions,data.trans.label) || ...
                    ~isequal(t,data.trans.t(:)')
                error('Stages in different files must be equal')
            end
            if isempty(datSTA)
                datSTA = NaN(noFIL,noSTA);
                datTRA = NaN(noFIL,noTRA,numel(t));
            end
            
            %APPEND
            datSTA(fil,:)   = data.stages.mean;
            datTRA(fil,:,:) = data.trans.mean';
        end %file loop
        DataSTA{axi,cmp} = datSTA;
        DataTRA{axi,cmp} = datTRA;
    end %compare loop
end %axis loop


%% AVERAGE DATA
tmp = cell(size(DataSTA));
mDat.stage.data  = tmp;
mDat.stage.mean  = tmp;
mDat.stage.std   = tmp;
mDat.stage.n     = tmp;
mDat.trans.mean  = tmp;
mDat.trans.std   = tmp;
mDat.trans.n     = tmp;
mDat.trans.pks.y = NaN(noAXI,noCMP,noTRA,2);
mDat.trans.pks.t = NaN(noAXI,noCMP,noTRA,2);
indL = t>=-15 & t<=0; %left, prior transition
indR = t>=0  & t<=15; %right, after transition
for axi = 1:noAXI
    for cmp = 1:noCMP
        %stage data
        data = DataSTA{axi,cmp};
        mDat.stage.data{axi,cmp} = data;
        mDat.stage.mean{axi,cmp} = nanmean(data,1);
        mDat.stage.std{axi,cmp}  = nanstd(data,[],1);
        mDat.stage.n{axi,cmp}    = sum(~isnan(data),1);
        %transition data
        data = shiftdim(DataTRA{axi,cmp},1);
        Mea = nanmean(data,3);
        mDat.trans.mean{axi,cmp} = Mea;
        mDat.trans.std{axi,cmp}  = nanstd(data,[],3);
        mDat.trans.n{axi,cmp}    = sum(~isnan(data),3);
        %peaks
        for tra = 1:noTRA
            mea =  Mea(tra,:);
            yPKS = NaN(1,2);
            indP = NaN(1,2);
            clear tmp
            if max(mea(indL))>max(mea(indR))
                tmp.L = mea; tmp.L(~indL) = -inf;
                tmp.R = mea; tmp.R(~indR) = inf;
                [yPKS(1),indP(1)] = max(tmp.L);
                [yPKS(2),indP(2)] = min(tmp.R);
            else
                tmp.L = mea; tmp.L(~indL) = inf;
                tmp.R = mea; tmp.R(~indR) = -inf;
                [yPKS(1),indP(1)] = min(tmp.L);
                [yPKS(2),indP(2)] = max(tmp.R);
            end
            tPKS = t(indP);
            mDat.trans.pks.y(axi,cmp,tra,:) = yPKS;
            mDat.trans.pks.t(axi,cmp,tra,:) = tPKS;
        end %loop transition
    end %looop cmp
end %loop axis
%additional (new) for plot 0
tmp = cell(size(DataSTA));
mDat.stageDiff.labels = tmp;
mDat.stageDiff.data   = tmp;
mDat.stageDiff.mean   = tmp;
mDat.stageDiff.std    = tmp;
mDat.stageDiff.n      = tmp;
for axi = 1:noAXI
    for cmp = 1:noCMP
        data = DataSTA{axi,cmp};
        %stage differences as index
        IND = [...
            2,1;...
            2,3;...
            1,3;...
            ];
        nn   = size(IND,1);
        dat  = NaN(size(data));
        mea  = NaN(1,nn);
        err  = NaN(1,nn);
        n    = NaN(1,nn);
        labs = cell(1,nn);
        for k = 1:nn
            ind = IND(k,:);
            %tmp = data(:,ind(1))-data(:,ind(2));
            tmp = abs(data(:,ind(1))-data(:,ind(2)));
            dat(:,k) = tmp;
            mea(k)   = nanmean(tmp);
            err(k)   = nanstd(tmp);
            n(k)     = sum(~isnan(tmp));
            labs{k}  = sprintf('%s - %s',stages{ind});
        end
        mDat.stageDiff.labels{axi,cmp} = labs;
        mDat.stageDiff.data{axi,cmp}   = dat;
        mDat.stageDiff.mean{axi,cmp}   = mea;
        mDat.stageDiff.std{axi,cmp}    = err;
        mDat.stageDiff.n{axi,cmp}      = n;
    end
end

%% PLOT I 
%clc; close all
fprintf('\nPLOT FIGURE I (stages)\n')
%figure % axes
dx = [55,50,20]; dy = [50,50,50];
rows = noAXI; cols = 1;
wid = 500; hig = 150; %axes width/height
hf = figure('visible','off','unit','pixel','position',...
    [200,200,cols*wid+[1,cols-1,1]*dx',rows*hig+[1,rows-1,1]*dy'],...
    props.fig{:});
movegui(hf,'center'); drawnow
set(hf,'visible','on');
ha = fig_createAxes(hf,[rows,cols],dx,dy,'pixel');
set(ha,'unit','normalized')
%init
noSTAd = numel(mDat.stageDiff.labels{1});
x0 = 1:noSTAd;
dx = linspace(-.5,.5,noCMP+2);
dl = mean(diff(dx));
dx([1,end]) = [];
maxY = -inf; minY = inf;
%axie loop
for axi = 1:noAXI
    labelAXI = strjoin(comp.axes{axi,2},', ');
    set(hf,'CurrentAxes',ha(axi))
    %compare loop
    legSTR = cell(1,noCMP);
    hp = NaN(1,noCMP);
    for cmp = 1:noCMP
        labelCMP = strjoin(comp.cmp{cmp,2},', ');
        %data
        labs = mDat.stageDiff.labels{axi,cmp};
        dat  = mDat.stageDiff.data{axi,cmp};
        mea  = mDat.stageDiff.mean{axi,cmp};
        err  = mDat.stageDiff.std{axi,cmp};
        n    = max(mDat.stageDiff.n{axi,cmp});
        switch lower(comp.error)
            case 'std'
            case 'sem'
                err = err/sqrt(n);
            otherwise
                error('uups')
        end
        legSTR{cmp} = sprintf('%s (n = %i)',labelCMP,n);
        x = x0+dx(cmp);
        maxY = max([maxY;mea(:)+err(:)]);
        minY = min([minY;mea(:)-err(:)]);
        %plot
        errorbar(x,mea,err,err ,'linestyle','none','color','k');
        hold on
        hp(cmp) = bar(x,mea,'barwidth',0.8*dl,'facecolor',colsCMP(cmp,:));
        for k = 1:numel(x)
            tmp = dat(:,k);
            h = plot(repmat(x(k),size(tmp)),tmp,'o',...
                'markersize',4,'linewidth',0.5,'color',[1,1,1]*0);
            %plot(x(k),nanmean(tmp),'ko') 
            minY = min([minY;min(tmp)]);
            maxY = max([maxY;max(tmp)]);
        end
    end %compare loop
    %text
    title(labelAXI)
    ylabel(sprintf('Mean + %s',comp.error))
    legend(hp,legSTR)
end %axie loop
%settings
yl = NaN(1,2);
if minY<0
    yl(1)=1.1*minY;
else
    yl(1)=0;
end
if maxY>0
    yl(2)=1.1*maxY;
else
    yl(2)=0;
end
set(ha,'xlim',[0.3,noSTAd+0.7],'xtick',x0,'xticklabel',labs,...
    'ylim',yl,props.axi.stages{:})
linkaxes(ha,'xy')

%SAVE
if opt.save.figs
    sname = fullfile(sPath,'DFF_diffStages');
    for k = 1:numel(opt.save.funsFig)
        fun = opt.save.funsFig{k};
        fun(hf,sname)
    end
    fprintf(' - Saved: %s.*\n',sname)
    %save individual data
    DD = mDat.stageDiff;
    tab = {'axis','data','stages'}';
    for axi = 1:noAXI
        labelAXI = strjoin(comp.axes{axi,2},', ');
        for cmp = 1:noCMP
            labelCMP = strjoin(comp.cmp{cmp,2},', ');
            %append
            D = DD.data{axi,cmp};
            [rows,cols] = size(D);
            rows = (1:rows)+3;
            cols = (1:cols)+size(tab,2);
            tab{1,cols(1)} = labelAXI;
            tab{2,cols(1)} = labelCMP;
            tab(3,cols)    = DD.labels{axi,cmp};
            tab(rows,cols) = num2cell(D);
        end
    end
    sname = [sname,'.xlsx'];
    if exist(sname,'file')==2
        delete(sname)
    end
    xlswrite(sname,tab)
    try
        xls_cellFit(sname)
    catch
    end
end


%% PLOT II (stages)
%clc; close all
fprintf('PLOT FIGURE II (stages)\n')
%figure % axes
dx = [55,50,20]; dy = [50,50,50];
rows = noAXI; cols = 1;
wid = 500; hig = 150; %axes width/height
hf = figure('visible','off','unit','pixel','position',...
    [200,200,cols*wid+[1,cols-1,1]*dx',rows*hig+[1,rows-1,1]*dy'],...
    props.fig{:});
movegui(hf,'center'); drawnow
set(hf,'visible','on');
ha = fig_createAxes(hf,[rows,cols],dx,dy,'pixel');
set(ha,'unit','normalized')
%init
x0 = 1:noSTA;
dx = linspace(-.5,.5,noCMP+2);
dl = mean(diff(dx));
dx([1,end]) = [];
maxY = 0;
%axie loop
for axi = 1:noAXI
    labelAXI = strjoin(comp.axes{axi,2},', ');
    set(hf,'CurrentAxes',ha(axi))
    %compare loop
    legSTR = cell(1,noCMP);
    for cmp = 1:noCMP
        labelCMP = strjoin(comp.cmp{cmp,2},', ');
        %data
        dat = mDat.stage.data{axi,cmp};
        mea = mDat.stage.mean{axi,cmp};
        err = mDat.stage.std{axi,cmp};
        n   = mDat.stage.n{axi,cmp};
        if strcmpi(comp.error,'SEM')
            err = err./sqrt(n);
        end
        legSTR{cmp} = sprintf('%s (n = %i)',labelCMP,max(n));
        x = x0+dx(cmp);
        maxY = max([maxY;mea(:)+err(:)]);
        %plot
        errorbar(x,mea,NaN(size(err)),err ,'linestyle','none','color','k');
        hold on
        hp(cmp) = bar(x,mea,'barwidth',0.8*dl,'facecolor',colsCMP(cmp,:));
        for k = 1:numel(x)
            tmp = dat(:,k);
            h = plot(repmat(x(k),size(tmp)),tmp,'o',...
                'markersize',4,'linewidth',0.5,'color',[1,1,1]*0);
            %plot(x(k),nanmean(tmp),'ko') 
            maxY = max([maxY;max(tmp)]);
        end
    end %compare loop
    %text
    title(labelAXI)
    ylabel(sprintf('Mean + %s',comp.error))
    legend(hp,legSTR)
end %axie loop
% %settings
set(ha,'xlim',[0.3,noSTA+0.7],'xtick',x0,'xticklabel',stages,...
    'ylim',[0,maxY*1.1],props.axi.stages{:})
linkaxes(ha,'xy')

%SAVE
if opt.save.figs
    sname = fullfile(sPath,'DFF_stages');
    for k = 1:numel(opt.save.funsFig)
        fun = opt.save.funsFig{k};
        fun(hf,sname)
    end
    fprintf(' - Saved: %s.*\n',sname)
    %save individual data
    DD = mDat.stage;
    tab = {'axis','data','stages'}';
    for axi = 1:noAXI
        labelAXI = strjoin(comp.axes{axi,2},', ');
        for cmp = 1:noCMP
            labelCMP = strjoin(comp.cmp{cmp,2},', ');
            %append
            D = DD.data{axi,cmp};
            [rows,cols] = size(D);
            rows = (1:rows)+3;
            cols = (1:cols)+size(tab,2);
            tab{1,cols(1)} = labelAXI;
            tab{2,cols(1)} = labelCMP;
            tab(3,cols)    = stages;
            tab(rows,cols) = num2cell(D);
        end
    end
    sname = [sname,'.xlsx'];
    if exist(sname,'file')==2
        delete(sname)
    end
    xlswrite(sname,tab)
    try
        xls_cellFit(sname)
    catch
    end
end
%% PLOT III (transition)
% clc; close all
% And min, max peak data
fprintf('PLOT FIGURE III (transitions_20)\n')
%figure/axes options
dx = [55,65,20]; dy = [45,65,50];
rows = ceil(sqrt(noTRA)); cols = ceil(noTRA/rows);
wid  = 200; hig = 200; %axes width/height
pos  = [200,200,cols*wid+[1,cols-1,1]*dx',rows*hig+[1,rows-1,1]*dy'];
%FIGURE LOOP (based on axes)
tabPKS = {'Peaks','','1st','2nd'};
printFiltered = true;
for axi = 1:noAXI
    labelAXI = strjoin(comp.axes{axi,2},', ');
    
    %FIGURE AXES
    hf = figure('visible','off','unit','pixel','position',pos,props.fig{:});
    movegui(hf,'center'); drawnow
    set(hf,'visible','on');
    ha = fig_createAxes(hf,[rows,cols],dx,dy,'pixel')';
    set(ha,'unit','normalized')
    MIMA = NaN(size(ha,1),size(ha,2),2); %init
    
    %CMP LOOP
    hp = NaN(1,noCMP);
    legSTR = cell(1,noCMP);
    for cmp = 1:noCMP
        labelCMP = strjoin(comp.cmp{cmp,2},', ');
        Mea = mDat.trans.mean{axi,cmp};
        Err = mDat.trans.std{axi,cmp};
        N   = mDat.trans.n{axi,cmp};
        if strcmpi(comp.error,'SEM')
            Err = Err./sqrt(N);
        end
        
        %TRANSITION LOOP
        maxN = 0;
        for tra = 1:noTRA
            trans = transitions{tra};
            set(hf,'CurrentAxes',ha(tra))
            %data
            n    = N(tra,:);
            mea  = Mea(tra,:);
            err  = Err(tra,:);
            xPks = squeeze(mDat.trans.pks.t(axi,cmp,tra,:));
            yPks = squeeze(mDat.trans.pks.y(axi,cmp,tra,:));
            
            %RESAMPLE (filter and interpolation)
            fsR = round(1/mean(diff(t))); %sampling rate raw data
            if ~isempty(fs) && fs<fsR
                if printFiltered
                    fprintf([' - Data filtered/interpolated, ',...
                        '%g --> %g Hz\n'],fsR,fs)
                    printFiltered = false;
                end
                t2 = t(1):1/fs:t(end);
                [b,a] = butter(6,fs/3/fsR);
                mea   = interp1(t,filtfilt(b,a,mea),t2);
                err   = interp1(t,filtfilt(b,a,err),t2);
            else
                t2 = t;
            end
            
            %PLOT
            [row,col] = find(ha==gca);
            %error area
            if true
                x = [t2,t2(end:-1:1)];
                y = [mea-err,mea(end:-1:1)+err(end:-1:1)];
                hh = fill(x,y,colsCMP(cmp,:),'edgecolor','none',...
                    'facealpha',0.3);
                uistack(hh,'bottom')
                hold on
                
                MIMA(row,col,:) = [min(mea-err),max(mea+err)];
            else
                err = NaN;
                MIMA(row,col,:) = [min(mea),max(mea)];
            end
            %mean
            hp(cmp) = plot(t2,mea,'color',colsCMP(cmp,:),'linewidth',2);
            hold on
            %min max
            for k = 1:2
                plot(xPks(k),yPks(k),'o','color',colsCMP(cmp,:),...
                    'linewidth',2);
            end
            % append
            [~,ind] = sort(xPks(:));
            tmp = num2cell(yPks(ind));
            tabPKS(end+1,:) = [trans,labelCMP,tmp'];
            %text
            if cmp==1
                title({trans})
                xlabel('Time [s]')
                ylabel(sprintf('Mean \\pm %s',comp.error))
            end
            maxN = max([maxN,max(n)]);
        end %loop cmp
        %text
        legSTR{cmp} = sprintf('%s (n = %i)',labelCMP,maxN);
        if cmp==noCMP
            set(hf,'CurrentAxes',ha(tra))
            legend(hp,legSTR)
        end
    end %loop transition
    sgtitle(labelAXI,'FontWeight','bold')
    
    %settings
    linkaxes(ha,'x')
    set(ha,'xlim',t([1,end]),'xtick',-20:5:20,'xgrid','on',...
        'unit','normalized',props.axi.trans1{:})
 
    set(ha(noTRA+1:end),'visible','off')
    
    
    set(ha,'ylim',[0,max(MIMA(:))*1.2])
   
    %% SAVE
    if opt.save.figs
        sname = fullfile(sPath,sprintf('DFF_trans_%s',...
            strrep(labelAXI,', ','-')));
        for k = 1:numel(opt.save.funsFig)
            fun = opt.save.funsFig{k};
            fun(hf,sname)
        end
        fprintf(' - Saved: %s.*\n',sname)
        %save individual data
        DD = mDat.trans;
        tab = {'data','stages'}';
        for axi = 1:noAXI
            labelAXI = strjoin(comp.axes{axi,2},', ');
            for cmp = 1:noCMP
                labelCMP = strjoin(comp.cmp{cmp,2},', ');
                %data
                L = mDat.stageDiff.labels{axi,cmp};
                D = DD.mean{axi,cmp};
                D = D';
                %append
                [rows,cols] = size(D);
                rows = (1:rows)+2;
                cols = (1:cols)+size(tab,2);
                tab{1,cols(1)} = labelCMP;
                tab(2,cols) = transitions;
                tab(rows,cols) = num2cell(D);
            end
        end
        %combine tables
        [r1,c1] = size(tab);
        [r2,c2] = size(tabPKS);
        Tab = cell(max([r1,r2]),c1+c2+1);
        Tab(1:r1,1:c1) = tab;
        [~,ind] = sort(tabPKS(2:end,1));
        tabPKS = tabPKS([1;ind+1],:);
        Tab(1:r2,c1+2:end) = tabPKS;
        %data
        sname = [sname,'.xlsx'];
        if exist(sname,'file')==2
            delete(sname)
        end
        xlswrite(sname,Tab)
        try
            xls_cellFit(sname)
        catch
        end
    end
end

%% PLOT IV (peak diff)
% clc; close all
fprintf('PLOT FIGURE IV (transition change)\n')
%figure/axes options
dx = [55,65,20]; dy = [45,65,50];
rows = ceil(sqrt(noTRA)); cols = ceil(noTRA/rows);
wid  = 200; hig = 200; %axes width/height
pos  = [200,200,cols*wid+[1,cols-1,1]*dx',rows*hig+[1,rows-1,1]*dy'];

%FIGURE LOOP
hf = figure('visible','off','unit','pixel','position',pos,props.fig{:});
movegui(hf,'center'); drawnow
set(hf,'visible','on');
ha = fig_createAxes(hf,[rows,cols],dx,dy,'pixel')';
set(ha,'unit','normalized')
mima = [inf,-inf]; %init
%init
x0 = 1:noAXI;
dx = linspace(-.5,.5,noCMP+2);
dl = mean(diff(dx));
dx([1,end]) = [];
hp = NaN(noCMP,noTRA);
mima = [inf,-inf]; %init
%init table
xlab = cellfun(@(x)strjoin(x,', '),comp.axes(:,2),'uniformoutput',false);
clab = cellfun(@(x)strjoin(x,', '),comp.cmp(:,2),'uniformoutput',false);
tab = cell(noCMP+2,noAXI*noTRA+1);
tab(:,1) = {''};
tab(3:end,1) = clab(:);
tmp = repmat(xlab(:),1,noTRA);
tab(1,2:end) = tmp(:);
tab(2,2:end) = repmat(transitions(:)',1,noAXI);

for cmp = 1:noCMP
    labelCMP = strjoin(comp.cmp{cmp,2},', ');
    for tra = 1:noTRA
        labelTRA = transitions{tra};
        set(hf,'CurrentAxes',ha(tra))
        %data
        x = x0+dx(cmp);
        tmp = shiftdim(mDat.trans.pks.y(:,cmp,tra,:),3);
        y = diff(tmp,[],1);
        mima = [min([mima(1);y(:)]),max([mima(2);y(:)])];
        
        %append
        for k = 1:numel(xlab)
            col = strcmpi(tab(1,:),xlab(k)) & strcmpi(tab(2,:),labelTRA);
            row = strcmpi(tab(:,1),labelCMP);
            tab(row,col) = {y(k)};
        end
        %plot
        hp(cmp,tra) = bar(x,y,'barwidth',0.8*dl,...
            'facecolor',colsCMP(cmp,:));
        hold on
        %text
        if cmp==1
            title(labelTRA);
        end
        if cmp==noCMP
            legend(hp(:,tra),cellfun(@(x)strjoin(x,', '),comp.cmp(:,2),...
                'uniformoutput',false))
        end
    end
end
%settings
linkaxes(ha,'xy')
xl = x0([1,end])+[-1,1]*2*dl;
yl = max(abs(mima))*[-1,1]*1.1;
set(ha,'xtick',x0,'xticklabel',xlab,'xlim',xl,'ylim',yl)
sgtitle('Transition change','FontWeight','bold')

%SAVE
if opt.save.figs
    sname = fullfile(sPath,'DFF_transDiff');
    for k = 1:numel(opt.save.funsFig)
        fun = opt.save.funsFig{k};
        fun(hf,sname)
    end
    fprintf(' - Saved: %s.*\n',sname)
    %save individual data
    sname = [sname,'.xlsx'];
    if exist(sname,'file')==2
        delete(sname)
    end
    xlswrite(sname,tab)
    try
        xls_cellFit(sname)
    catch
    end
end

%%
%% SAVE DATA
if opt.save.dat
    clear info %just in case
    sname = fullfile(sPath,sprintf('DFF_Sleep'));
    info.dimensions  = 'condition x area x stages/transition';
    info.sensors     = sensors(:)';
    info.areas       = areas(:)';
    info.stages      = stages(:)';
    info.transitions = trans(:)';
    save([sname,'.mat'],'info','mDat')
    fprintf(' - Saved: %s.mat \n\n',sname)
end

%%
if ~opt.save.figs
    fprintf('[\bFigures NOT saved]\b\n')
end