% Creates and operates the GUI
function [Tsidx,startIdx] = SplineGridGUI_bunt(Tsidx, data)
colact=[0.55 0.20 0.20]; % Bordeaux dark red
colina=[0 0 0 0.1]; % Transparent black
colsel='r'; % Light red
Tsidx = Tsidx(:)'; % Always got on my nerves

tf=figure('visible','off');

%% Create figure
uf = figure;
uf.Name=['Spline Grid Selection (Close to proceed)']; %#ok<NBRAK>
uf.Position = [100 100 1000 700];
uf.MenuBar='none';
uf.ToolBar='none';
set(uf,'units','pixel')

% Adds first and last spectrum
Tsidx=sort(unique([Tsidx 1 length(data.T)]));

fdata=guidata(tf);
fdata.ax1=axes('parent',uf,'units','pixel','position',[50 150 400 500]);
fdata.ax2=axes('parent',uf,'units','pixel','position',[550 150 400 500]);
fdata.panel=uipanel('parent',uf,'units','pixel','position',[50 20 900 90]);

str={'Click on lines to (de)select spline grid points (dark red).',...
    'Ctrl+Click to select the spectrum to start the optimization from (light red).',...
    'There has to be at least on selected start spectrum'};
fdata.label= uicontrol('style','text','string',str,'parent',fdata.panel,'position',[20 10 860 70]);

fdata.colact=colact;
fdata.colina=colina;
fdata.colsel=colsel;

guidata(tf,fdata);

%% Plots data

axes(fdata.ax1)
b = data.check(); % Checks if data is rectangular
fdata.b=b;

if b
    Dm=zeros(length(data.X{1}),length(data.T));
    for j=1:length(data.T)
        Dm(:,j)=data.D{j};
    end

    contourf(fdata.ax1,data.T,data.X{1},Dm,20,'Edgecolor','none'); % 20 = #contour lines

    % Visibility improvement with 'clim'
    set(gca, 'CLim' , data.clim);
    hold on
    % shading interp % If desired
    view([-90 90]);

    p1arr=zeros(length(data.T),1);
    for j=1:length(data.T)
        p1arr(j)=plot3(fdata.ax1,data.T(j)*[1 1],data.X{j}([1 end]),[1 1],'color', colina,'linewidth',3);
        hold on
    end
else
    set(fdata.ax1,'visible',off);
    set(fdata.ax2,'position',[50 150 900 500]);
end

axes(fdata.ax2);
p2arr=zeros(length(data.T),1);

%FancySchmancy Colors ?
Kolores = jet(ceil(length(data.T)*1.7)); % produces a colour gradient
colour_crop_ind = floor(length(data.T)*0.25); % crops the colour gradient
for j=1:length(data.T)
    p2arr(j)=plot3(fdata.ax2,data.X{j},data.T(j)*ones(size(data.X{j})),data.D{j},'color', Kolores(colour_crop_ind + j,:),'linewidth',2);
    hold on
end
set(gca, 'XDir', 'reverse');
axis tight;
%% modifies color of Tsidx

pact=zeros(length(data.T),1);
pact(Tsidx)=1;
pact(Tsidx(1))=2;

for j=1:length(data.T)
    if pact(j)==1
        if b
            set(p1arr(j),'color',colact);
        end
        set(p2arr(j),'color',colact);
    elseif pact(j)==2
        if b
            set(p1arr(j),'color',colsel);
        end
        set(p2arr(j),'color',colsel);
    else
        if b
            set(p1arr(j),'color', colina);
        end
        %         set(p2arr(j),'color',colina);
    end
end

fdata.p1arr=p1arr;
fdata.p2arr=p2arr;
fdata.pact=pact;
guidata(tf,fdata);


% set callbacks for line click
if b
    set(p1arr,'ButtonDownFcn',{@my_change,tf,uf})
end
set(p2arr,'ButtonDownFcn',{@my_change,tf,uf})



%%
set([fdata.ax1 fdata.ax2],'units','normalized')

axes(fdata.ax1)
box off
%disable zoom, rotate etc %Wäre es nicht anders besser?
fdata.ax1.Toolbar.Visible = 'off';

% wait for figure to be closed
uiwait(uf);

%write output
fdata=guidata(tf);
Tsidx=find(fdata.pact);
startIdx=find(nonzeros(fdata.pact)==2,1);
disp(Tsidx);
% obj=getappdata(tf,'content');
delete(tf);

end



function my_change(hobj,evt,tf,uf) %#ok<INUSL>

fdata=guidata(tf);
%         if b
%             set(p1arr(j),'color',colina)
%         end
%         set(p2arr(j),'color',colina);


modifiers = get(uf,'CurrentModifier');
CtrlPressed  = ismember('control', modifiers);

idx=[];
if fdata.b
    idx=find(fdata.p1arr==hobj);
end
if isempty(idx)
    idx=find(fdata.p2arr==hobj);
end

idxsel=find(fdata.pact==2);


% set states
if CtrlPressed
    if idx~=idxsel
        fdata.pact(idx)=2;
        fdata.pact(idxsel)=1;
    end
else
    if fdata.pact(idx)==1
        fdata.pact(idx)=0;
    elseif fdata.pact(idx)==2
        % selected spec not deselectable
    elseif fdata.pact(idx)==0
        fdata.pact(idx)=1;
    end
end


% FancyColors 
Kolores = jet(ceil(length(fdata.pact)*1.7));
colour_crop_ind = floor(length(fdata.pact)*0.25);

% Set colors
for j=1:length(fdata.pact)
    if fdata.pact(j)==1
        if fdata.b
            set(fdata.p1arr(j),'color',fdata.colact);
        end
        set(fdata.p2arr(j),'color',fdata.colact);
    elseif fdata.pact(j)==2
        if fdata.b
            set(fdata.p1arr(j),'color',fdata.colsel);
        end
        set(fdata.p2arr(j),'color',fdata.colsel);
    else
        if fdata.b
            set(fdata.p1arr(j),'color', fdata.colina);
        end
        set(fdata.p2arr(j),'color',Kolores(colour_crop_ind + j,:));
    end
end
guidata(tf,fdata)

end