function funwithcurves
% basic interface for manipulating control polygons and curves
% you need to call


% create a new figure 
figure('Name','Fun with Curves','NumberTitle','off','color','white')


%% define variables used and stroe them in handle 
myhandles=[];
myhandles.tsampling=linspace(0,1,250)';% number of parameter values used for drawing the curves 
myhandles.controlPts=zeros(0,2);
myhandles.control_polygon=0;
myhandles.curve=0;
myhandles.dragStyle=0;
myhandles.dragedControPoint=0;
myhandles.axis=0;
myhandles.checked=0;
myhandles.t=1/2; %for intermediate deCasteljau


%% design a basic  interface
% seup axis with a grid layout for drawing 
hax = gca; %get handle of current axes
set(hax,'Unit','normalized','Position',[0 0 1 1]); % make axis fill figure
axis manual; % disable automatic resizing of axis
grid on 
grid minor; % draw a fine grid 

hold on; % keep waiting for more input in the figure

myhandles.axis=hax;
%guidata(gcf,myhandles) 



% Parameter t slider
myhandles.t_slider=uicontrol('Units','normalized', ...
	'BackgroundColor',get(gcf,'Color'), ...
    'position',[.25 .05 .5 .05],...
	'Min',0,'Max',1,'Value',myhandles.t,...
	'Sliderstep',[.05,.05],...
	'Callback',@slider_update,...
	'Style','slider','Enable','off');

% Text for displaying hte value of the slider on the figure
myhandles.t_slider_text=uicontrol('Units','normalized', ...
	'BackgroundColor','white', ...
	'HorizontalAlignment','center', ...
	'Position',[.45 .1 .1 .05],...
	'String','t=0.50', ...
	'Style','text');

%CHECKBOX , when checked the slider for drawing de intermediate steps of
%decastljau algorithm at time t given by the slider
checkboxTag = sprintf('Check for intermediate \nTooltip decasteljau steps');
myhandles.t_slider_chekbox=uicontrol('Units','normalized', ...
	'BackgroundColor','white', ...
	'HorizontalAlignment','center', ...
	'Position',[.2 .06 .04 .04],...
    'TooltipString',checkboxTag,... 
	'Style','checkbox','Callback',@checkbox_Callback);


% help text 
uicontrol('Units','normalized', ...
	'BackgroundColor','white', ...
    'ForegroundColor',[.94,.5,.94],...
	'HorizontalAlignment','center', ...
	'FontWeight','bold', ...
	'Position',[.15,.0,.7,.03],...
	'String','mouse button left=Adds, right=drags, middle = removes, control points', ...
	'Style','text');

% reset clear the screen and re-initialize data
uicontrol('Units','normalized', ...
	'BackgroundColor','white', ...
    'ForegroundColor','red',...
	'String','Reset', ...
	'Position',[.9,.0,.1,.04], ...
	'Callback',@reset,...
	'Style','pushbutton');

%save myhandle using guidata 
guidata(gcf,myhandles) 


%d%efine Mouse callbacks
set(gcf,'WindowButtonDownFcn',@mousePress);
set(gcf,'WindowButtonMotionFcn',@mouseDrag);
set(gcf,'WindowButtonUpFcn',@mouseRelease);
%%% Callback Functions %%%

% Reset button callback
function reset(varargin) 
%global myhandles;
myhandles = guidata(gcbo);

myhandles.controlPts=zeros(0,2);
myhandles.dragStyle=0;
myhandles.dragedControPoint=0;




%deletes curve from figure
if myhandles.curve~=0
    delete(myhandles.curve)
    myhandles.curve=0;
end


%deletes control polygon from figure
if myhandles.control_polygon~=0
    delete(myhandles.control_polygon)
    myhandles.control_polygon=0;
end

guidata(gcbo,myhandles) 
%%

%checkbox callback
function checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

myhandles = guidata(gcbo);
if (get(hObject,'Value') == get(hObject,'Max'))
    set(myhandles.t_slider,'Enable','on')
    myhandles.checked=1;
	
else
	set(myhandles.t_slider,'Enable','off')
    myhandles.checked=0;
end
guidata(gcbo,myhandles) 
 

%slider callback
function slider_update(varargin) 
%global myhandles;
myhandles = guidata(gcbo);

%% Read t value
myhandles.t = get(myhandles.t_slider,'Value');
set(myhandles.t_slider,'Value',myhandles.t); %update the value in handles
% display current value in the text box
set(myhandles.t_slider_text,'String',sprintf('t=%4.2f',myhandles.t)); 

%% make the text box follow the slider :)
txtpos=get(myhandles.t_slider_text,'position');
sliderpos=get(myhandles.t_slider,'position');
txtpos(1)=sliderpos(1)+myhandles.t*(sliderpos(3)-.05);
set(myhandles.t_slider_text,'position',txtpos);
%% update plot
guidata(gcbo,myhandles) 
updatePlot;



%%=======================================================================%%
%%                            Mouse actions                              %%
%%=======================================================================%%


%% Mouse button pressed
function mousePress(varargin)

% get the variables stored in handles
myhandles = guidata(gcbo);
closeEnough=.02;

%decide action based on which mouse button is pressed
%case 1
if strcmp(get(gcf,'SelectionType'),'normal') % left button: new control point being added 
  currentpos=get(myhandles.axis,'CurrentPoint'); %get the location of the current mouse position
  currentpos=currentpos(1,1:2);   % keep only what we need (see matlab help for currentpoint under figure properties)
  if min(currentpos) >= 0 & max(currentpos) <=1  % within the axis    
    myhandles.controlPts(end+1,:)=currentpos;% append to the tail of the array
    n=size(myhandles.controlPts,1);
    myhandles.dragStyle=1;
    myhandles.dragedControPoint=n;
    
    guidata(gcbo,myhandles) %update data in handles    
  end
end

% case 2
if strcmp(get(gcf,'SelectionType'),'alt') % right button: control point being dragged
  currentpos=get(myhandles.axis,'CurrentPoint');
  currentpos=currentpos(1,1:2);  
  if min(currentpos) >= 0 && max(currentpos) <=1 % within the axis
    % get closest control point to the current mous position
    [mindist,minindx]=min((currentpos(1)-myhandles.controlPts(:,1)).^2+(currentpos(2)-myhandles.controlPts(:,2)).^2);
    if mindist < closeEnough % only update when the the user clicks neer control points
      myhandles.dragStyle=1; 
      myhandles.dragedControPoint=minindx;
      myhandles.controlPts(minindx,:)=currentpos;
      
      guidata(gcbo,myhandles)       % update data in handles
    end
  end
end

% case 3
if strcmp(get(gcf,'SelectionType'),'extend') % middle button: clocked point being deleted 
  currentpos=get(myhandles.axis,'CurrentPoint');
  currentpos=currentpos(1,1:2);  
  if min(currentpos) >= 0 && max(currentpos) <=1 % within the axis
    % get closest control point to the current mous position
    [mindist,minindx]=min((currentpos(1)-myhandles.controlPts(:,1)).^2+(currentpos(2)-myhandles.controlPts(:,2)).^2);
    if mindist < closeEnough 
      myhandles.controlPts(minindx,:)=[];
      
      guidata(gcbo,myhandles)  % update data in handles      
    end
  end
end
%% update the plot in the figure based on handles updates from the cases above
updatePlot;



% Mouse button released
function mouseRelease(varargin)
%global myhandles;
myhandles = guidata(gcbo);

myhandles.dragStyle=0;
myhandles.dragedControPoint=0;
guidata(gcbo,myhandles) 



% Mouse being moved (after clicking)
function mouseDrag(varargin)
%global myhandles;
myhandles = guidata(gcbo);

if myhandles.dragStyle
  p=get(myhandles.axis,'CurrentPoint');
  p=p(1,1:2);
  if min(p) >= 0 & max(p) <=1 
    myhandles.controlPts(myhandles.dragedControPoint,:)=p;
    guidata(gcbo,myhandles) 

    updatePlot;
  end 
end 

%%=======================================================================%%
%%                             update plot                               %%
%%=======================================================================%%
  
%%% Internal Functions %%%
function updatePlot
%global myhandles;
myhandles = guidata(gcbo);

ncp=size(myhandles.controlPts,1);%number of control points

%decide when to draw/update the control polygon and the curve
if ncp>0
  % draw polygon
  if myhandles.control_polygon==0  %draw polygon for the first time
    myhandles.control_polygon=plot(myhandles.controlPts(:,1),myhandles.controlPts(:,2),'-o',...
        'color',[.4,.4,.4],...
        'LineWidth',2,...
    'MarkerSize',15,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor',[0.5,0.5,0.75]);
  else
    set(myhandles.control_polygon,'Xdata',myhandles.controlPts(:,1));
    set(myhandles.control_polygon,'Ydata',myhandles.controlPts(:,2));
  end
end

if ncp>1
 
  % draw curve
   %%============================your code==================================== 
  %% call your function for computing the curve coordinates here
  % it takes as inpust the control points: myhandles.controlPts 
  % and the parametric sampling: myhandles.tsampling
  %and outputs the coorinates of the whole curves in the array: mycurve
  
  mycurve= my_nice_function(myhandles.controlPts,myhandles.tsampling);
  %%
  
  if myhandles.curve==0 %draw curve for the first time  
    myhandles.curve=plot(mycurve(:,1),mycurve(:,2),'-r','Linewidth',2);  
  else %update curve
    set(myhandles.curve,'Xdata',mycurve(:,1));
    set(myhandles.curve,'Ydata',mycurve(:,2));
  end
end

%% two special case which happen when deleting points
if ncp<2
    if myhandles.curve~=0
        delete(myhandles.curve)
        myhandles.curve=0;
    end
end

if ncp<1
    if myhandles.control_polygon~=0
        delete(myhandles.control_polygon)
        myhandles.control_polygon=0;
    end
end



if myhandles.checked % if the checkbox is checked show the intermediate steps
 
 %%============================your code==================================== 
 %%Your code for drawing the intermediate decasteljau steps (4.1.b) can be put here   
 
 
 %%=========================================================================
end

%%

%update the handles data
guidata(gcbo,myhandles) 

 

%%
function mycurve= my_nice_function(controlPts,tsampling);
%define your function here for 4.1 a
 %%============================your code==================================== 