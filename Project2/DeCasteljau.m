clf('reset');
hold on;

%Input points
global points;
points = [0 1; 1 5; 5 6; 7 5; 9 1];
t = 0.5;

%set axis limits
setAxisLimits(points);

%bezier curve smoothness
global detail;
detail = 100;


%draw Bezier Curve
plotBezierCurve();


%visualize deCasteljau
%visualizeCasteljau(points(:,1), points(:,2), t);



% ============ a) ============

function [] = plotBezierCurve()
    global detail;
    curve = zeros(detail,2);
    step = 1/detail;
    
    %calculate the bezier values in [0,1] with steps of ratio 1/detail
    for k = 0:step:1
        p = deCasteljau(k);
        curve(round(k*detail+1),1) = p(1);
        curve(round(k*detail+1),2) = p(2);
    end
    
    %plot the calculated points on the bezier curve
    plot(curve(:,1), curve(:,2));
     
end


% ============ b) ============

function [] = visualizeCasteljau(X,Y,t)

    global points;
    
    %draw initial points before first iteration
    drawIterPointsAndLines(X,Y, false);
    
    %draw iteration steps
    npoints = size(X,1);
    savepoints = points;
    
    for i =1:npoints
        numnewp = npoints-i;
        itpoints = zeros(numnewp-i, 2);  %save the points for each iteration
        for l = 1:(numnewp)
            itpoints(l,1) = convexCombination(savepoints(l,1),savepoints(l+1,1),t);
            itpoints(l,2) = convexCombination(savepoints(l,2),savepoints(l+1,2),t);
        end
        drawIterPointsAndLines(itpoints(:,1),itpoints(:,2), true);
        
        savepoints = itpoints;
    end
    
    %draw bezier curve
    plotBezierCurve();
end


% =========== De Casteljau Algorithmn ======
function [p] = deCasteljau(t)
    global points;
    X = points(:,1);
    Y = points(:,2);
    
    npoints = size(X,1);
    savepoints = points;
    
    p = zeros(2,1);
    
    for i =1:(npoints-1)
        numnewp = npoints-i;
        itpoints = zeros(numnewp-i, 2);  %save the points for each iteration
        for l = 1:(numnewp)
            itpoints(l,1) = convexCombination(savepoints(l,1),savepoints(l+1,1),t);
            itpoints(l,2) = convexCombination(savepoints(l,2),savepoints(l+1,2),t);
        end
        
        %finished? then return last the determined point
        if(i == npoints-1)
          p(1) = itpoints(1);
          p(2) = itpoints(2);
        end
        
        savepoints = itpoints;
    end
end


%convex combination of control points
function [mean] = convexCombination(a1, a2, t)
    mean = (1-t)*a1 + t * a2;
end



% =========== Helper Functions ==========

%Draws the given points and lines between them
%The points have to be in the correct oder
%bool b: change point style: true = no dots    
function drawIterPointsAndLines(X,Y,b)
   plot(X,Y);
   if(b)
       scatter(X,Y, 'x');
   else
       scatter(X,Y, 'filled');
   end
   
end

%redefines the axis limits to display the graph better
function [] = setAxisLimits(P)
    X = P(:,1);
    Y = P(:,2);
    
    minx = min(X);
    maxx = max(X);
    miny = min(Y);
    maxy = max(Y);
    
    xlim([(minx-1) (maxx+1)])
    ylim([(miny-1) (maxy+1)])
end

