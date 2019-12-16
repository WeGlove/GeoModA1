global points;
points = [0 0;9 9;0 9; 9 0;5 1; 4 3];
%points = [0 5;4 5;6 8;7 5];
global plot_progress;
plot_progress = 0;
clf('reset');
hold on;

global gradient;
gradient =0.5;

curve = bezier(0:0.01:1);
ha = plot(curve(:,1), curve(:,2),  '.-', 'LineWidth',10, 'DisplayName',' 0.5');
plot_progress = 1;

curvVector = get_all_curvature(0:0.01:1);
cd = uint8(interpolate(curvVector)');
% colors are stored in the first three rows, each column corresponds to a 
% the color of a sampling point on the curves [r g b]'
cd(4,:)=0;
drawnow
set(ha.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

function colors = interpolate(vector)
    global  gradient
    colors = [];
    for i = 1:length(vector)
        if 0 <= vector(i) &&  vector(i) < gradient
            reverse  = 255*(gradient-vector(i))/gradient;
            colors = [colors; 255*vector(i)/gradient+reverse reverse reverse 1];
       elseif vector(i) > gradient
            colors = [colors; 255 0 0 1];
        elseif -gradient < vector(i) && vector(i) < 0
            reverse  = 255*(gradient+vector(i))/gradient;
            colors = [colors; reverse reverse (255*-vector(i))/gradient+reverse 1];
        else
            colors = [colors; 0 0 255 1];
        end
    end
end

function curve = bezier(tsampling)
    global points;
    global plot_progress;
    
    curve = [];
    for j = 1:length(tsampling)
        controlPts = points;
        for i = 1:(length(controlPts)-1)
            if plot_progress ~= 0
                plot(controlPts(:,1),controlPts(:,2),'x');
                line(controlPts(:,1),controlPts(:,2))
            end
            controlPts = de_casteljau(controlPts, tsampling(j));
        end
        if plot_progress ~= 0
                plot(controlPts(:,1),controlPts(:,2),'x');
                line(controlPts(:,1),controlPts(:,2))
            end
        curve = [curve; controlPts];
    end
end

function curvature = get_all_curvature(tsampling)
    curvature = [];
    for i = 1:length(tsampling)
        curvature = [curvature; get_curvature(tsampling(i))];
    end
end

function curvature = get_curvature(t)
    first  = get_derivative(1,t);
    second = get_derivative(2,t);
    curvature = first(1) * second(2) - second(1) * first(2);
    curvature = curvature / (first(1)^2 + first(2)^2)^(3/2);
end

function derivative = get_derivative(order, t)
    global points;
    ctlPts = points;
    for i = 1:(length(points)-1-order)
        ctlPts = de_casteljau(ctlPts,t);
    end
    for i = 1:order
        ctlPts = de_casteljau_der(ctlPts);
    end
    n = 1;
    for i = length(points)-(order-1) : length(points)
        n = n*i;
    end
    derivative = ctlPts .* n;
end

function intermediates = de_casteljau_der(points)
     intermediates = zeros(length(points)-1,2);
     for i = 1:(length(points)-1)
        intermediates(i,1) = points(i,1) .* (-1) + points(i+1,1) .* (1);
        intermediates(i,2) = points(i,2) .* (-1) + points(i+1,2) .* (1);
     end
end
     
function intermediates = de_casteljau(points,t)
     intermediates = zeros(length(points)-1,2);
     for i = 1:(length(points)-1)
        intermediates(i,1) = points(i,1) .* (1-t) + points(i+1,1) .* (t);
        intermediates(i,2) = points(i,2) .* (1-t) + points(i+1,2) .* (t);
     end
end
