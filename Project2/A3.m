%points = [9 9;9 3;0 9;2 0];
points = [2 0; 0 9; 9 3; 9 9];
%points = [1 1; 3 1];
global plot_progress; 
plot_progress = 1;

clf('reset');
hold on;

t1 = 0.1;
t2 = 0.8;

curve = bezier(points, 0:0.01:1, 'k',2);

[lower_points, points1] = subdivide(points, t1);
[points3, points2] = subdivide(lower_points, t2);

curve1 = bezier(points1, 0:0.01:1, 'r',1);
curve2 = bezier(points2, 0:0.01:1, 'g',1);
curve3 = bezier(points3, 0:0.01:1, 'b',1);

plot(curve(:,1), curve(:,2), 'k', 'LineWidth',6);
plot(curve1(:,1), curve1(:,2), 'r', 'LineWidth',4);
plot(curve2(:,1), curve2(:,2), 'g', 'LineWidth',4);
plot(curve3(:,1), curve3(:,2), 'b', 'LineWidth',4);


function curve = bezier(points, tsampling, color, width)
    global plot_progress;
    
    curve = [];
    for j = 1:length(tsampling)
        controlPts = points;
        plot(controlPts(:,1),controlPts(:,2),'x');
        line(controlPts(:,1),controlPts(:,2),'Color',color, 'LineWidth', width);
        for i = 1:(length(controlPts)-1)
            controlPts = de_casteljau(controlPts, tsampling(j));
        end
        curve = [curve; controlPts];
    end
end

function [lower_points, upper_points] = subdivide(points, t)
    controlPts = points;
    lower_points = [];
    upper_points = [];
    for i = 1:(length(controlPts)-1)
            upper_points = [upper_points; controlPts(1,1) controlPts(1,2)];
            lower_points = [lower_points; controlPts(length(controlPts(:,1)),1) controlPts(length(controlPts(:,1)),2)];
            controlPts = de_casteljau(controlPts, t);
    end
    upper_points = [upper_points; controlPts(1,1) controlPts(1,2)];
    lower_points = [lower_points; controlPts(length(controlPts(:,1)),1) controlPts(length(controlPts(:,1)),2)];
    lower_points = flip(lower_points);
end
     
function intermediates = de_casteljau(points,t)
     intermediates = zeros(length(points)-1,2);
     for i = 1:(length(points)-1)
        intermediates(i,1) = points(i,1) .* (1-t) + points(i+1,1) .* (t);
        intermediates(i,2) = points(i,2) .* (1-t) + points(i+1,2) .* (t);
     end
end
