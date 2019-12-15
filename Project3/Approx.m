curvy_factor = 3;
factor = 12;

content = get_data("SineRandom.txt");
%content = get_data("CamelHeadSilhouette.txt");
%content = get_data("MaxPlanckSilhouette.txt");

clf('reset');
hold on;

plot(content(:,1),content(:,2), 'Color', 'b'); %Draw he original curve

content = mva(content, factor);
points = getCtrlPts(content, curvy_factor);
plot(points(:,1),points(:,2),'Color', 'r');
draw_curves(points, 'g', 2);

%SAMPLING

function sample = mva(points, factor)
    sample = [];
    for i = 1 : length(points)/factor
        avg = 0;
        k=0;
        for j = 1 : factor
            if i+j < length(points)
                avg = avg + points((i-1)*factor+j,:);
                k = k+1;
            end
        end
        avg = avg / k;
        sample = [sample; avg];
    end
end

%DATA IO
function content = get_data(file_name)
    fileID = fopen(file_name,'r'); 
    formatSpec = '%f'; %Specifies Data Type f == Float
    size = [2,Inf]; %Specifies to make it a 2 column vector
    content = fscanf(fileID, formatSpec, [2,Inf]);
    content = content'; %Transpose to get correct vector
end

%BEZIER CURVES
function draw_curves(points, color, width)
    tsampling = 0:0.01:1;
    for i = 0 : (length(points)+4)/4 - 1 - 1
        curve = bezier(points(i*4+1:i*4+4,:), tsampling, color, width);
        plot(curve(:,1),curve(:,2), 'Color', color, 'LineWidth', width);
    end
end

function curve = bezier(points, tsampling, color, width)    
    curve = [];
    for j = 1:length(tsampling)
        controlPts = points;
        %plot(controlPts(:,1),controlPts(:,2),'x');
        %line(controlPts(:,1),controlPts(:,2),'Color',color, 'LineWidth', width);
        for i = 1:(length(controlPts)-1)
            controlPts = de_casteljau(controlPts, tsampling(j));
        end
        curve = [curve; controlPts];
    end
end

function intermediates = de_casteljau(points,t)
     intermediates = zeros(length(points)-1,2);
     for i = 1:(length(points)-1)
        intermediates(i,1) = points(i,1) .* (1-t) + points(i+1,1) .* (t);
        intermediates(i,2) = points(i,2) .* (1-t) + points(i+1,2) .* (t);
     end
end

function ctrlPts = getCtrlPts(intersections, curvy_factor)
    ctrlPts = [];
    for i = 1 : length(intersections)
        if i == 1
            interToRight = intersections(i+1,:) - intersections(i,:);
            R = intersections(i,:) + interToRight / curvy_factor;
            ctrlPts = [ctrlPts; intersections(i,:); R];
        elseif i == length(intersections)
            interToLeft = intersections(i-1,:) - intersections(i,:);
            L = intersections(i,:) + interToLeft / curvy_factor;
            ctrlPts = [ctrlPts; L; intersections(i,:)];
        else
            [L,R] = get_par_points(intersections(i,:), intersections(i-1,:), intersections(i+1,:), curvy_factor);
            ctrlPts = [ctrlPts; L; intersections(i,:); intersections(i,:); R];
        end
    end
end

function [L,R] = get_par_points(intersection, left, right, curvy_factor)
    lefttoright = right - left;
    DistLToI = norm(intersection - left);
    DistRToI = norm(intersection - right);
    if DistLToI < DistRToI
        Dist = DistLToI / curvy_factor;
    else
        Dist = DistRToI / curvy_factor;
    end
    DistLToR = norm(lefttoright);
    R = lefttoright / DistLToR * Dist + intersection;
    L = -lefttoright / DistLToR * Dist + intersection;
end