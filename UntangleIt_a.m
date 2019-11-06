n = 25;
[x,y] = createPolygon(n, true); %True = without intersections
fun = @fun_ex1b; % choose fun : Ex1a = @fun_ex1a or Ex1b = @fun_ex1b
amount = 100;
delay = 0.1;
clf("reset")
hold on;
show(x, y, fun, amount, delay)


%Plots multiple iterations of X and Y
function [] = show(x, y, fun, amount, delay)
    for i = 1:amount
        pause(delay);
        [x, y] = fun(x,y);
        PlotPolygon(x,y);  
    end
end

%Input: n = edges Bool: True = without intersections
%create a Polygon with x,y are vectors
function [x,y] = createPolygon(n, intersections)
tripel = zeros(n,3);
if intersections == true
%random_int
    tripel(:,1) = randi([1,100],n,1); %random x
    tripel(:,2) = randi([1,100],n,1); %random y
% In general, N random numbers in the interval (a,b)
% r = a + (b-a)*rand(N,1).
    % tripel(:,1) = 1+(100-1)*rand(n,1);
    % tripel(:,2) = 1+(100-1)*rand(n,1);
    center_point = (sum(tripel))/n;
    tripel(:,3) = atan2(tripel(:,1)-center_point(1),tripel(:,2)-center_point(2)); %angle    
    tripel = sortrows(tripel,3); %sort the points by angle
else 
   %random_int
    tripel(:,1) = randi([1,100],n,1); %random x
    tripel(:,2) = randi([1,100],n,1); %random y
% In general, N random numbers in the interval (a,b)
% r = a + (b-a)*rand(N,1).
    % tripel(:,1) = 1+(100-1)*rand(n,1);
    % tripel(:,2) = 1+(100-1)*rand(n,1);
end
x = tripel(:,1);
y = tripel(:,2);
end

%Plots the Polygon
function [] = PlotPolygon(x,y)
polygon = polyshape(x,y);
plot(polygon);
end

%Does One iteration step of 1) a)
function [x,y] = fun_ex1a(x,y)
x=subfun_ex1a(x);
y=subfun_ex1a(y);
end

% function 1 a)
function [y] = subfun_ex1a(x)
n = length(x);
y = zeros(n,1);
for k = 1:n
    if k < n
        y(k) = 1/2 * (x(k)+x(k+1));
    else
        y(k) = 1/2 * (x(k)+x(1));
    end
end
end

%shift and normalize 1) b)
function [x,y] = fun_ex1b(x,y)
x = subfun_ex1a(x)- mean(x);
x = x/norm(x);
y = subfun_ex1a(y)- mean(y);
y = y/norm(y);
end