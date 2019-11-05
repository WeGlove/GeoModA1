n = 20;
X = (rand(n,1)+1)*1;
Y = (rand(n,1)+1)*1;
k=500;
delay = 0;
%TODO How to get arbitrary X and Y

hold on;

show(X,Y,n,k,delay);

function [] = show(X, Y, verticies, amount, delay)
    for i = 1:amount
        pause(delay);
        [X, Y] = iterate(X,Y,verticies);
        plotThis(X,Y);  
    end
end

function  [] = plotThis(X, Y)
    pgon = polyshape(X,Y);
    plot(pgon);
end

function [X, Y] = iterate(XIn, YIn, n)
    X = zeros(n,0);
    Y = zeros(n,0);
    
    for i = 1:n
        if i < n
            X(i) = (XIn(i) + XIn(i+1))*0.5;
            Y(i) = (YIn(i) + YIn(i+1))*0.5;
        else
            X(i) = (XIn(i) + XIn(1))*0.5;
            Y(i) = (YIn(i) + YIn(1))*0.5;
        end
        
    end
end