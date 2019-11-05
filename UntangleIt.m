n = 20;
X = (rand(n,1)+1)*100;
Y = (rand(n,1)+1)*100;
k=500;
delay = 0.5;

clf("reset")
hold on;

show(X,Y,n,k,delay);

%Plots multiple iterations of X and Y
function [] = show(X, Y, verticies, amount, delay)
    for i = 1:amount
        pause(delay);
        [X, Y] = shift(X,Y,verticies);
        plotThis(X,Y);  
    end
end

%Plots the Polygon
function  [] = plotThis(X, Y)
    pgon = polyshape(X,Y);
    plot(pgon);
end

%Does One iteration step of 1) a)
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

%Does the shifting and normalization in b)
function [X,Y] = shift(XIn, YIn,n)
    [X,Y] = iterate(XIn,YIn,n);
    meanX = mean(X);
    meanY = mean(Y);
    for i = 1:n
        X(i) = X(i) - meanX;
        Y(i) = Y(i) - meanY;
    end
    
    for i = 1:n
        X(i) = X(i) / norm(X);
        Y(i) = Y(i) / norm(Y);
    end
end
