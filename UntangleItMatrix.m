n = 4; % Number of vericies
k= 500; % Number of iterations
delay = 1; % Delay in showing new plots 

clf("reset")
hold on;

[X,Y] = getRandomGon(n);

matrix = getIterMatrix(n);
display(matrix);

show(X, Y, k, delay, matrix);

%Returns the matrix needed for a)
function [A] = getIterMatrix(n)
    A = sparse(n);
    for i = 1:n
        if i < n
            A(i,i) = 0.5;
            A(i,i+1) = 0.5;
        else
            A(i,i) = 0.5;
            A(i,1) = 0.5;
        end
    end
end

%Generates a new random polygon
function [X,Y] = getRandomGon(n)
    X = (rand(n,1)+1)*100;
    Y = (rand(n,1)+1)*100;
end

%Plots multiple iterations of X and Y
function [] = show(X, Y, amount, delay, matrix)
    for i = 1:amount
        pause(delay);
        X = matrix * X;
        Y = matrix * Y;
        plotThis(X,Y);  
    end
end

%Plots the Polygon
function  [] = plotThis(X, Y)
    pgon = polyshape(X,Y);
    plot(pgon);
end

%Does the shifting and normalization in b)
function [X,Y] = shift(XIn, YIn,n)
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
