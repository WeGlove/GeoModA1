clf("reset")
hold on;

start = -2;
stop = 2;

%Global Variables
global a;
global m;
a = 1;
m = 3;
global X;
X = sample(10, start, stop);
global Y;
Y = power_fun(X);
global base;
base = get_vandermonde(X,Y);


%Plot the function
fplot(@power_fun, [-3,3]);

%Plot the markers

plotMarkers(X,Y);

%Plot Lagrange

fplot(@lagrange, [-2,2])

%Plot everything else

fplot(@power_base, [-2,2]);


function [L] = power_base(X)
    global base;
    L = zeros(length(X));
    for i = 1:length(base)
        L = L + base(i) * (X.^(i-1));
    end
end


function [L] = get_vandermonde(X, Y)
    vandermonde = zeros(length(X));
    for i = 1:length(X)
        for j = 1:length(X)
            vandermonde(i,j) = X(i)^(j-1);
        end
    end
    L = linsolve(vandermonde,Y);
end

function [l] = lagrange_polynom(x, i)
    global X;
    l = 1;
    for j = 1: length(X)
        if j ~= i
            l = l .*((x - X(j)) ./ (X(i) - X(j)));
        end
    end
end

function [out] = lagrange(x)
    global Y;
    out = zeros(length(x),1);
    for i = 1 : length(Y)
        out = out + Y(i) .* lagrange_polynom(x,i);
    end
end

function [y] = power_fun(x)
    global a;
    global m;
    y = a.*(x.^m);
end

function [out] = sample(n, start, stop)
    out = rand(n,1)*(stop-start) + start;
    out = sort(out);
end

function [] = plotThis(X,Y)
    plot(X,Y);
end

function [] = plotMarkers(X,Y)
    plot(X,Y,'o');
end