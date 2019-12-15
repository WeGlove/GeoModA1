content = get_data("SineRandom.txt");
%content = get_data("CamelHeadSilhouette.txt");
%content = get_data("MaxPlanckSilhouette.txt");

plot(content(:,1),content(:,2))

function content = get_data(file_name)
    fileID = fopen(file_name,'r'); 
    formatSpec = '%f'; %Specifies Data Type f == Float
    size = [2,Inf]; %Specifies to make it a 2 column vector
    content = fscanf(fileID, formatSpec, [2,Inf]);
    content = content'; %Transpose to get correct vector
end

function [X, Y] = format(vector)
    for i = 1:length(vector)
        if mod(i,2) == 1
            X = [X;vector(i)];
        else
            Y = [Y;vector(i)];
        end
    end
end