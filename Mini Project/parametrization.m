clear, close all;

%[p,t]=loadmesh('samplemeshes/mushroom.off');%mac

[points,triangles]=loadmesh('samplemeshes\pig.off');%windows

plot_mesh_collected(points,triangles)

boundary_loops = get_boundary_loops(points,triangles);

C = get_boundary_conditions(@sample_unit_circle, boundary_loops{1}, points);
[points,triangles] = rearange_for_boundary(boundary_loops{1}, points, triangles);

% Rearanges the points such that the points corresponding to the
% <boundary> loop are at the very bottom. Updates the triangles accordingly
function [points, triangles] = rearange_for_boundary(boundary, points, triangles)
    for index = 1: length(boundary)
        [points,triangles] = swap_points(boundary(index), length(points)-length(boundary)+index, points, triangles);
    end
end

% Swaps two points in points and updates the triangles accordingly
function [points, triangles] = swap_points(id_first, id_second, points, triangles)
    temp = points(:,id_second);
    points(:,id_second) = points(:,id_first);
    points(:,id_first) = temp;
    for tri_index = 1:length(triangles)
        triangle = triangles(:, tri_index);
        triangles(:,tri_index) = update_index_in_triangle(id_first, id_second, triangle);
    end
end

% Helper function for <swap_points>. Updates a single triangle according to
% the swap
function triangle = update_index_in_triangle(id_first, id_second, triangle)
    for index = 1 : length(triangle)
        if triangle(index) == id_first
            triangle(index) = id_second;
        elseif triangle(index) == id_second
            triangle(index) = id_first;
        end
    end
end

%Creates the boundary conditions (the C vector in the formular). The
%conditions are set by a <sampler> function which is expected to return a
%2D vector.
function C = get_boundary_conditions(sampler, boundary, points)
    C = [zeros(length(points)-length(boundary),2) ;sampler(length(boundary))];
end

%Function to retrieve the points from a specific trinagle. Also works with
%index vectors!
function [indexes, tri_points] = get_points_from_triangle(triangle_index, points, triangles)
        indexes = triangles(:, triangle_index);
        tri_points = points(indexes);
end

%Returns boundary loops of the mesh
function boundaryloops = get_boundary_loops(points, triangles)
    nv=size(points,2);

    A=sparse(triangles([1 2 3],:),triangles([2 3 1],:),1,nv,nv);

    % The neighbors of any given vertex vi can be obtained
    % from the connectivity matrix C as
    %find neighbors
    %nvi= find(C(i,:));
    % A simple way to check if the mesh has a boundary would
    % be to compute the matrix B
    B=A-A';
    Bplus=spones(B==1);
    %Bminus=spones(B==-1);

    [ri,rj]=find(Bplus);%find returns row and col indices
    %[rl,rk]=find(Bminus); 

    bndry=rj;
    %rbndry=rl;
    %lbndry=ri;
    %rbl=[rl rj ri];


    boundaryloops=[];
    %if bndry
    %boundary=cell(1,size(bndry,2));
    ib=1;
    tic
    while bndry
        %bbegin=1;
        [~,bloop]=predecessor(Bplus,bndry(1));
        bndry=setdiff(bndry,bloop);
    
        boundaryloops{1,ib}=bloop';
        ib=ib+1
    end  
end


%An even simpler plotting function than the one given
function plot_mesh_collected(points, triangles)
    points=points'; %flip to colums for convenience in matrix multiplication

    figure, plotmesh(points,triangles');
    rzview('on') %allows for rotating zoom and move with the mouse (left-right-middle)

end

%Samples <amount> number of points on a unit circle
function samples = sample_unit_circle(amount)
    samples = [];
    for i = 0 : amount-1
        samples = [samples; cos(360/amount*i) sin(360/amount*i)];
    end 
end
