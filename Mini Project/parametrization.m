clear, close all;

%[p,t]=loadmesh('samplemeshes/mushroom.off');%mac

[points,triangles]=loadmesh('samplemeshes\pig.off');%windows

plot_mesh_collected(points,triangles)

circle = sample_unit_circle(10);

boundary_loops = get_boundary_loops(points,triangles);

function C = get_boundary_conditions(amount, boundary_amount, sample_function)
    C = [zeros(amount-boundary_amount) sample_function(boundary_amount)];
end

function [U C] = sort_boundary_points(boundary,points)
    new_points = zeros(length(points))
    for i = 1:length(points)
        if new_points
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
