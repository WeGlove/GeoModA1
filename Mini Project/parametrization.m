clear, close all;

[points,triangles]=loadmesh('samplemeshes/head.off');%windows

normals = (get_normals(points, triangles) + 1)/2;
%plot_mesh_collected(points,triangles,normals)

boundary_loops = get_boundary_loops(points,triangles);

[C,M, arr_points, arr_triangles] = get_linear_equation(boundary_loops, points, triangles);
spy(M);
M = full(M);
boundary_loops = get_boundary_loops(arr_points, arr_triangles);

U = linsolve(M,full(C));
c_col = get_chess_colors(0.5, U, arr_triangles) * 255;
out = to_mesh(U);
%plot_mesh_collected(arr_points,arr_triangles, c_col')
plot_mesh_collected(out',arr_triangles,c_col')
%to_off("yee.off", out, arr);

function colors = get_chess_colors(mo, points, triangles)
    colors = [];
    for i = 1 : length(triangles)
        triangle_point = get_points_from_triangle(i,points, triangles);
        
        centroid = points(triangle_point(1),:) + points(triangle_point(2),:) + points(triangle_point(3),:);
        centroid = centroid / 3;
        if (mod(centroid(1),mo)/mo < 0.5 && mod(centroid(2),mo)/mo < 0.5) || (mod(centroid(1),mo)/mo > 0.5 && mod(centroid(2),mo)/mo > 0.5)
        %if mod(centroid(1) + centroid(2), mo)/mo < 0.5
            colors = [colors; [0 0 0]];
        else
            colors = [colors; [1 1 1]];
        end
    end
end


function normals = get_normals(points, triangles)
    normals = zeros(3,length(triangles));
    for i = 1 : length(triangles)
        triangle = triangles(:,i);
        a = points(:,triangle(1));
        b = points(:,triangle(2));
        c = points(:,triangle(3));
        normal = cross(b - a, c - a);
        normal = normal / norm(normal);
        normals(:,i) = normal;
    end
end

function to_off(file_name,points, triangles)
    fileID = fopen(file_name, 'w');
    points = points';
    fprintf(fileID, "OFF\n")
    fprintf(fileID, "%d %d %d\n",[length(points), length(triangles), 0])
    fprintf(fileID,'%f %f %f\n',points);
    fprintf(fileID,'3 %d %d %d\n',triangles);
end

function points = to_mesh(points)
    new = zeros(length(points),3);
    for i = 1 : length(points)
        new(i,:) = [points(i,:) 0];
    end
    points = new;
end

function [C,M,arr_points,arr_triangles] = get_linear_equation(boundary_loops, points, triangles)
    loop = 1;
    C = get_boundary_conditions(@sample_unit_circle, boundary_loops{loop}, points);
    [arr_points, arr_triangles] = rearange_for_boundary(boundary_loops{loop}, points, triangles);
    M = get_M(boundary_loops{loop}, arr_points, arr_triangles, 0, 1);
end

function M = get_M(boundary, points, triangles, lambda, mu)
    rings = form_triangle_rings(points,triangles);
    for i = 1 : length(rings)
        ring = rings{i};
        isRing = ring(2,1) == ring(3,length(ring(1,:)));
    end
    M_alpha = get_M_alpha(rings, boundary, points);
    display("Got Alpha");
    M_chi = get_M_chi(rings, boundary, points);
    E = full(M_chi);
    display("Got Chi");
    M = M_alpha * lambda + M_chi * mu; 
    for i=1:length(boundary)
        M(length(points)-length(boundary)+i,length(points)-length(boundary)+i) = 1;
    end
end

function M_alpha = get_M_alpha(rings, boundary, points)
    M_alpha = sparse(length(points),length(points));
    for index_row = 1 : length(points) - length(boundary)
        ring = rings{index_row};
            
        accumulate = 0;
        prime_point = points(:, index_row);
        
        for index_ring = 1 : length(ring(1,:))
            middle_point_index = ring(2,index_ring);
            alpha_point_index = ring(3,index_ring);
            if index_ring == 1
                a = length(ring(1,:));
            else
                a = index_ring-1;
            end
            if index_ring == length(ring(1,:))
                b = 1;
            else
                b = index_ring+1;
            end
            beta_point_index = ring(2,a);
            
            middle_point = points(:,middle_point_index);
            alpha_point = points(:,alpha_point_index);
            beta_point = points(:,beta_point_index);

            alpha = get_angle(alpha_point, prime_point, middle_point);
            beta = get_angle(beta_point, prime_point, middle_point);
            
            if ring(3,index_ring) ~= ring(2,b)
                coa = 0;
            else
                coa = cot(alpha);
            end
            if ring(2,index_ring) ~= ring(3,a)
                cob = 0;
            else
                cob = cot(beta);
            end
            value = coa + cob;
            M_alpha(index_row,middle_point_index) = value;
            accumulate = accumulate + value;
        end
            M_alpha(index_row, index_row) = - accumulate;
    end
end

function M_chi = get_M_chi(rings, boundary, points)
    M_chi = sparse(length(points),length(points));
    for index_row = 1 : length(points) - length(boundary)
        ring = rings{index_row};
        accumulate = 0;
        prime_point = points(:, index_row);
        
        for index_ring = 1 : length(ring(1,:))
            middle_point_index = ring(2,index_ring);
            alpha_point_index = ring(3,index_ring);
            if index_ring == 1
                a = length(ring(1,:));
            else
                a = index_ring-1;
            end
            if index_ring == length(ring(1,:))
                b = 1;
            else
                b = index_ring+1;
            end
            beta_point_index = ring(2,a);
            
            middle_point = points(:,middle_point_index);
            alpha_point = points(:,alpha_point_index);
            beta_point = points(:,beta_point_index);
            
            isRing = ring(2,1) == ring(3,length(ring(1,:)));
            
            gamma = get_angle(middle_point, prime_point, beta_point);
            delta = get_angle(middle_point, prime_point, alpha_point);
            
            h = (prime_point(1) - middle_point(1))^2 + (prime_point(2) - middle_point(2))^2 + (prime_point(3) - middle_point(3))^2; 
            
            if ring(3,index_ring) ~= ring(2,b)
                cod = 0;
            else
                cod = cot(delta);
            end
            if ring(2,index_ring) ~= ring(3,a)
                cog = 0;
            else
                cog = cot(gamma);
            end
            value = (cog + cod)/h;
            M_chi(index_row, middle_point_index) = value;
            accumulate = accumulate + value;
        end
            M_chi(index_row, index_row) = - accumulate;
    end
end

%Returns a sparse matrix representation of <triangles>
function matrix = get_sparse_triangle_matrix(points, triangles)
    matrix = sparse(length(points),length(triangles));
    a = [];
    for triangle_index = 1 : length(triangles)
        matrix(triangles(1,triangle_index), triangle_index) = 1;
        matrix(triangles(2,triangle_index), triangle_index) = 1;
        matrix(triangles(3,triangle_index), triangle_index) = 1;
        if ismember(2009,triangles(:,triangle_index))
            a = [a triangles(:,triangle_index)];
        end
    end
end

%Returns triangle indexes in <triangles_matrix> connected to <point_index>
function connected_triangles = get_connected_triangles(point_index, triangle_matrix)
    connected_triangles = [];
    searchVector = zeros(length(triangle_matrix(:,1)),1);
    searchVector(point_index) = 1;
    triangles = searchVector' * triangle_matrix;
    for triangle_index = 1 : length(triangles)
        if triangles(triangle_index) == 1
            connected_triangles = [connected_triangles; triangle_index];
        end
    end
    
end

% if <point_index> is in a triangle in <triangles> swaps position of
% <point_index> with value at <position>
function organized_triangles = organize_triangles(triangles, point_index, position)
    organized_triangles = [];
    for index = 1 : length(triangles(1,:))
        triangle = triangles(:, index);
        organized_triangles = [organized_triangles organize_triangle(triangle, point_index, position)];
    end
end

% if <point_index> is in <triangle> swaps position of <point_index> with 
% value at <position>
function organized_triangle = organize_triangle(triangle, point_index, position)
    if ismember(point_index, triangle)
        point_tri_index = find(triangle == point_index);
        triangle(point_tri_index) = triangle(position);
        triangle(position) = point_index;
        organized_triangle = triangle;
    end
end

function rings = form_triangle_rings(points, triangles)
    rings = [];
    triangle_mat = get_sparse_triangle_matrix(points, triangles);
    for i= 1:length(points)
        connected_triangles = get_connected_triangles(i,triangle_mat); 
        explicit_triangles = zeros(length(connected_triangles),3);
        for index = 1 : length(connected_triangles)
            explicit_triangles(index,:) = triangles(:,connected_triangles(index));
        end
        rings{i} = form_triangle_ring(explicit_triangles', i);
    end
end

% Forms a organized Ring out of a triangle vector <triangles>
% Requires the triangles to acutally form a ring
% <prime_point_index> is the common middle point
function ring = form_triangle_ring(triangles, prime_point_index)
    triangles = organize_triangles(triangles, prime_point_index, 1);
    ring = triangles(:,1);
    unsorted = triangles(:,2:length(triangles(1,:)));
    while ~isempty(unsorted)
        toSort = unsorted(:,1);
        new_ring = sort_triangle_into_ring(ring, toSort);
        if ~isempty(new_ring)
            if length(new_ring) <= 1
                unsorted = [];
            else
                unsorted = unsorted(:,2:length(unsorted(1,:)));
            end
            ring = new_ring;
        else
            unsorted = [unsorted(:,2:length(unsorted(1,:))) unsorted(:,1)];
        end
    end
end

%Sorts one common point first triangle <triangle> into a common point first triangle
%ring <ring>. The second index will be the one connecting left and the
%thrid the one connecting right
%This code is really suboptimal, please fix
function ring = sort_triangle_into_ring(ring, triangle)
    left = ring(:,1);
    right = ring(:,length(ring(1,:)));
    
    
        
    if ismember(triangle(2),left)
        if length(ring(1,:)) == 1
            if left(3) == triangle(2)
                ring = organize_triangle(left, left(3),2);
            end
        end
        ring = [organize_triangle(triangle, triangle(2), 3) ring];
    elseif ismember(triangle(3),left)
        if length(ring(1,:)) == 1
            if left(3) == triangle(3)
                ring = organize_triangle(left, left(3), 2);
            end
        end
        ring = [organize_triangle(triangle, triangle(3), 3) ring];
    elseif ismember(triangle(2),right)
        if length(ring(1,:)) == 1
            if right(2) == triangle(2)
                ring = organize_triangle(right, right(2),3);
            end
        end
        ring = [ring organize_triangle(triangle, triangle(2), 2)];
    elseif ismember(triangle(3),right)
        if length(ring(1,:)) == 1
            if right(2) == triangle(3)
                ring = organize_triangle(left, left(2),3);
            end
        end
        ring = [ring organize_triangle(triangle, triangle(3), 2)];
    else
        ring = [];
    end
end

% Returns the angle at <angle_point> between the vector spanned with
% <second_point> and <third_point>
% Formular used is inverse cosine on scalar product
function angle = get_angle(angle_point, second_point, third_point)
    from_second = (second_point - angle_point);
    from_second = from_second / sqrt(from_second(1)^2 + from_second(2)^2 + from_second(3)^2);
    from_third = (third_point - angle_point);
    from_third = from_third / sqrt(from_third(1)^2 + from_third(2)^2 + from_third(3)^2);
    a = dot(from_second, from_third);
    angle = acos(a);
end

% Rearanges the points such that the points corresponding to the
% <boundary> loop are at the very bottom. Updates the triangles accordingly
function [points, triangles] = rearange_for_boundary(boundary, points, triangles)
    for index = 1: length(boundary)
        [points,triangles] = swap_points(boundary(index), length(points)-length(boundary)+index, points, triangles);
        for index_rec = index+1: length(boundary)
            if boundary(index_rec) == boundary(index)
                boundary(index_rec) = length(points)-length(boundary)+index;
            else
                if boundary(index_rec) == length(points)-length(boundary)+index
                    boundary(index_rec) = boundary(index);
                end
            end
            
        end
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
function plot_mesh_collected(points, triangles, face_colors)
    points=points'; %flip to colums for convenience in matrix multiplication
    
    options = {};
    options.edge_color = [255,255, 255];
    options.face_vertex_color = face_colors';
    options.face_color = [0,0,0];
    figure, plotmesh(points,triangles', options);
    rzview('on') %allows for rotating zoom and move with the mouse (left-right-middle)

end

%Samples <amount> number of points on a unit circle
function samples = sample_unit_circle(amount)
    samples = [];
    for i = 0 : amount-1
        samples = [samples; cos(2*pi /amount*i - pi) sin(2*pi /amount*i - pi)];
    end 
end
