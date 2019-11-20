%Here are two examples of coloring curves (this is a hidden matlab functionality)
% use it in combination with your colormap for showing curvature.
t=pi:.01:2*pi;
mx=2.5+cos(t);
my=3.5+sin(t); 

figure
ha = plot(mx,my, '.-', 'LineWidth',10, 'DisplayName',' 0.5');

%example 1
n=length(t);
cd= uint8(randi(255,4,n));
% colors are stored in the first three rows, each column corresponds to a 
% the color of a sampling point on the curves [r g b]'
cd(:,4)=0;
drawnow
set(ha.Edge, 'ColorBinding','interpolated', 'ColorData',cd)


%example 2 using the parula colormap from matlab
figure
ha = plot(mx,my, '.-', 'LineWidth',10, 'DisplayName',' 0.5');

n=length(t);
cd= uint8(randi(255,4,n));
cd = [uint8(parula(n)*255) uint8(ones(n,1))].'; %'
cd(:,4)=0;
drawnow
set(ha.Edge, 'ColorBinding','interpolated', 'ColorData',cd)