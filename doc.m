%% discrete_contour_mesh_patch
%
% Function to mesh one discrete 2D or 3D contour
% composed of sorted or disordered 3D points.
%
% Author and support : nicolas.douillet (at) free.fr, 2020.
%
%% Syntax
%
% [C, T] = discrete_contour_mesh_patch(V);
% [C, T] = discrete_contour_mesh_patch(V, mode);
% [C, T, N] = discrete_contour_mesh_patch(V);
%
%% Description
%
% [C, T] = discrete_contour_mesh_patch(V) meshes the discrete contour
% defined by the vertex set V and store the resulting triangulation in T
% after having sorted -default behaviour- V elements by their angle
% relatively to their isobarycentre, and whom result is stored in C.
%
% [C, T] = discrete_contour_mesh_patch(V, mode) uses the mode parameter to
% sort (mode = 'raw') or not (mode = 'sorted') the vertex set V.
%
% [C, T, N] = discrete_contour_mesh_patch(V) also computes and returns the
% vertex normals set N.
%
%% See also
%
% <https://fr.mathworks.com/help/matlab/ref/mesh.html?s_tid=srchtitle mesh> |
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html?s_tid=srchtitle trimesh> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/77004-mesh-processing-toolbox?s_tid=prof_contriblnk mesh_processing_toolbox>
%
%% Input arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the point set, size(V) = [nb_vertices,3].
%        [ |  | | ]
%
% - mode : character string in the set {'raw','sorted'},
%          the type of contour considered. Case insensitive.
%
%% Output arguments
%
%        [ |  |  |]
% - C = [Cx Cy Cz], real matrix double, the output sorted point set, size(C) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the resulting triangulation, size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - N : [Nx Ny Nz], real matrix double, the vertex normal vectors, size(N) = [nb_vertices,3].
%        [ |  |  |]
%
%% Example #1 : convex hull of random 2D point set
n = 32;
V = 2*(rand(n,2)-0.5);
H_raw = convhull(V);
V = cat(2,V,zeros(size(V,1),1));        

Rmx = @(theta)[1 0          0;
               0 cos(theta) -sin(theta);
               0 sin(theta) cos(theta)];

V = (Rmx(0.25*pi)*V')'; % pi/4 rotation around X axe, to check real 3D
V = V(unique(H_raw,'stable'),:);
V = V([end,end-1,1,2:end-2],:); % disorder

[C,T,N] = discrete_contour_mesh_patch(V);

figure;
plot3(C(:,1),C(:,2),C(:,3),'bo','LineWidth',4,'MarkerSize',6,'MarkerFaceColor', [0 0 1],'MarkerEdgeColor', [0 0 1]), hold on;
line(cat(1,C(:,1),C(1,1)),cat(1,C(:,2),C(1,2)),cat(1,C(:,3),C(1,3)),'Color',[0 1 0],'LineWidth',6), hold on;
quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3),'Color',[1 0.5 0],'Linewidth',2), hold on;
trisurf(T,C(:,1),C(:,2),C(:,3)), shading faceted, hold on;
colormap([0 1 1]);        
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal, axis tight;
set(gcf,'Color',[0 0 0]), set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);

%% Example #2 : 2D discrete ellipse
n = 16;
angl_step = 2*pi/n;
alpha = 0:angl_step:2*pi-angl_step;
X = 0.5*cos(alpha)';
Y = 2*sin(alpha)';
V = cat(2,X,Y);
H_raw = convhull(V);
V = cat(2,V,zeros(size(V,1),1));        

Rmx = @(theta)[1 0          0;
               0 cos(theta) -sin(theta);
               0 sin(theta) cos(theta)];

V = (Rmx(0.25*pi)*V')'; % pi/4 rotation around X axe, to check real 3D
V = V(unique(H_raw,'stable'),:);
V = V([end,end-1,1,2:end-2],:); % disorder

[C,T,N] = discrete_contour_mesh_patch(V);

figure;
plot3(C(:,1),C(:,2),C(:,3),'bo','LineWidth',4,'MarkerSize',6,'MarkerFaceColor', [0 0 1],'MarkerEdgeColor', [0 0 1]), hold on;
line(cat(1,C(:,1),C(1,1)),cat(1,C(:,2),C(1,2)),cat(1,C(:,3),C(1,3)),'Color',[0 1 0],'LineWidth',6), hold on;
quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3),'Color',[1 0.5 0],'Linewidth',2), hold on;
trisurf(T,C(:,1),C(:,2),C(:,3)), shading faceted, hold on;
colormap([0 1 1]);        
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal, axis tight;
set(gcf,'Color',[0 0 0]), set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);

%% Example #3 : 3D 'flower' test
n = 16;
angl_step = 2*pi/n;
alpha = 0:angl_step:2*pi-angl_step;
X = cos(alpha)';
Y = sin(alpha)';
R = sqrt(X.^2 + Y.^2);
X = X + (0.2*R.*sin(3*alpha)'.*cos(3*alpha)');
Y = Y + (0.2*R.*sin(3*alpha)'.*sin(3*alpha)');
Z = sqrt(X.^2 + Y.^2);
V = cat(2,X,Y,Z);

V = V([end,end-1,1,2:end-2],:); % disorder

[C,T,N] = discrete_contour_mesh_patch(V);

figure;
plot3(C(:,1),C(:,2),C(:,3),'bo','LineWidth',4,'MarkerSize',6,'MarkerFaceColor', [0 0 1],'MarkerEdgeColor', [0 0 1]), hold on;
line(cat(1,C(:,1),C(1,1)),cat(1,C(:,2),C(1,2)),cat(1,C(:,3),C(1,3)),'Color',[0 1 0],'LineWidth',6), hold on;
quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3),'Color',[1 0.5 0],'Linewidth',2), hold on;
trisurf(T,C(:,1),C(:,2),C(:,3)), shading faceted, hold on;
colormap([0 1 1]);        
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal, axis tight;
set(gcf,'Color',[0 0 0]), set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);