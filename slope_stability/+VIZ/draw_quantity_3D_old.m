function draw_quantity_3D_old(coord,surf,U,Q_node,x1,x2,x3,y1,y2,z)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord     - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%                number of nodes
%    surf      - array containing numbers of nodes defining each surface
%                element, size(surf)=(n_p,n_s), n_s = number of surface elements
%    U         - nodal displacements, size(U)=(3,n_n) to catch deformed shape
%                if the deformed shape is not required then set 0*U
%    Q_node    - prescribed nodal quantity, size(Q_node)=(1,n_n)
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%    size_xy   - size of the body in x and y direction (integer)
%    size_z    - size of the body in z-direction (integer) 
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%    body=(0,size_xy)  x(0,size_xy)  x(0,size_z)\setminus
%         (0,size_hole)x(0,size_hole)x(0,size_z)
%
% ======================================================================
%

  figure;
  hold on;
  
  % visualization of the quantity
  patch('Faces',surf(1:3,:)','Vertices',coord'+U',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
  colorbar;
  
  % undeformed shape of the body
  plot3([0,x1+x2+x3],[0,0],[0,0])
  plot3([x1+x2,x1+x2+x3],[y1,y1],[0,0])
  plot3([0,x1],[y1+y2,y1+y2],[0,0])
  plot3([0,x1+x2+x3],[0,0],[z,z])
  plot3([x1+x2,x1+x2+x3],[y1,y1],[z,z])
  plot3([0,x1],[y1+y2,y1+y2],[z,z])
  plot3([x1,x1+x2],[y1+y2,y1],[0,0])
  plot3([x1,x1+x2],[y1+y2,y1],[z,z])
  plot3([0,0],[0,y1+y2],[0,0])
  plot3([0,0],[0,y1+y2],[z,z])
  plot3([x1+x2+x3,x1+x2+x3],[0,y1],[0,0])
  plot3([x1+x2+x3,x1+x2+x3],[0,y1],[z,z])
  plot3([0,0],[0,0],[0,z])
  plot3([x1+x2+x3,x1+x2+x3],[0,0],[0,z])
  plot3([0,0],[y1+y2,y1+y2],[0,z])
  plot3([x1,x1],[y1+y2,y1+y2],[0,z])
  plot3([x1+x2,x1+x2],[y1,y1],[0,z])
  plot3([x1+x2+x3,x1+x2+x3],[y1,y1],[0,z]) 
  
  %
  box on
  view([0.5 1 -2]);     
  axis equal;  % real ratios
  hold off;
 % axis off;
end