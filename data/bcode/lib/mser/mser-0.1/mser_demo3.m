% MSER_DEMO3  Demonstrates MSER on a volumetric image

% --------------------------------------------------------------------
%                                                          Create data
% --------------------------------------------------------------------

% volumetric coordinate (x,y,z)
x       = linspace(-1,1,50) ;
[x,y,z] = meshgrid(x,x,x) ;

% create funny volumetric image
I   = sin(4*x).*cos(4*y).*sin(z) ;
I   = I-min(I(:)) ;
I   = I/max(I(:)) ;

% quantize the image in 10 levels
lev = 10 ;
I   = lev*I ;
Ir  = round(I) ;

% --------------------------------------------------------------------
%                                                      Compute regions
% --------------------------------------------------------------------
[idx,ell,p] = mser(uint8(Ir),uint8(1));

% --------------------------------------------------------------------
%                                                                Plots
% --------------------------------------------------------------------

% The image is quantized; store in LEV its range.
lev = unique(Ir(idx)) ;

figure(100); clf;
K=min(length(lev),4) ;

r=.99 ;

% one level per time
for k=1:K
  tightsubplot(K,k) ;
  [i,j,m] = ind2sub(size(I), idx(Ir(idx)==lev(k)) ) ;
  
  % compute level set of level LEV(k)
  Is = double(Ir<=lev(k)) ;
    
  p1 = patch(isosurface(Is,r), ...
             'FaceColor','blue','EdgeColor','none') ;
  p2 = patch(isocaps(Is,r),...
             'FaceColor','interp','EdgeColor','none') ;  
  isonormals(I,p1)
  hold on ;
%  plot3(j,i,m,'g*','MarkerSize',50) ;
  
%  [i,j,m]=ind2sub(size(I), find(Ir<=lev(k))) ;
%  plot3(j,i,m,'b.') ;
  
  view(3); axis vis3d tight
  camlight; lighting phong ;
  
  % find regions that have this level
  sel = find( Ir(idx) == lev(k) ) ;
  
  % plot fitted ellipsoid
  for r=sel'
    E = ell(:,r) ;
    c = E(1:3) ;
    A = zeros(3) ;
    A(1,1) = E(4) ;
    A(1,2) = E(5) ;
    A(2,2) = E(6) ;
    A(1,3) = E(7) ;
    A(2,3) = E(8) ;
    A(3,3) = E(9) ;
    
    A = A + A' - diag(diag(A)) ;

    % correct var. order
    perm = [0 1 0 ; 1 0 0 ; 0 0 1] ;
    A = perm*A*perm ;
    
    [V,D] = eig(A) ;
    A = 2.5*V*sqrt(D) ;
    
    [x,y,z]=sphere ;
    [P,Q]=size(x) ;    
    X=A*[x(:)';y(:)';z(:)'] ;
    x=reshape(X(1,:),P,Q)+c(2) ;
    y=reshape(X(2,:),P,Q)+c(1) ;
    z=reshape(X(3,:),P,Q)+c(3) ;
    surf(x,y,z,'FaceAlpha',.5) ;
    
%    plot3(c(2),c(1),c(3),'r*','MarkerSize',70) ;

    
  end
end
