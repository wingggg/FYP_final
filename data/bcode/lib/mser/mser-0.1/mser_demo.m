function mser_demo
% MSER_DEMO  Demonstrates MSER 

which_image = 3 ;

% --------------------------------------------------------------------
%                                                          Create data
% --------------------------------------------------------------------
switch which_image
  case 1
    I   = rand(200,200) ;
    I   = imsmooth(I,10) ;
    I   = I-min(I(:)) ;
    I   = I/max(I(:)) ;
    lev = 10 ;
    I   = uint8(round(I*lev)) ;
  
  case 2    
    I = zeros(200,200) ;
    I(50:150,50:150)=5 ;
    I = imsmooth(I,10) ;
    I = uint8(round(I)) ;
    
  case 3
    I = imreadbw('data/car1.jpg') ;
    I = uint8(255*I) ;
end

% --------------------------------------------------------------------
%                                                        Compute MSERs
% --------------------------------------------------------------------
[idx,ell,p,a] = mser(I, uint8(2)) ;

% --------------------------------------------------------------------
%                                                                Plots
% --------------------------------------------------------------------
[i,j] = ind2sub(size(I),idx) ;

figure(100) ; clf ; imagesc(I) ; hold on ;
set(gca,'Position',[0 0 1 1]) ;
plot(j,i,'g*') ; colormap gray ;

% swap x with y
ell = ell([2 1 5 4 3],:) ;

for k=1:size(ell,2)
  E = ell(:,k) ;
  c = E(1:2) ;
  A = zeros(2) ;
  A(1,1) = E(3) ;
  A(1,2) = E(4) ;
  A(2,2) = E(5) ;  
  A = A + A' - diag(diag(A)) ;
  
  [V,D] = eig(A) ;
  A = 2.5*V*sqrt(D) ;
  
  X = A*[cos(linspace(0,2*pi,30)) ; sin(linspace(0,2*pi,30)) ;] ;
  X(1,:) = X(1,:) + c(1) ;
  X(2,:) = X(2,:) + c(2) ;

  plot(X(1,:),X(2,:),'r-','LineWidth',2) ;
  plot(c(1),c(2),'r.') ;
  plot(j(k),i(k),'g*') ;
end

line([j'; ell(1,:)],[i'; ell(2,:)],'color','b') ;
