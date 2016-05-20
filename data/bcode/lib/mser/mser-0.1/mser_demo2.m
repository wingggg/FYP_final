% MSER_DEMO2  Demonstrate MSER code

I = load('clown') ; I = uint8(I.X) ;
figure(1) ; imagesc(I) ; colormap gray; hold on ;

i = double(i) ;
j = double(j) ;

[r,ell] = mser(I,uint8(5)) ;

r=double(r) ;

[i,j]=ind2sub(size(I),r) ;                    
plot(j,i,'r*') ;

ell = ell([2 1 5 4 3],:) ;
plotframe(ell); 
return


figure(2) ; 

clear MOV ;
K = size(ell,2) ;
for k=1:K
  clf ;
  Ib = I <= I(i(k),j(k)) ;
  mask = bwselect(Ib,j(k),i(k),8) ;
  imagesc(cat(3,I,255*uint8(mask),I)) ; colormap gray ; hold on ;
  %set(gca,'position',[0 0 1 1]) ; axis off ; axis equal ;
  plot(j(k),i(k),'r*') ;
  plotframe(ell(:,k),'color','r') ;
  MOV(k) = getframe(gca) ;
end
