function [Darea,Garea,DGoverlap] = MCAREA(dx,dy,gx,gy,n_mc)
%*****************************************************
%* Calculates area of overlap of two polygons
%*  using a Monte Carlo approach; note - if
%*  there is no overlap, the algorithm doesn't
%*  calculate the areas for the two polygons (rather
%*  it returns 0 for those)
%* 
%*  dx,dy - x and y coordinates of "data" polygon
%*  gx,gy - x and y coordinates of "grid" polygon
%*  n_mc - number of Monte Carlo samples
%*
%****************************************************

dxmin = min(dx);
dxmax = max(dx);
dymin = min(dy);
dymax = max(dy);

gxmin = min(gx);
gxmax = max(gx);
gymin = min(gy);
gymax = max(gy);


%Do a quick check for overlap; if bounding squares don't
% intersect, then polygons can't intersect
dsx = [dxmin dxmin dxmax dxmax dxmin];
dsy = [dymin dymax dymax dymin dymin];
gsx = [gxmin gxmin gxmax gxmax gxmin];
gsy = [gymin gymax gymax gymin gymin];

tstint = polybool('intersection',dsx,dsy,gsx,gsy);
if isempty(tstint)
    DGoverlap = 0;
    Darea = 0;
    Garea = 0;
else

   %area of overlap
   xmin = min([dxmin,gxmin]);
   xmax = max([dxmax,gxmax]);
   ymin = min([dymin,gymin]);
   ymax = max([dymax,gymax]);
   sq_ar = abs(xmax-xmin) * abs(ymax-ymin);
   xsamp = unifrnd(xmin,xmax,n_mc,1);
   ysamp = unifrnd(ymin,ymax,n_mc,1);
   din2 = inpolygon(xsamp,ysamp,dx,dy);
   gin2 = inpolygon(xsamp,ysamp,gx,gy);
   dgprod = din2 .* gin2;
   DGoverlap = (sum(dgprod)/n_mc)*sq_ar;

   if DGoverlap ~= 0
      %area of d polygon
      xsamp = unifrnd(dxmin,dxmax,n_mc,1);
      ysamp = unifrnd(dymin,dymax,n_mc,1);
      dsq_ar = abs(dxmax-dxmin) * abs(dymax-dymin);
      din = inpolygon(xsamp,ysamp,dx,dy);
      Darea = (sum(din)/n_mc)*dsq_ar;

      %area of g polygon
      xsamp = unifrnd(gxmin,gxmax,n_mc,1);
      ysamp = unifrnd(gymin,gymax,n_mc,1);
      gsq_ar = abs(gxmax-gxmin) * abs(gymax-gymin);
      gin = inpolygon(xsamp,ysamp,gx,gy);
      Garea = (sum(gin)/n_mc)*gsq_ar;
   else
       Darea = 0;
       Garea = 0;
   end
   end
end
