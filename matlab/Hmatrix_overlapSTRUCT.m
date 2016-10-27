function [H] = Hmatrix_overlapSTRUCT(D,G,denomflag,Dxyflag,Gxyflag,n_mc)
%HMATRIX_OVERLAP.M 
% calculate area of overlap of data polygons and grid polygons
% to form a change of support matrix, H
%
%   D - structural array of length nd that has .X and .Y locations for 
%        data polygons
%   
%   G - structural array of length ng that has .X and .Y locations for 
%         grid polygons [NOTE: for school districts, use .LON and .LAT 
%         instead of .X and .Y
%  
%   denomflag: 1 - denom = 1
%              2 - D area in the denom
%              3 - G area in the denom
%
%   Dxyflag: 1 - use .X and .Y, else use .LON and .LAT
%   Gxyflag: 1 - use .X and .Y, else use .LON and .LAT
%
%   n_mc: number of Monte Carlo samples
%
%  modification of 7/25/2012 version for structural
%    arrays; 01/09/2013 ckw
%*******************************************

[nd,tmp1] = size(D);
[ng,tmp2] = size(G);

H = sparse(zeros(nd,ng));
for k=1:nd
    if mod(k,1)==0
        disp(sprintf('***Location: %d',k))
    end
    
    if Dxyflag == 1
       dx = D(k).X;
       dy = D(k).Y;
    else
       dx = D(k).LON;
       dy = D(k).LAT;
    end
    for j=1:ng
        if Gxyflag == 1
            gx = G(j).X;
            gy = G(j).Y;
        else
            gx = G(j).LON;
            gy = G(j).LAT;
        end
        [Darea,Garea,DGoverlap]=MCAREA(dx,dy,gx,gy,n_mc);
        if denomflag == 1   %1 in denom
            H(k,j) = DGoverlap;
        elseif denomflag == 2  %D area in denom
            if DGoverlap == 0
                H(k,j) = 0;
            else
               tmp = DGoverlap/Darea;
               if tmp > 1
                   tmp = 1;
               end
               H(k,j) = tmp;
            end
        elseif denomflag == 3   %G area in denom
            if DGoverlap == 0
                H(k,j) = 0;
            else
               tmp = DGoverlap/Garea;
               if tmp > 1
                  tmp = 1;
               end
               H(k,j) = tmp;
            end
        end
    end
end


end

