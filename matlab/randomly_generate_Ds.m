function points2 = randomly_generate_Ds(LU,B)

points = [];
        n = size(LU,1);
        upp2 = max(LU(1).Lat);
low2 = min(LU(1).Lat);
upp1 = max(LU(1).Lon);
low1 = min(LU(1).Lon);
for j = 1:n
        upp1=max([upp1,LU(j).Lon]);
        low1=min([low1,LU(j).Lon]);
        upp2=max([upp2,LU(j).Lat]);
        low2=min([low2,LU(j).Lat]);

end
        
        s1=unifrnd(low1,upp1,[B,1]);
        s2=unifrnd(low2,upp2,[B,1]);
for j = 1:n
        check=inpolygon(s1,s2,LU(j).Lon,LU(j).Lat);
        points = cat(1,points,cat(2,s1(check==1),s2(check==1)));
end
while size(points,1)<B
        s1=unifrnd(low1,upp1,[B,1]);
        s2=unifrnd(low2,upp2,[B,1]);
for j = 1:n
        check=inpolygon(s1,s2,LU(j).Lon,LU(j).Lat);
        points = cat(1,points,cat(2,s1(check==1),s2(check==1)));
end
end
points2(:,2) = points(1:B,2);
points2(:,1) = points(1:B,1);

end