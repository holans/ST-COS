function S = ArealBi2_spacetime(LU,times,level1,level2,level3,B,srad,trad)
r = size(level1,1) + size(level2,1) + size(level3,1);
n = size(LU,1);
S = zeros(n,r);
T = max(size(times));

for t = 1:T
    for j = 1:n
        points2 = randomly_generate_Ds(LU(j),B);
        s2{t}(j,:)=points2(:,1)';
        s1{t}(j,:)=points2(:,2)';
    end
end

for t = 1:T
    for i = 1:B
	    Sq = Create_S_time([s2{t}(:,i),s1{t}(:,i)],times(t),level1,level2,level3,srad,trad);
		S = Sq + S;
        if mod(i,1)==0
			disp(strcat('***Iteration ',num2str(i))); 
        end
    end
end
S = S/(B*T);


end
