
clear all; 

%%%%%%%%%%%%%%%%2013 period 1%%%%%%%%%%%%%%%%%%%%
county1 = shaperead('period1');
names2 = {county1.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1 = county1(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1,1)
    [lat lon] = minvtran(mstruct,county1(k).X,county1(k).Y);
    county1(k).Lat = lat;
    county1(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1.csv');
% Zagg = sqrt(county(:,1));
eZagg1 = countyd1(indexConus2,1);
%known survey variances
esigmavar1 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1 = log(eZagg1);
sigmavar1 = esigmavar1.*(1./eZagg1).^2;



%%%%%%%%%%%%%%%%2013 period 2%%%%%%%%%%%%%%%%%%%%
county2 = shaperead('period2');
names2 = {county2.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2 = county2(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2,1)
    [lat lon] = minvtran(mstruct,county2(k).X,county2(k).Y);
    county2(k).Lat = lat;
    county2(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2.csv');
% Zagg = sqrt(county(:,1));
eZagg2 = countyd2(indexConus2,1);
%known survey variances
esigmavar2 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2 = log(eZagg2);
sigmavar2 = esigmavar2.*(1./eZagg2).^2;

%%%%%%%%%%%%%%%%2013 period 3%%%%%%%%%%%%%%%%%%%%
county3 = shaperead('period3');
names2 = {county3.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county3);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county3 = county3(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county3,1)
    [lat lon] = minvtran(mstruct,county3(k).X,county3(k).Y);
    county3(k).Lat = lat;
    county3(k).Lon = lon;
 end
 

%data
countyd3=dlmread('period3.csv');
% Zagg = sqrt(county(:,1));
eZagg3 = countyd3(indexConus2,1);
%known survey variances
esigmavar3 = (countyd3(indexConus2,2)/1.645).^2;
Zagg3 = log(eZagg3);
sigmavar3 = esigmavar3.*(1./eZagg3).^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county2;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2 = H1;
 for k=1:ng
     tsum = sum(H2(:,k));
     if tsum ~= 0
        H2(:,k) = H2(:,k)/tsum;
     end
 end
 
 G = county1;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H);
 H1 = H;
 for k=1:ng
     tsum = sum(H1(:,k));
     if tsum ~= 0
        H1(:,k) = H1(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%2012 period 1%%%%%%%%%%%%%%%%%%%%
county1_2012 = shaperead('period1_2012');
names2 = {county1_2012.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2012);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2012 = county1_2012(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2012,1)
    [lat lon] = minvtran(mstruct,county1_2012(k).X,county1_2012(k).Y);
    county1_2012(k).Lat = lat;
    county1_2012(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2012.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2012 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2012 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2012 = log(eZagg1_2012);
sigmavar1_2012 = esigmavar1_2012.*(1./eZagg1_2012).^2;


%%%%%%%%%%%%%%%%2012 period 2%%%%%%%%%%%%%%%%%%%%
county2_2012 = shaperead('period2_2012');
names2 = {county2_2012.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2012);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2012 = county2_2012(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2012,1)
    [lat lon] = minvtran(mstruct,county2_2012(k).X,county2_2012(k).Y);
    county2_2012(k).Lat = lat;
    county2_2012(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2012.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2012 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2012 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2012 = log(eZagg2_2012);
sigmavar2_2012 = esigmavar2_2012.*(1./eZagg2_2012).^2;

%%%%%%%%%%%%%%%%2012 period 3%%%%%%%%%%%%%%%%%%%%
county3_2012 = shaperead('period3_2012');
names2 = {county3_2012.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county3_2012);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county3_2012 = county3_2012(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county3_2012,1)
    [lat lon] = minvtran(mstruct,county3_2012(k).X,county3_2012(k).Y);
    county3_2012(k).Lat = lat;
    county3_2012(k).Lon = lon;
 end
 

%data
countyd3=dlmread('period3_2012.csv');
% Zagg = sqrt(county(:,1));
eZagg3_2012 = countyd3(indexConus2,1);
%known survey variances
esigmavar3_2012 = (countyd3(indexConus2,2)/1.645).^2;
Zagg3_2012 = log(eZagg3_2012);
sigmavar3_2012 = esigmavar3_2012.*(1./eZagg3_2012).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county3_2012;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H3_2012 = H1;
 for k=1:ng
     tsum = sum(H3_2012(:,k));
     if tsum ~= 0
        H3_2012(:,k) = H3_2012(:,k)/tsum;
     end
 end

 G = county2_2012;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2012 = H1;
 for k=1:ng
     tsum = sum(H2_2012(:,k));
     if tsum ~= 0
        H2_2012(:,k) = H2_2012(:,k)/tsum;
     end
 end
 
 G = county1_2012;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H);
 H1_2012 = H;
 for k=1:ng
     tsum = sum(H1_2012(:,k));
     if tsum ~= 0
        H1_2012(:,k) = H1_2012(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%2011 period 1%%%%%%%%%%%%%%%%%%%%
county1_2011 = shaperead('period1_2011');
names2 = {county1_2011.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2011);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2011 = county1_2011(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2011,1)
    [lat lon] = minvtran(mstruct,county1_2011(k).X,county1_2011(k).Y);
    county1_2011(k).Lat = lat;
    county1_2011(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2011.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2011 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2011 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2011 = log(eZagg1_2011);
sigmavar1_2011 = esigmavar1_2011.*(1./eZagg1_2011).^2;



%%%%%%%%%%%%%%%%2011 period 2%%%%%%%%%%%%%%%%%%%%
county2_2011 = shaperead('period2_2011');
names2 = {county2_2011.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2011);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2011 = county2_2011(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2011,1)
    [lat lon] = minvtran(mstruct,county2_2011(k).X,county2_2011(k).Y);
    county2_2011(k).Lat = lat;
    county2_2011(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2011.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2011 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2011 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2011 = log(eZagg2_2011);
sigmavar2_2011 = esigmavar2_2011.*(1./eZagg2_2011).^2;

%%%%%%%%%%%%%%%%2011 period 3%%%%%%%%%%%%%%%%%%%%
county3_2011 = shaperead('period3_2011');
names2 = {county3_2011.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county3_2011);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county3_2011 = county3_2011(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county3_2011,1)
    [lat lon] = minvtran(mstruct,county3_2011(k).X,county3_2011(k).Y);
    county3_2011(k).Lat = lat;
    county3_2011(k).Lon = lon;
 end
 

%data
countyd3=dlmread('period3_2011.csv');
% Zagg = sqrt(county(:,1));
eZagg3_2011 = countyd3(indexConus2,1);
%known survey variances
esigmavar3_2011 = (countyd3(indexConus2,2)/1.645).^2;
Zagg3_2011 = log(eZagg3_2011);
sigmavar3_2011 = esigmavar3_2011.*(1./eZagg3_2011).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county3_2011;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H3_2011 = H1;
 for k=1:ng
     tsum = sum(H3_2011(:,k));
     if tsum ~= 0
        H3_2011(:,k) = H3_2011(:,k)/tsum;
     end
 end

 G = county2_2011;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2011 = H1;
 for k=1:ng
     tsum = sum(H2_2011(:,k));
     if tsum ~= 0
        H2_2011(:,k) = H2_2011(:,k)/tsum;
     end
 end
 
 G = county1_2011;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H);
 H1_2011 = H;
 for k=1:ng
     tsum = sum(H1_2011(:,k));
     if tsum ~= 0
        H1_2011(:,k) = H1_2011(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%2010 period 1%%%%%%%%%%%%%%%%%%%%
county1_2010 = shaperead('period1_2010');
names2 = {county1_2010.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2010);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2010 = county1_2010(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2010,1)
    [lat lon] = minvtran(mstruct,county1_2010(k).X,county1_2010(k).Y);
    county1_2010(k).Lat = lat;
    county1_2010(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2010.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2010 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2010 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2010 = log(eZagg1_2010);
sigmavar1_2010 = esigmavar1_2010.*(1./eZagg1_2010).^2;



%%%%%%%%%%%%%%%%2010 period 2%%%%%%%%%%%%%%%%%%%%
county2_2010 = shaperead('period2_2010');
names2 = {county2_2010.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2010);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2010 = county2_2010(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2010,1)
    [lat lon] = minvtran(mstruct,county2_2010(k).X,county2_2010(k).Y);
    county2_2010(k).Lat = lat;
    county2_2010(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2010.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2010 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2010 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2010 = log(eZagg2_2010);
sigmavar2_2010 = esigmavar2_2010.*(1./eZagg2_2010).^2;

%%%%%%%%%%%%%%%%2010 period 3%%%%%%%%%%%%%%%%%%%%
county3_2010 = shaperead('period3_2010');
names2 = {county3_2010.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county3_2010);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county3_2010 = county3_2010(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county3_2010,1)
    [lat lon] = minvtran(mstruct,county3_2010(k).X,county3_2010(k).Y);
    county3_2010(k).Lat = lat;
    county3_2010(k).Lon = lon;
 end
 

%data
countyd3=dlmread('period3_2010.csv');
% Zagg = sqrt(county(:,1));
eZagg3_2010 = countyd3(indexConus2,1);
%known survey variances
esigmavar3_2010 = (countyd3(indexConus2,2)/1.645).^2;
Zagg3_2010 = log(eZagg3_2010);
sigmavar3_2010 = esigmavar3_2010.*(1./eZagg3_2010).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county3_2010;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H3_2010 = H1;
 for k=1:ng
     tsum = sum(H3_2010(:,k));
     if tsum ~= 0
        H3_2010(:,k) = H3_2010(:,k)/tsum;
     end
 end

 G = county2_2010;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2010 = H1;
 for k=1:ng
     tsum = sum(H2_2010(:,k));
     if tsum ~= 0
        H2_2010(:,k) = H2_2010(:,k)/tsum;
     end
 end
 
 G = county1_2010;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H);
 H1_2010 = H;
 for k=1:ng
     tsum = sum(H1_2010(:,k));
     if tsum ~= 0
        H1_2010(:,k) = H1_2010(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%2009 period 1%%%%%%%%%%%%%%%%%%%%
county1_2009 = shaperead('period1_2009');
names2 = {county1_2009.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2009);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2009 = county1_2009(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2009,1)
    [lat lon] = minvtran(mstruct,county1_2009(k).X,county1_2009(k).Y);
    county1_2009(k).Lat = lat;
    county1_2009(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2009.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2009 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2009 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2009 = log(eZagg1_2009);
sigmavar1_2009 = esigmavar1_2009.*(1./eZagg1_2009).^2;



%%%%%%%%%%%%%%%%2009 period 2%%%%%%%%%%%%%%%%%%%%
county2_2009 = shaperead('period2_2009');
names2 = {county2_2009.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2009);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2009 = county2_2009(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2009,1)
    [lat lon] = minvtran(mstruct,county2_2009(k).X,county2_2009(k).Y);
    county2_2009(k).Lat = lat;
    county2_2009(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2009.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2009 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2009 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2009 = log(eZagg2_2009);
sigmavar2_2009 = esigmavar2_2009.*(1./eZagg2_2009).^2;

%%%%%%%%%%%%%%%%2009 period 3%%%%%%%%%%%%%%%%%%%%
county3_2009 = shaperead('period3_2009');
names2 = {county3_2009.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county3_2009);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county3_2009 = county3_2009(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county3_2009,1)
    [lat lon] = minvtran(mstruct,county3_2009(k).X,county3_2009(k).Y);
    county3_2009(k).Lat = lat;
    county3_2009(k).Lon = lon;
 end
 

%data
countyd3=dlmread('period3_2009.csv');
% Zagg = sqrt(county(:,1));
eZagg3_2009 = countyd3(indexConus2,1);
%known survey variances
esigmavar3_2009 = (countyd3(indexConus2,2)/1.645).^2;
Zagg3_2009 = log(eZagg3_2009);
sigmavar3_2009 = esigmavar3_2009.*(1./eZagg3_2009).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county3_2009;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H3_2009 = H1;
 for k=1:ng
     tsum = sum(H3_2009(:,k));
     if tsum ~= 0
        H3_2009(:,k) = H3_2009(:,k)/tsum;
     end
 end

 G = county2_2009;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2009 = H1;
 for k=1:ng
     tsum = sum(H2_2009(:,k));
     if tsum ~= 0
        H2_2009(:,k) = H2_2009(:,k)/tsum;
     end
 end
 
 G = county1_2009;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H);
 H1_2009 = H;
 for k=1:ng
     tsum = sum(H1_2009(:,k));
     if tsum ~= 0
        H1_2009(:,k) = H1_2009(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%2006 period 1%%%%%%%%%%%%%%%%%%%%
county1_2006 = shaperead('period1_2006');
names2 = {county1_2006.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2006);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2006 = county1_2006(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2006,1)
    [lat lon] = minvtran(mstruct,county1_2006(k).X,county1_2006(k).Y);
    county1_2006(k).Lat = lat;
    county1_2006(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2006.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2006 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2006 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2006 = log(eZagg1_2006);
sigmavar1_2006 = esigmavar1_2006.*(1./eZagg1_2006).^2;

%%%%%%%%%%%%%%%%2007 period 1%%%%%%%%%%%%%%%%%%%%
county1_2007 = shaperead('period1_2007');
names2 = {county1_2007.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2007);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2007 = county1_2007(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2007,1)
    [lat lon] = minvtran(mstruct,county1_2007(k).X,county1_2007(k).Y);
    county1_2007(k).Lat = lat;
    county1_2007(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2007.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2007 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2007 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2007 = log(eZagg1_2007);
sigmavar1_2007 = esigmavar1_2007.*(1./eZagg1_2007).^2;

%%%%%%%%%%%%%%%%2008 period 1%%%%%%%%%%%%%%%%%%%%
county1_2008 = shaperead('period1_2008');
names2 = {county1_2008.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county1_2008);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county1_2008 = county1_2008(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county1_2008,1)
    [lat lon] = minvtran(mstruct,county1_2008(k).X,county1_2008(k).Y);
    county1_2008(k).Lat = lat;
    county1_2008(k).Lon = lon;
 end
 

%data
countyd1=dlmread('period1_2008.csv');
% Zagg = sqrt(county(:,1));
eZagg1_2008 = countyd1(indexConus2,1);
%known survey variances
esigmavar1_2008 = (countyd1(indexConus2,2)/1.645).^2;
Zagg1_2008 = log(eZagg1_2008);
sigmavar1_2008 = esigmavar1_2008.*(1./eZagg1_2008).^2;

%%%%%%%%%%%%%%%%2007 period 2%%%%%%%%%%%%%%%%%%%%
county2_2007 = shaperead('period2_2007');
names2 = {county2_2007.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2007);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2007 = county2_2007(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2007,1)
    [lat lon] = minvtran(mstruct,county2_2007(k).X,county2_2007(k).Y);
    county2_2007(k).Lat = lat;
    county2_2007(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2007.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2007 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2007 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2007 = log(eZagg2_2007);
sigmavar2_2007 = esigmavar2_2007.*(1./eZagg2_2007).^2;


%%%%%%%%%%%%%%%%2008 period 2%%%%%%%%%%%%%%%%%%%%
county2_2008 = shaperead('period2_2008');
names2 = {county2_2008.STATE};

%  convert school district locs to lat/lon
mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);%note: if you want to convert lat,lon to x,y, don't do this

indexHawaii2 = strcmp('02',names2);
indexAlaska2 = strcmp('15',names2);
indexPR2 = strcmp('72',names2);
indexConus2 = 1:numel(county2_2008);
indexConus2(indexHawaii2|indexAlaska2|indexPR2) = []; 
county2_2008 = county2_2008(indexConus2);


mstruct = defaultm('eqaconicstd');
mstruct.geoid = almanac('earth','grs80','foot');
mstruct.origin = [24.395833 -96.000000];
mstruct.falseeasting = 0.0;
mstruct.falsenorthing = 0.0;
mstruct.mapparallels = [29.5 45.5];
mstruct = defaultm(mstruct);
 for k=1:size(county2_2008,1)
    [lat lon] = minvtran(mstruct,county2_2008(k).X,county2_2008(k).Y);
    county2_2008(k).Lat = lat;
    county2_2008(k).Lon = lon;
 end
 

%data
countyd2=dlmread('period2_2008.csv');
% Zagg = sqrt(county(:,1));
eZagg2_2008 = countyd2(indexConus2,1);
%known survey variances
esigmavar2_2008 = (countyd2(indexConus2,2)/1.645).^2;
Zagg2_2008 = log(eZagg2_2008);
sigmavar2_2008 = esigmavar2_2008.*(1./eZagg2_2008).^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 G = county1_2006;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H1_2006 = H1;
 for k=1:ng
     tsum = sum(H1_2006(:,k));
     if tsum ~= 0
        H1_2006(:,k) = H1_2006(:,k)/tsum;
     end
 end

 G = county1_2007;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H1_2007 = H1;
 for k=1:ng
     tsum = sum(H1_2007(:,k));
     if tsum ~= 0
        H1_2007(:,k) = H1_2007(:,k)/tsum;
     end
 end
 
 G = county1_2008;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H1_2008 = H1;
 for k=1:ng
     tsum = sum(H1_2008(:,k));
     if tsum ~= 0
        H1_2008(:,k) = H1_2008(:,k)/tsum;
     end
 end
 
 
 G = county2_2007;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2007 = H1;
 for k=1:ng
     tsum = sum(H2_2007(:,k));
     if tsum ~= 0
        H2_2007(:,k) = H2_2007(:,k)/tsum;
     end
 end
 
 G = county2_2008;
 D = county3;
 %use MC area program: denomarea = 3 (G); use D x,y; use G lat lon
 H1 = Hmatrix_overlapSTRUCT(D,G,2,1,1,100);
 [nd,ng]=size(H1);
 H2_2008 = H1;
 for k=1:ng
     tsum = sum(H2_2008(:,k));
     if tsum ~= 0
        H2_2008(:,k) = H2_2008(:,k)/tsum;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save -v7.3 prep_sptm_alldataH;



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xi H
sigmavar = cat(1,esigmavar1,esigmavar1_2012,esigmavar1_2011,esigmavar1_2010,esigmavar1_2009,esigmavar1_2008,esigmavar1_2007,esigmavar1_2006,...
    esigmavar2_2012,esigmavar2_2011,esigmavar2_2010,esigmavar2_2009,esigmavar2_2008,esigmavar2_2007,...
    esigmavar3,esigmavar3_2012,esigmavar3_2011,esigmavar3_2010,esigmavar3_2009);
% 
% sigmavar = cat(1,sigmavar1,sigmavar1_2012,sigmavar1_2011,sigmavar1_2010,sigmavar1_2009,sigmavar1_2008,sigmavar1_2007,sigmavar1_2006,...
%     sigmavar2_2012,sigmavar2_2011,sigmavar2_2010,sigmavar2_2009,sigmavar2_2008,sigmavar2_2007,...
%     sigmavar3,sigmavar3_2012,sigmavar3_2011,sigmavar3_2010,sigmavar3_2009);

H = cat(1,H1',H1_2012',H1_2011',H1_2010',H1_2009',H1_2008',H1_2007',H1_2006',...
    H2_2012',H2_2011',H2_2010',H2_2009',H2_2008',H2_2007',...
    eye(3109),H3_2012',H3_2011',H3_2010',H3_2009');

Vinv = diag(1./sigmavar);

HpVinv = H'*Vinv; 
HpinvVH = HpVinv*H;

clear Vinv;

[EigHinvVHp,LamHpinvVH] = eig(HpinvVH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get MCMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%knots at each time point
level2=dlmread('knots250_ACS.csv');
t1 = [level2, 2012.5*ones(size(level2,1),1)];
t3 = [level2, 2011*ones(size(level2,1),1)]; 
t1_2012 = [level2, 2011.5*ones(size(level2,1),1)];
t2_2012 = [level2, 2011*ones(size(level2,1),1)]; 
t3_2012 = [level2, 2010*ones(size(level2,1),1)]; 
t1_2011 = [level2, 2010.5*ones(size(level2,1),1)];
t2_2011 = [level2, 2010*ones(size(level2,1),1)]; 
t3_2011 = [level2, 2009*ones(size(level2,1),1)]; 
t1_2010 = [level2, 2009.5*ones(size(level2,1),1)];
t2_2010 = [level2, 2009*ones(size(level2,1),1)]; 
t3_2010 = [level2, 2008*ones(size(level2,1),1)]; 
t1_2009 = [level2, 2008.5*ones(size(level2,1),1)];
t2_2009 = [level2, 2008*ones(size(level2,1),1)]; 
t3_2009 = [level2, 2007*ones(size(level2,1),1)]; 
t1_2008 = [level2, 2007.5*ones(size(level2,1),1)]; 
t2_2008 = [level2, 2007*ones(size(level2,1),1)]; 
t1_2007 = [level2, 2006.5*ones(size(level2,1),1)]; 
t2_2007 = [level2, 2006*ones(size(level2,1),1)]; 
t1_2006 = [level2, 2005.5*ones(size(level2,1),1)]; 

level = cat(1,t1,t3,t1_2012,t2_2012,t3_2012,t1_2011,t2_2011,t3_2011,...
    t1_2010,t2_2010,t3_2010,t1_2009,t2_2009,t3_2009,t1_2008,...
    t2_2008,t1_2007,t2_2007,t1_2006);


%intercept-only model
n = size(Zagg1,1)+size(Zagg3,1)+size(Zagg1_2012,1)+size(Zagg2_2012,1)+size(Zagg3_2012,1)+...
    size(Zagg1_2008,1)+size(Zagg1_2007,1)+size(Zagg1_2006,1)+...
    size(Zagg2_2008,1)+size(Zagg2_2007,1)+...
    size(Zagg1_2011,1)+size(Zagg2_2011,1)+size(Zagg3_2011,1)+...
    size(Zagg1_2010,1)+size(Zagg2_2010,1)+size(Zagg3_2010,1)+...
    size(Zagg1_2009,1)+size(Zagg2_2009,1)+size(Zagg3_2009,1);

X0 = ones(n,1);
X1 = ones(size(Zagg1,1)+size(Zagg1_2012,1)+size(Zagg1_2011,1)+size(Zagg1_2010,1)+size(Zagg1_2009,1)+...
    size(Zagg1_2008,1)+size(Zagg1_2007,1),1);
X2 = ones(size(Zagg2_2012,1)+size(Zagg2_2011,1)+size(Zagg2_2010,1)+size(Zagg2_2009,1)+size(Zagg2_2008,1)+size(Zagg2_2007,1),1);
X3 = ones(size(Zagg3,1)+size(Zagg3_2012,1)+size(Zagg3_2011,1)+size(Zagg3_2010,1)+size(Zagg3_2009,1),1); 
X = cat(2,X0,cat(1,X1,zeros(n-size(X1,1),1)),cat(1,zeros(size(X1,1),1),X2,zeros(n-size(X1,1)-size(X2,1),1)));
% 
% for i = 1:3109
% [ geom, iner, cpmo ] = polygeom( county3(i).Lat(isnan(county3(i).Lat)==0), county3(i).Lon(isnan(county3(i).Lon)==0));
% XA(i,1) = geom(1);
% XLA(i,1) = geom(2);
% XLO(i,1) = geom(3);
% end
% 
% X = cat(2,X,XA,XLA,XLO);

%Basis Functions

SGBF1= ArealBi2_spacetime(county1,2013,level,[],[],100,0.5,0.5);
SGBF3= ArealBi2_spacetime(county3,2009:2013,level,[],[],100,0.5,0.5);
SGBF1_2012= ArealBi2_spacetime(county1_2012,2012,level,[],[],100,0.5,0.5);
SGBF2_2012= ArealBi2_spacetime(county2_2012,2010:2012,level,[],[],100,0.5,0.5);
SGBF3_2012= ArealBi2_spacetime(county3_2012,2008:2012,level,[],[],100,0.5,0.5);
SGBF1_2011= ArealBi2_spacetime(county1_2011,2011,level,[],[],100,0.5,0.5);
SGBF2_2011= ArealBi2_spacetime(county2_2011,2009:2011,level,[],[],100,0.5,0.5);
SGBF3_2011= ArealBi2_spacetime(county3_2011,2007:2011,level,[],[],100,0.5,0.5);
SGBF1_2010= ArealBi2_spacetime(county1_2010,2010,level,[],[],100,0.5,0.5);
SGBF2_2010= ArealBi2_spacetime(county2_2010,2008:2010,level,[],[],100,0.5,0.5);
SGBF3_2010= ArealBi2_spacetime(county3_2010,2006:2010,level,[],[],100,0.5,0.5);
SGBF1_2009= ArealBi2_spacetime(county1_2009,2009,level,[],[],100,0.5,0.5);
SGBF2_2009= ArealBi2_spacetime(county2_2009,2007:2009,level,[],[],100,0.5,0.5);
SGBF3_2009= ArealBi2_spacetime(county3_2009,2005:2009,level,[],[],100,0.5,0.5);
SGBF1_2008= ArealBi2_spacetime(county1_2008,2008,level,[],[],100,0.5,0.5);
SGBF2_2008= ArealBi2_spacetime(county2_2008,2006:2008,level,[],[],100,0.5,0.5);
SGBF1_2007= ArealBi2_spacetime(county1_2007,2007,level,[],[],100,0.5,0.5);
SGBF2_2007= ArealBi2_spacetime(county2_2007,2005:2007,level,[],[],100,0.5,0.5);
SGBF1_2006= ArealBi2_spacetime(county1_2006,2006,level,[],[],100,0.5,0.5);

% S = cat(1,SGBF1,SGBF3,SGBF1_2012,SGBF2_2012,SGBF3_2012,...
%     SGBF1_2011,SGBF2_2011,SGBF3_2011,SGBF1_2010,SGBF2_2010,SGBF3_2010,...
%     SGBF1_2009,SGBF2_2009,SGBF3_2009);

S = cat(1,SGBF1,SGBF1_2012,SGBF1_2011,SGBF1_2010,SGBF1_2009,SGBF1_2008,SGBF1_2007,SGBF1_2006,...
    SGBF2_2012,SGBF2_2011,SGBF2_2010,SGBF2_2009,SGBF2_2008,SGBF2_2007,...
    SGBF3,SGBF3_2012,SGBF3_2011,SGBF3_2010,SGBF3_2009);

[S1,idx]=licols(S);

Sconnector1 = ArealBi2_spacetime(county3,2005,level,[],[],100,0.5,0.5); 
Sconnector2 = ArealBi2_spacetime(county3,2006,level,[],[],100,0.5,0.5); 
Sconnector3 = ArealBi2_spacetime(county3,2007,level,[],[],100,0.5,0.5); 
Sconnector4 = ArealBi2_spacetime(county3,2008,level,[],[],100,0.5,0.5); 
Sconnector5 = ArealBi2_spacetime(county3,2009,level,[],[],100,0.5,0.5); 
Sconnector6 = ArealBi2_spacetime(county3,2010,level,[],[],100,0.5,0.5); 
Sconnector7 = ArealBi2_spacetime(county3,2011,level,[],[],100,0.5,0.5); 
Sconnector8 = ArealBi2_spacetime(county3,2012,level,[],[],100,0.5,0.5); 
Sconnector9 = ArealBi2_spacetime(county3,2013,level,[],[],100,0.5,0.5); 

Sconnector = cat(1,Sconnector1,Sconnector2,Sconnector3,Sconnector4,Sconnector5,...
    Sconnector6,Sconnector7,Sconnector8,Sconnector9);
Sconnectorf = Sconnector(:,idx);

% Full Model
load('countAdj.txt')
for j = 1:3109
if sum(countAdj(j,:))>0
countAdj(j,:) = countAdj(j,:)/sum(countAdj(j,:));
end
end
Q = eye(3109) - 0.9*countAdj;
Qinv = pinv(Q);

[PQ,LQ] = eig(Q);

% FullCovar = blkdiag(Qinv,Qinv,Qinv,Qinv,Qinv);
% K = pinv(Sconnector'*Sconnector)*Sconnector'*FullCovar*Sconnector*pinv(Sconnector'*Sconnector);

%Moran's I Propagator
B = [PQ'*ones(3109,1),eye(3109)];
BpB = B'*B;
BpBinv = pinv(BpB);
GBI = B*BpBinv*B';
[M,LM] = eig(GBI);
M = real(M);

%Target Covariance
Kinv=make_full_model_sptcovar_9(Q,M,Sconnectorf,3109);
[P,D] = eig(Kinv);
D = real(diag(D));
D(D<0) = 0;
Dinv = D;
Dinv(D>0) = (1./D(D>0));
K = real(P)*diag(Dinv)*real(P');
Kinv = real(P)*diag(D)*real(P');
save Ktarget K;

%get givens angles of K
% De=logit_linear_design(size(S,2));

%Full data and variance vectors
% Zagg = cat(1,eZagg1,eZagg3,eZagg1_2012,eZagg2_2012,eZagg3_2012,...
%     eZagg1_2011,eZagg2_2011,eZagg3_2011,eZagg1_2010,eZagg2_2010,eZagg3_2010,...
%     eZagg1_2009,eZagg2_2009,eZagg3_2009);
% sigmavar = cat(1,esigmavar1,esigmavar3,esigmavar1_2012,esigmavar2_2012,esigmavar3_2012,...
%     esigmavar1_2011,esigmavar2_2011,esigmavar3_2011,esigmavar1_2010,esigmavar2_2010,esigmavar3_2010,...
%     esigmavar1_2009,esigmavar2_2009,esigmavar3_2009);

Zagg = cat(1,eZagg1,eZagg1_2012,eZagg1_2011,eZagg1_2010,eZagg1_2009,eZagg1_2008,eZagg1_2007,eZagg1_2006,...
    eZagg2_2012,eZagg2_2011,eZagg2_2010,eZagg2_2009,eZagg2_2008,eZagg2_2007,...
    eZagg3,eZagg3_2012,eZagg3_2011,eZagg3_2010,eZagg3_2009);
% 
% Zagg = cat(1,Zagg1,Zagg1_2012,Zagg1_2011,Zagg1_2010,Zagg1_2009,Zagg1_2008,Zagg1_2007,Zagg1_2006,...
%     Zagg2_2012,Zagg2_2011,Zagg2_2010,Zagg2_2009,Zagg2_2008,Zagg2_2007,...
%     Zagg3,Zagg3_2012,Zagg3_2011,Zagg3_2010,Zagg3_2009);

save -v7.3 prep_sptm_alldataH_7815;