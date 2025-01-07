%% patch tilt and nearest neighbor analysis

%%% This script take a scaled and averaged image of the patch, data, as the
%%% main input file for analysis. User should pick particles from the 'data'
%%% file and determine the particle x-, y-coordinates as well as the
%%% molecular symmetry value (saved as 'X', 'Y', and 'nf' array variables)

%%% MIJI is recommanded for visualizing the analysis (MIJ) but not require.

data;  % Input scaled HS-AFM image of the patch (averaged frame)
Y;  % Input particle coordinate X-axis, a N-by-1 array for N particles
X;  % Input particle coordinate Y-axis, a N-by-1 array for N particles
nf;  % Input particle molecular symmetry array, a N-by-1 array for N particles

%%%% define useful parameters
% enable analysis of oligomer mixture
scale = 5;    % input data scale factor
radius = 35;  % circular mask outer radius
radius2 = 5;  % circular mask inner radius
allowance = 2;   % particle center adjustment allowance
res = 0.5;   % spatial resolution, unit: nm/pixel

close all
% MIJ.closeAllWindows();   % MIJI


%% determine the protomer locations
[xx, yy] = ndgrid(-radius:radius, -radius:radius);
circle_sel = xx.^2 + yy.^2 < radius.^2;    % outer circular mask
circle2_sel = xx.^2 + yy.^2 > radius2.^2;     % inner circular mask

XX = X;
YY = Y;

data = data-min(data(:));   % all values should be non-negative
data2 = data.*0;   % record particle center
data3 = data.*0;   % record protomer center

particles = [];   % collect symmetriezed particles  (not very useful)
particles2 = [];   % collect symmetriezed and masked particles

AA = zeros(numel(XX), 1);   % record protomer angle
RR = zeros(numel(XX), 1);   % record protomer radius

for i = 1:numel(XX)
    xi = round(XX(i));
    yi = round(YY(i));
    nfi = nf(i);

    frame = data(yi-radius: yi+radius, xi-radius: xi+radius);    % crop particle i
    frame2 =  nfold(frame, nfi);   % apply n-fold symmetry to particle i
    frame2 = frame2.*circle_sel;   % apply circular mask to the symmetrized particle i

    particles{i} = frame;
    particles2{i} = frame2;

    %%%% partile center adjustment (optional): 
    % must have the highest internal symmtery score 
    % (cross-correlation to the symmetrized particle)

    score = 0;  
    yjm = 0;
    xjm = 0;
    for yj = -allowance:allowance
        for xj = -allowance:allowance
            frame2 = data(yi-radius+yj: yi+radius+yj, xi-radius+xj: xi+radius+xj);   % crop particle i 
            frame3 = nfold(frame2, nfi);   % apply n-fold symmetry to particle i
            frame2 = frame2.*circle_sel;   % apply circular mask to particle i
            frame3 = frame3.*circle_sel;   % apply circular mask to the symmetrized particle i
            scorej = corr2(frame2,frame3);   % calculate cross-correlation (internal symmetry value)
            %%%%%% update new center position adjustment value if reaching a higher score
            if scorej > score
                score = scorej;
                particles2{i} = frame3;     % update symmetriezed and masked particles
                yjm = yj;
                xjm = xj;
            end
        end
    end

    XX(i) = xi + xjm;   % adjust center x coordinate
    YY(i) = yi + yjm;   % adjust center y coordinate

    %%%% protomer center search
    frame2 = particles2{i};
    frame2 = frame2 .* circle2_sel;   % apply (inner) circular mask to the circular and symmetrized particle i

    %%%%%%  denoise and find local maxima
    frame3 = double(frame2 > prctile(frame2(:), 95));   % find local maxima
    frame3 = double(imgaussfilt(frame3) > 0.0001);
    D = bwdist(~frame3);
    [ym, xm] = find(D == max(D, [], "all"), 1);   % obtain (one) local maximum coordinates (cartesian)
    [angle, rho] = cart2pol(xm - radius - 1, ym - radius - 1);   % convert to polar coordinates
    angle = mod(angle + 2*pi, 2*pi/nfi);   % normalize angle value to 0 - (360/nf) degrees
    if angle >= pi/nf
        angle = angle - 2*pi/nfi;
    end

    AA(i) = angle;    % protomer angle
    RR(i) = rho;    % protomer radius
    %%%% display particle center and protomer center 
    data2(YY(i), XX(i)) = 1;
    for j = 1:nfi
        data3(YY(i) + round(radius * sin(angle + j*(2*pi/nfi))/2), XX(i) + round(radius * cos(angle + j*(2*pi/nfi))/2)) = 1;
    end
end

data2 = 50*double(imgaussfilt(data2, 2));   % particle center
data3 = 10*double(imgaussfilt(data3, 1));   % protomer center

%  MIJ.createImage(data + data2 + data3);  % MIJI, display the particle/protomer assignement results


%%%% update 'XYA' file recording particle info
XYA = [XX/scale YY/scale AA RR/scale];   


%% determine particle pair based on delaunay triangulation
%%%% obtain particle info from 'XYA' file
XX = XYA(:, 1);
YY = XYA(:, 2);
AA = XYA(:, 3);
RR = XYA(:, 4);

%%%% delaunay triangulation 
DT = delaunay(XX, YY);

%%%% plot delaunay triangulation results
label = [1:numel(XX)]';
label = num2str(label);
label = cellstr(label);
figure(20)
hold on
scatter(XX, YY);
text(XX+1, YY+1, label);
triplot(DT,XX,YY);
hold off
axis ij

%%%% find all paired particles
pairs = zeros(numel(DT), 2);
it = 0;
for i = 1:numel(DT)/3
    DT(i, :) = sort(DT(i, :));
    it = it + 1;
    pairs(it, 1) = DT(i, 1);
    pairs(it, 2) = DT(i, 2);
    it = it + 1;
    pairs(it, 1) = DT(i, 1);
    pairs(it, 2) = DT(i, 3);
    it = it + 1;
    pairs(it, 1) = DT(i, 2);
    pairs(it, 2) = DT(i, 3);
end

pairs_dt = unique(pairs,'rows');

%% Manual! Remove wrong pairs
%%%% Note that wrong pairs should be MANNUALLY REMOVED at this step
%%% User should work on 'pairs_update'
pairs_update = pairs_dt;

%%% Here we provide a ground true 'pairs' variable for the TEST run only!
pairs = pairs_update;
pairs = test_data_pairs;   % User should comment out this line (for TEST run only)
%% particle pair analysis

sz = size(pairs);
P = zeros(sz(1), 7);   % output file recording pair info

data4 = data.*0;   % record particle pairs
for i = 1:sz(1)
    %%%% read particle pair coordinates
    x1 = XX(pairs(i, 1));    % particle 1 x coordinate
    x2 = XX(pairs(i, 2));    % particle 2 x coordinate
    y1 = YY(pairs(i, 1));    % particle 1 y coordinate
    y2 = YY(pairs(i, 2));    % particle 2 y coordinate
    a1 = AA(pairs(i, 1));    % particle 1 angle
    a2 = AA(pairs(i, 2));    % particle 2 angle
    nf1 = nf(pairs(i, 1));    % particle 1 molecular symmetry
    nf2 = nf(pairs(i, 1));    % particle 1 molecular symmetry
    %%%% display particle pairs
    x3 = linspace(x1, x2, 20*scale);
    y3 = linspace(y1, y2, 20*scale);
    for k = 1:numel(x3)
        xx3 = round(x3(k)*scale);
        yy3 = round(y3(k)*scale);
        data4(yy3, xx3) = 1;
    end

    %%%% record particle pair info
    P(i, 1:2) = pairs(i, :);    % particle pair id
    %%%%%% relative spatial relationship between two paired particles
    P(i, 3) = sqrt((x1-x2)^2 + (y1-y2)^2);    % particle pair distance (ctr-to-ctr)
    P(i, 4) = atan2(y2-y1, x2-x1);    % particle pair angle (alpha), 2 relative to 1
    P(i, 5) = atan2(y1-y2, x1-x2);    % particle pair angle (alpha), 1 relative to 2, should be pi - 'particle pair angle, 2 relative to 1'%

    %%%%%% relative spatial relationship between the protomer in one paired
    %%%%%% particle to the position of the other particle 
    P(i, 6) = P(i, 4) - a1;    % particle 1, protomer angle (theta) relative to particle pair
    P(i, 7) = P(i, 5) - a2;    % particle 2, protomer angle (theta) relative to particle pair
    
    %%%%%% normalize angle to -180/nf to 180/nf degrees
    if P(i, 6) >= 0
        P(i, 6) = mod(P(i, 6),2*pi/nf1);
    elseif P(i, 6) < 0
        P(i, 6) = -mod(-P(i, 6),2*pi/nf1);
    end

    if P(i, 6) > pi/nf1
        P(i, 6) = P(i, 6) - 2*pi/nf1;
    elseif P(i, 6) < -pi/nf1
        P(i, 6) = P(i, 6) + 2*pi/nf1;
    end

    if P(i, 7) >= 0
        P(i, 7) = mod(P(i, 7),2*pi/nf2);
    elseif P(i, 7) < 0
        P(i, 7) = -mod(-P(i, 7),2*pi/nf2);
    end

    if P(i, 7) > pi/nf2
        P(i, 7) = P(i, 7) - 2*pi/nf2;
    elseif P(i, 7) < -pi/nf2
        P(i, 7) = P(i, 7) + 2*pi/nf2;
    end
end


%MIJ.createImage(data + data2 + data3 + data4);   % MIJI, display particle/protomer/pair results

disp("particle nearest neighbor analysis done")

%% measure particle tilt

sz = size(XYA);
omega = zeros(sz(1), 1);    % record particle tilt

for i = 1:sz(1)
    %%%%%% read particle info
    xi = XYA(i, 1) * scale;   % particle x coordinate
    yi = XYA(i, 2) * scale;   % particle y coordinate
    thetai = XYA(i, 3);    % particle protomer angle
    rhoi = XYA(i, 4) * scale;   %  particle protomer radius
    nfi = nf(i);   %  particle molecular symmetry
    %%%%%% collect protomer coordinates
    coors = zeros(nfi, 3);
    particle = data(yi-radius:yi+radius, xi-radius:xi+radius);   % crop particle
    for j = 1:nfi
        BW = particle.*0;

        xj = rhoi*cos(thetai + (j-1)*2*pi/nfi);   % protomer j, x coordinate, unit: pixel
        yj = rhoi*sin(thetai + (j-1)*2*pi/nfi);   % protomer j, y coordinate, unit: pixel
        coors(j, 1) = res*xj/scale;     % protomer j, x coordinate, unit: nm
        coors(j, 2) = res*yj/scale;     % protomer j, y coordinate, unit: nm

        xj = radius + round(xj);   % convert to integer
        yj = radius + round(yj);  % convert to integer
        BW(yj, xj) = 1;

        D = bwdist(BW);     % calculate pixel distance to the protomer location
        D = D < 3;         % all pixels within 3 pixel of the protomer location are included in height (z-value) determination
        coors(j, 3) = mean(particle(D));    % protomer height (z-value) is the average of all included pixels 
    end

    C = nchoosek(1:nfi,3);    % find all combinations of vector pairs (three coordinates required for vector pairs)
    omegaj = zeros(numel(C)/3, 1);
    for j = 1:numel(C)/3
        sel = C(j, :);
        coorsj = coors(sel, :);
        v1 = coorsj(3, :) - coorsj(2, :);   % define vector 1 in vector pair j
        v2 = coorsj(2, :) - coorsj(1, :);   % define vector 2 in vector pair j
        n0 = cross(v1, v2);   % cross product of the vectors (vector 3) gives another vector normal to the vector 1-2 plane
        omegaj(j) = acos(dot(n0, [0 0 1])./norm(n0));    % tilt of vector pair j is the projection of vector 3 to the z-axis
    end
    omega(i) = mean(omegaj);     % tilt is the averaged tilt values from all vector pairs
end


%%%% update 'XYA' file recording particle info and molecular symmetry
XYA(:, 5) = omega;
XYA(:, 6) = nf;

disp("particle tilt analysis done")


%% Output
%%%% The main outputs from this script is matrix 'XYA' that records all
%%%% single-partile related measurements and matrix 'P' that records all
%%%% nearest neighbor pair measurements.
%%%% 'XYA' has a dimension of N-by-6, for N particles (row). 
%%%% Values stored in each column are as follows:
%%%%   col1: particle coordinate x (unit: pixel, unscaled)
%%%%   col2: particle coordinate y (unit: pixel, unscaled)
%%%%   col3: protomer coordinate angle (unit: radian)
%%%%   col4: protomer coordinate radius (unit: pixel, unscaled)
%%%%   col5: particle tilt omega (unit: radian, not normalized)
%%%%   col6: particle molecular symmetry
%%%% 'P' has a dimension of M-by-7, for M particle (nearest neighbor) pairs (row). 
%%%% Values stored in each column are as follows:
%%%%   col1: particle pair id p1
%%%%   col2: particle pair id p2
%%%%   col3: particle pair distance (ctr-to-ctr)
%%%%   col4: particle pair angle, p2 relative to p1, unit: radian (not normalized)
%%%%   col5: particle pair angle, p1 relative to p2, unit: radian (not normalized)
%%%%   col6: particle 1, protomer angle (theta) relative to particle pair, unit: radian (normalized)
%%%%   col7: particle 2, protomer angle (theta) relative to particle pair, unit: radian (normalized)

%% helper functions

%%% n-fold symmetry 
function out = nfold(in, nf)
out = in;

for i = 1:nf-1
    out = out + imrotate(in, i*360/nf, "bicubic", "crop");

end

out = out/nf;
end