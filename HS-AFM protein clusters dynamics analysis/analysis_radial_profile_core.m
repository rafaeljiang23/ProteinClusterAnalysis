%% radial profile analysis part a: core analysis
%%% This script takes a pre-processed mask file, BW2, from the main script
%%% for the core patch dynamics analysis. The outputs include the patch
%%% line tension and leading velocity measurements.


%%%% add empty pixels to the x- and y- dimension to the mask file
[d1, d2, d3] = size(BW2);

BW2a = zeros(d1+20, d2+20, d3);
BW2a(11:end-10, 11:end-10, :) = BW2;
BW2 = BW2a;

%% convert to radial profile a (probability)
%%%% the (probability) radial profile a should be a R(radias)-by-A(angle)-by-T(time) matrix
%%%% of doubles indicating the probability of the pixels being at the patch
%%%% circumference

[d1, d2, d3] = size(BW2);

%%%% radial profile: patch radius should be the pixels having 0 distance to the patch mask
D = BW2.*0;
for i = 1:d3
    D(:, :, i) = bwdist(double(BW2(:, :, i) > 0));
end
circum = double(D == 1);    % radial profile a
circum = imgaussfilt(circum);    % apply a gaussian filter to radial profile a

% MIJ.createImage(circum);   % MIJI, display radial profile a
circum_BW = double(circum > 0);

%%%% resize the radial profile a
d4 = 1000;   % d4: new angular dimension (number of pixels corresponding t0 360 degrees)
d5 = 100;   % new radial dimension (number of pixezls corresponding the the radius)
circum2 = zeros(d5, d4, d3);   % resized radial profile data

[xx, yy] = ndgrid(1:d1, 1:d2);
ctr = zeros(d3, 2);    % normalized center of mass of the patch
for i = 1:d3
    frame2 = bwmorph(squeeze(circum_BW(:, :, i)), "clean");
    sel = squeeze(BW2(:, :, i)) > 0;
    xc = mean(xx(sel));    % center of mass (x) of the patch
    yc = mean(yy(sel));    % center of mass (y) of the patch
    ctr(i, 1) = xc-d1/2;
    ctr(i, 2) = yc-d2/2;
    frame3 = im_cart2pol(frame2, [xc, yc]);   % cartesian to polar conversion of the patch pixels
    sz = size(frame3);   % find the dimension sizes of the patch radial profile
    frame3 = imresize(frame3, [sz(1), d4]);   % resize patch radial profile
    circum2(1:sz(1), :, i) = frame3;
end
circum2 = imgaussfilt(circum2);   % apply a gaussian filter to the resized radial profile a
% MIJ.createImage(circum2);   % MIJI, display resized radial profile data a


%%  create radial kymograph and radial profile b (absolute)
%%%% the kymograph should be a A(angle)-by-T(time) matrix of doubles
%%%% indicating R(A) - radius as a funcion of angle - at each time
%%%% point

%%%% the (absolute) radial profile b should be a
%%%% R(radias)-by-A(angle)-by-T(time) matrix of values 0 and 1, indicting
%%%% the absolute radius value at each angle point at each time point

kymo = zeros(d4, d3);   % kymograph data
circum3 = circum2.*0;   %  radial profile b (absolute)
for i = 1:d4
    for t = 1:d3
        line = squeeze(circum2(:, i, t));
        l = find(line > 0);
        %%% probability weighted radius at time t and angle i
        kymo(i, t) = sum(line(l).*l./sum(line(l)));
        r = sum(line(l).*l./sum(line(l)));
        circum3(round(r), i, t) = 1;
    end
end

%MIJ.createImage(kymo');   % MIJI, display kymograph of the radial profile
%MIJ.createImage(circum3);   % MIJI, display radial profile b

%% Fourier analysis
%%%% Fourier analysis of the kymograph to characterize the patch contour
%%%% fluctuation

%%%% Fourier fitting setup: Fourier analysis with 16 frequency
fittypef16 = fittype(['a0 + a1*cos(x) + b1*sin(x) +a2*cos(2*x) + b2*sin(2*x) + ' ...
    'a3*cos(3*x) + b3*sin(3*x) + a4*cos(4*x) + b4*sin(4*x) + ' ...
    'a5*cos(5*x) + b5*sin(5*x) + a6*cos(6*x) + b6*sin(6*x) + ' ...
    'a7*cos(7*x) + b7*sin(7*x) + a8*cos(8*x) + b8*sin(8*x) + ' ...
    'a9*cos(9*x) + b9*sin(9*x) + a10*cos(10*x) + b10*sin(10*x) + ' ...
    'a11*cos(11*x) + b11*sin(11*x) + a12*cos(12*x) + b12*sin(12*x) + '...
    'a13*cos(13*x) + b13*sin(13*x) + a14*cos(14*x) + b14*sin(14*x) + '...
    'a15*cos(15*x) + b15*sin(15*x) + a16*cos(16*x) + b16*sin(16*x)'],...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a0' 'a1' 'b1' 'a2' 'b2' 'a3' 'b3' 'a4' 'b4'...
    'a5' 'b5' 'a6' 'b6' 'a7' 'b7' 'a8' 'b8'...
    'a9' 'b9' 'a10' 'b10' 'a11' 'b11' 'a12' 'b12' ...
    'a13' 'b13' 'a14' 'b14' 'a15' 'b15' 'a16' 'b16'});


optionsf = fitoptions('method','NonlinearLeastSquares', 'normalize', 'off');
Nfreq = 16;


%%%% Fourier fitting of radius profile (absolute)
aa = zeros(Nfreq, d3);
bb = zeros(Nfreq, d3);
radius = res*sqrt(sum(BW2(:))/(pi*d3));

figure(1)
hold on
t0 = 1;
iii = 1;
ttt = randi(t-1, 9, 1) + 1;
A = linspace(-pi, pi, d4)';

warning('off')
for t = 2:d3
    disp(100*t/d3 + "%")

    Xt1 = kymo(:, t);    % kymograph at time t
    Xt0 = kymo(:, t-1);      % kymograph at time t-1
    Xt1 = movmean(Xt1, 10);    % walking averaged kymograph at time t
    Xt0 = movmean(Xt0, 10);    % walking averaged kymograph at time t-1

    X = (Xt1 -Xt0)./Xt0;   % relative delta radius

    if ismember(t, ttt)
        iii = iii + 1;
        plot(X, "Color", [iii/10 iii/10 iii/10]);   % plot relative delta radius
    end
    %%% normalize relative delta radius (for fitting)
    Xmean = mean(X);
    Xstd = std(X);
    X = (X - Xmean)./Xstd;

    %%% Fourier fitting
    fourierseries = fit(A,X,fittypef16, optionsf);
    coeffvals = coeffvalues(fourierseries);
    coeffvals = coeffvals*Xstd;


    %%% prefactors for sine and cosine part
    for i = 1:Nfreq
        aa(i, t) = coeffvals(2*(i-1) + 2);
        bb(i, t)= coeffvals(2*(i-1) + 3);
    end
end
warning('on')
hold off
%% tension calculation
%%%% For the Fourier fitting, calculate the line tension for the patch
%%%% see tension determination equation in the Methods section

nn = 1:Nfreq;
xx = 1./(nn'.^2 - 1);   % x values
yy = xx.*0;   % y values
ee = xx.*0;   % error values
for i = 1:Nfreq
    aai = aa(i, :);
    bbi = bb(i, :);
    yyi = aai.^2+bbi.^2;
    radiusi = res*sqrt(squeeze(sum(BW2, [1 2]))/pi);
    yyi = yyi .* radiusi';
    yy(i) = mean(yyi);
    ee(i) = std(yyi)/sqrt(numel(yyi));
end

%%%% line fitting setup
%%% select frequency 3-8, Note that when nn = 1, xx = Inf and thus must be
%%% ignored. In practice, ignoring frequency nn = 2 gives better linear
%%% fitting.

xx1 = xx(3:8); 
yy1 = yy(3:8);
ee1 = ee(3:8);

foption = fitoptions('method','NonlinearLeastSquares', 'normalize', 'off',...
    'Weight', 1./ee1,'Exclude', []);
fittype1 = fittype('a*x + b',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a' 'b'});

%%%% normalize x and y values (for fitting)
xmean = mean(xx1);
ymean = mean(yy1);
xstd = std(xx1);
ystd = std(yy1);
xx2 = (xx1 - xmean)/xstd;
yy2 = (yy1 - ymean)/ystd;

%%%% line fitting
linefit = fit(xx2, yy2,fittype1,foption);
coeffline = coeffvalues(linefit);

%%%% plotting fitting result
figure(2)
hold on
errorbar(xx1, radius*yy1, -radius*ee1, radius*ee1, "k.",  "LineWidth", 1, "MarkerSize", 10);
xx3 = linspace(xx2(1), xx2(end), numel(xx1)*100);
yy3 = xx3 .* coeffline(1) + coeffline(2);
xx3 = xx3*xstd + xmean;
yy3 = yy3*ystd + ymean;
plot(xx3, radius*yy3, "k", "LineWidth", 1);
xlabel("1/(n^2-1)");
ylabel("(<a_n^2>+<b_n^2>)r0 (nm)");
set(gca, "LineWidth", 1)
set(gca, "FontSize", 12)


%%%% determine tension values
slope = coeffline(1)*ystd/xstd;
lambda = 2/(slope*pi);    % tension, unit: kBT/nm

disp("Line tension measurement: " + lambda + " kBT/nm");
%% Edge leading velocity
%%% Calculate the edge leading velocity. Note: User should consider the
%%% HS-AFM imaging rate (unit: frame/s) when reporting the leading veloctiy
%%% measurement in physical unit (nm/s).

Vlead = abs(kymo(:, 2:end) - kymo(:, 1:end-1));    % leading velocity, unit: pixel/frame
Vlead = mean(Vlead, "all") * res;    % leading velocity, unit: nm/frame

disp("Edge leading velocity measurement: " + Vlead + " nm/frame");
