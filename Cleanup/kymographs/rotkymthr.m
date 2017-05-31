%% ROTary KYMograph THResholder

I = imread('unprocessedkym.tif');
T = zeros(size(I));
F = zeros(size(I));

m = mean(I,2);
med  = median(I,2);
d = sqrt(var(I,0,2));

for rr = 1:size(I,1)
    T(rr,:) = I(rr,:) > (m(rr)+2*d(rr));
    
    % filter subtract mean of each row from the row
    F(rr,:) = I(rr,:) - m(rr);
end

% go through each gaussian blurred column of filtered kymograph and return max pixel index
thetas = zeros(1,size(I,2));
for cc = 1:size(F,2)
    [x,ind] = max(smooth(F(:,cc),6));
    thetas(cc) = ind*360/72;
end
thetas = thetas*pi/180;

%%
% subplot(311)
% plot(cos(unwrap(thetas(1:100))));
% subplot(312)
% plot(sin(unwrap(thetas(1:100))));
% subplot(313)
% plot(unwrap(thetas(1:100)));

%% plots
subplot(311)
u =smooth(unwrap(thetas(1:300)),10)*180/pi;
plot(u);
xlabel('frame')
title('unwrapped angle (degrees)')
subplot(312)
plot(diff(u));
title('angular velocity (degrees/sec)')
xlabel('frame');

T = diff(u)>0;
subplot(313)
% 100 ms median filter, since tumbles are 100-300ms
plot(medfilt1(double(T),8))
title('direction (100 msec median filtered)')
xlabel('frame')

figure(2)
subplot(211)
imagesc(I(:,1:300))
title('raw kymograph')
colormap gray
subplot(212)
imagesc(F(:,1:300))
title('filtered kymograph')