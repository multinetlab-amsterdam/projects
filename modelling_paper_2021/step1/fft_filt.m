function [f, rel_power] = fft_filt(timeseries, Fs,thres1, thres2, nr_rois)

% This function filters timeseries between the given frequency
% bands (thres1 and thres2)

roi = 1:nr_rois;
f = zeros(numel(timeseries(:,1)), numel(roi));
abs_power = zeros(numel(roi),1); 
power_05_48 = zeros(numel(roi),1);
rel_power = zeros(numel(roi),1); 

freq = 1/length(timeseries):Fs/length(timeseries):Fs;
% include only freq. band of interest
[~, inx_thr1] = min(abs(freq-thres1));
[~, inx_thr2] = min(abs(freq-thres2)); 
        
for roi = 1:nr_rois
    xdft = fft(timeseries(:,roi));
    abs_power(roi) = sum(abs(xdft(inx_thr1:inx_thr2)));
    power_05_48(roi) = sum(abs(xdft(8:630)));
    rel_power(roi) = abs_power(roi)/power_05_48(roi);
                
    % remove all other frequencies
    xdft(1:inx_thr1) = 0;
    xdft(inx_thr2:length(xdft)-inx_thr2) = 0;
    xdft(length(xdft)-inx_thr1:length(xdft)) = 0;                
              
    f(:,roi)= real(ifft2(xdft, 'symmetric'));
end

end 