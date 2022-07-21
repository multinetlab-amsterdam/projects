function ENV_corr = envelopecorr(data,sm_param)
% compute the amplitude envelope correlation 
% pairwise leakage correction
% sm_param is smoothing parameter in samples
ENV_corr = zeros(size(data,2),size(data,2));
for i = 1: size(data,2)
    
    % data and obtain envelope for region i
    x = data(:,i);
    htx = hilbert(x); 
    htx = smooth(htx,sm_param);
    envelope_x = sqrt(real(htx).^2 + imag(htx).^2);
    
    for j = 1: size(data,2)
    if i~=j
     
     y = data(:,j);
     [b,bint,r] = regress(y,x);

    % step 1, compute the hilbert transform and bring him to the origin
    hty = hilbert(r);
    hty = smooth(hty,sm_param);
        
    % step 2, compute envelope data
    envelope_y = sqrt(real(hty).^2 + imag(hty).^2);
    
    % correlation between envelopes
    ENV_corr(i,j) = corr(envelope_x,envelope_y);
            
    end
%     disp(j)
    end
%    fprintf('done for outer loop %d \n',i)

end
ENV_corr = (ENV_corr + ENV_corr')/2;
end