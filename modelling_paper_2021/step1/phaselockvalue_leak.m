function PLV = phaselockvalue_leak(data)
% compute PLV
PLV = zeros(size(data,2),size(data,2));

for m = 1: size(data,2)
    
    % data and obtain phase for region m
    x = data(:,m);
    htx=hilbert(x); htx=htx-mean(htx);
    phase1 = angle(htx);
    phase1(1:10)=[];phase1(end-10:end)=[];
    
    for j = 1: size(data,2)
        
    if m<j

    % data
    y = data(:,j);
    [b,bint,r] = regress(y,x);
    
    % step 1, compute the hilbert transform and bring him to the origin
    hty=hilbert(r); hty=hty-mean(hty);
    
    % step 2, compute the instantenous phase
    phase2 = angle(hty);
    phase2(1:10)=[];phase2(end-10:end)=[];   
    
    % compute plv
    RP = phase1 - phase2;                 % relative phase
    PLV(m,j) = abs(sum(exp(1i*RP))/length(RP));
    clear phase2 RP
    end
    
    end
end
    
PLV = PLV + PLV';

end