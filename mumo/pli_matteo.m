function m=pli_matteo(a)
% a is a filtered multichannel signal (time x channels)
% hilbert(a) calculates analytic signal (complex valued) of each 
% column of a. Phase Lag Index between channel i and j averaged over 
% time bins is stored in m(i,j) 
% number of channels 
N=size(a,2);
nch=size(a,2);
m(1:N,1:N)=0; 
complex_a=hilbert(a); 
for i=1:nch
    for j=1:nch
        if i<j
            
            %idx=aux<0.01
            %aux(idx)=0;
            %aux=floor(aux);
            m(i,j)=abs( mean (sign( (imag(complex_a(:,i)./complex_a(:,j)) ) ) ) ) ;
            
        end 
    end 
end 
m=m+m';
