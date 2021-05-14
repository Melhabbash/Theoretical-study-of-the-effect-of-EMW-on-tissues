clear all
clc

    
permeability=4*pi*1e-7;  %Free Space
permitivity=8.85e-12;    %Free Space
c=1/sqrt(permeability*permitivity);           
                      
m=45; 
ex=zeros(1,m);
hy=zeros(1,m);
cb=zeros(1,m);
eaf=zeros(1,m);
s=zeros(1,m);
dt=zeros(1,m);
p=zeros(1,m);
SAR=zeros(1,m);

dx = .01;
dt = dx/(2*c);
d0 = 0.3;                       
                                  
%region

k1=fix(d0/dx);% blood



% free space
permitivity(1,1:k1)=1;
dt(1,1:k1) = dx/(2*c);  
s(1,1:k1) = 0;

% blood
permitivity(1,k1:m)=  59.372  ;
dt(1,k1:m) = dx/(2*c); 
s(1,k1:m) = 2.0435;






ex_left_m1 = 0;
ex_left_m2 = 0;
ex_right_m1 = 0;
ex_right_m2 = 0;
ex_left_b1_m1 = 0;
ex_left_b1_m2 = 0;


iteration = 500;                       

for n = 1:iteration
    for k = 2:45
        eaf =dt(1,k)*s(1,k)/(2*permitivity(1,k)*8.85e-12);
        ca=(1-eaf)/(1+eaf);
        cb=0.5/(permitivity(1,k)*(1+eaf));
        
        ex(1,k) = ca*ex(1,k) + cb*(hy(1,k-1) - hy(1,k));
        p(1,k)=(ex(1,k)^2)*s(1,k)/2;
        SAR(1,k)=((ex(1,k)^2)*s(1,k))/(2*(1060));
    end
    source = sin(2*pi*1800*10^(6)*dt(1,1)*n);       
     ex(1,2) = source;
    
    % Boundaries
    if n>2
        KE=45;
        kc=30;
        % Left boundary
        ex(1,1) = ex_left_m2;      % set left boundary of ex(1) by the value of ex(2)@n-2
        ex_left_m2 = ex_left_m1;   % set the value of ex(2)@n-2 by the value of ex(2)@n-1
        ex_left_m1 = ex(1,2);      % set the value of ex(2)@n-1 by the value of ex(2)@n (current value)
        % Right boundary
        ex(1,KE) = ex_right_m2;    % set right boundary of ex(201) by the value of ex(200)@n-2
        ex_right_m2 = ex_right_m1; % set the value of ex(200)@n-2 by the value of ex(200)@n-1
        ex_right_m1 = ex(1,KE-1);  % set the value of ex(200)@n-1 by the value of ex(200)@n (current value)
        
        % dielectric edge
        ex(1,kc) = ex_left_b1_m2;
        ex_left_b1_m2 = ex_left_b1_m1;
        ex_left_b1_m1 = ex(1,kc-1);
    end
   
    for k = 1:44
         hy(1,k) = hy(1,k) + .5*(ex(1,k) - ex(1,k+1));  
    end
 
    plot(SAR);grid
    
    axis([0 m -0.001 0.001]); 
    line([k1,k1],[-.001,.001],'Color','r');
     
  
    frame = getframe();
    
    
end
