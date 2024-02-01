function [pos_array,hat_x]=cs_somp(y,T_Mat,t,L) 
[m,n]=size(y);
s=L;  
hat_x=zeros(t,n);                   
Aug_t=[];        
r_n=y;   

for times=1:s;  
    pro=T_Mat'*r_n;
    %[~,g]=size(T_Mat);
    for i=1:t
       product(i)=sum(abs(pro(i,:)));
    end
   
    [val,pos]=max(product); 
    
    Aug_t=[Aug_t,T_Mat(:,pos)];   
    T_Mat(:,pos)=zeros(m,1);  
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;   
    r_n=y-Aug_t*aug_x;    
    norm(r_n);
    pos_array(times)=pos;   
    
end

hat_x(pos_array,:)=aug_x;  %  重构的向量 


end