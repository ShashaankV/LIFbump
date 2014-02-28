%edited on 02-27-14 by Shashaank Vattikuti

function pv=popvec(f)

[N,ibins]=size(f);
theta=zeros(N,3);
theta(:,1)=360*[1:N]/N; %assign degree coordinate to neuron
theta(:,2)=cosd(theta(:,1)); 
theta(:,3)=sind(theta(:,1));
max_mag=sum(f); 

for j=1:ibins %go through time bins
    for i=1:N
        mag=f(i,j)/max_mag(j); %normalize firing rate by max value in time block
        rcos(i)=theta(i,2).*mag;
        rsin(i)=theta(i,3).*mag;
     end
     pv_rcos=sum(rcos); 
     pv_rsin=sum(rsin);
     PV=atand(pv_rsin/pv_rcos);
     if pv_rsin>0 & pv_rcos>0
         pv(j)=PV;
     elseif pv_rsin<0 & pv_rcos>0
         pv(j)=360+PV;
     else
         pv(j)=180+PV;
     end
 end
end