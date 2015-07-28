% This script is the main program.
clear all
global dt Nb N Nw h rho mu ip im a lambda;
global kp km dtheta dWin dBone K Xo Zo w0 s0 amp dir freq;

tic;
initialize
init_a

BMposition= zeros(Nb,2,clockmax);
BM_displacement=zeros(Nb,2,clockmax);
Bone_position= zeros(length(Zo),2,clockmax);
Bone_displacement= zeros(length(Zo),2,clockmax);
velocity= zeros(N,N,2,floor(clockmax/1000));
vorticity= zeros(N,N,floor(clockmax/1000));
i=1;

for clock=1:clockmax
 
  XX= X+(dt/2)*interp(u,X); %membrane
  ZZ= Z+(dt/2)*interpB(u,Z); %bone
  ff= spreadWin(f_source(clock*dt),w0)+spread(Force(XX),XX)+spreadB(BForce(ZZ),ZZ);
  
  [u,uu]=fluid(u,ff);
 
  X= X+dt*interp(uu,XX);
  Z= Z+dt*interpB(uu,ZZ);
  
  if mod(clock,1000)==0
      velocity(:,:,:,i)=u;
      vorticity(:,:,i)=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
      i=i+1;
  end
  
  BMposition(:,:,clock)=X;
  BM_displacement(:,:,clock)=X-Xo;
  Bone_position(:,:,clock)=Z;
  Bone_displacement(:,:,clock)=Z-Zo;
  
  %animation:
  %{
  plot(Xo(:,1),BM_displacement(:,2,clock))
  axis([0,7,-alpha,alpha])
  drawnow
  %}
  display(clock)
end

elapsed_time=clock;
save('data')
