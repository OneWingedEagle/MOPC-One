

function FEM1D
  tic
clear


EPS0=8.854E-12;
XMU0=1.2566306E-6;
T=0;

global Nelm;
global Nnodes;

cryst_const=1;
x0=0;
x1=20;
Nelm=100;
dx=(x1-x0)/Nelm;

eps1=1+0*1j;
eps2=2+0*1j;



Nnodes=Nelm+1;


global indices;
indices=zeros(Nnodes,1);

ix=0;


for i=1:Nnodes
if(i==1 || i==Nnodes)
indices(i)=0;
 continue;
 end;
ix=ix+1;
indices(i)=ix;
end

global xcoord;
xcoord=linspace(x0,x1,Nnodes);
global elms;
elms=zeros(Nelm,2);
global epss;
epss=zeros(Nelm,1);

for i=1:Nelm
  elms(i,1)=i;
  elms(i,2)=i+1;
end

%c=1/sqrt(EPS0*XMU0);
%lam=0.5e-6;
wn=1;
k1=2*pi*wn/cryst_const;
f=k1;
lam=1./f;

w2c2=k1^2;
w2c2

for ie=1:Nelm
  
  if(ie<Nelm/4 || ie>Nelm*3/4)
  epss(ie)=eps1;
  else
 epss(ie)=eps2;
 end

end

global neq;
neq=0;
for i=1:Nnodes
if(indices(i)>0)
 neq=neq+1;
 end;
end


global K;

SetK_Mat(w2c2);


global b;
b=zeros(neq,1);

SetRHS(k1);

sol=K\b;

for n=1:Nnodes
 
  index= indices(n);
  if(index==0) continue; end;
   x=xcoord(n);
  sol(index)= sol(index)+cos(k1*x)+1i*sin(k1*x);
  Ezr(n)=real(sol(index));
  Ezm(n)=imag(sol(index));
  Ez(n)= sqrt(Ezr(n)^2+ Ezm(n)^2);
end

Ezr=zeros(Nnodes,1);
Ezm=zeros(Nnodes,1);
Ez=zeros(Nnodes,1);
for n=1:Nnodes
 
  index= indices(n);
  if(index==0) continue; end;
  x=xcoord(n);
  Ezr(n)=real(sol(index));
  Ezm(n)=imag(sol(index));
  Ez(n)= sqrt(Ezr(n)^2+ Ezm(n)^2);
end

 colorRe='r';
 colorIm='b';
 
    figure(1)
    
    plot(xcoord,Ezr,colorRe);
     hold on
     plot(xcoord,Ezm,colorIm);
      hold on
     plot(xcoord,Ez,'k');
    hold on

toc

end

function SetK_Mat(w2c2)
global K;
global neq;
global Nelm;
global elms;
global indices;
global epss;
K=zeros(neq,neq);
S=zeros(neq,neq);
for ie=1:Nelm
    gg=grdNi_gradNj(ie);
    
    NiNj=Ni_dot_Nj(ie);

    w2c2eps=epss(ie)*w2c2;
   for j=1:2
    n1=elms(ie,j); 
    indx1=indices(n1);
     if(indx1>0)
      for k=1:2
       n2=elms(ie,k);
       indx2=indices(n2);
       if(indx2>0)
         K(indx1,indx2)=  K(indx1,indx2)+gg(j,k);
         S(indx1,indx2)=  S(indx1,indx2)+w2c2eps*NiNj(j,k);
       end
      end
     end
    end
  end

K=K+S;
end



function [gg]=grdNi_gradNj(ie)
gg=zeros(2,2);
global elms;
global xcoord;

  n1=elms(ie,1);
  n2=elms(ie,2);
  dx=xcoord(n2)-xcoord(n1);
  g=1./dx;
  
  gg(1,1)=g*g*dx;
  gg(1,2)=-g*g*dx;
  gg(2,1)=-g*g*dx;
  gg(2,2)=g*g*dx;
  
end

function [NiNj]=Ni_dot_Nj(ie)
NiNj=zeros(2,2);
global elms;
global xcoord;

  n1=elms(ie,1);
  n2=elms(ie,2);
  dx=xcoord(n2)-xcoord(n1);

  NiNj(1,1)=dx/2;
  NiNj(1,2)=dx/6;
  NiNj(2,1)=dx/6;
  NiNj(2,2)=dx/2;

  
end

function SetRHS(k)
global b;
global elms;
global xcoord;
global Nelm;
global indices;

k2=k*k;

for ie=1:Nelm
  
    n1=elms(ie,1); 
    n2=elms(ie,2);
    
    x1=xcoord(n1);
    x2=xcoord(n2);
    
    d=x2-x1;
    
    br1=(cos(k*(d+x1))-cos(k*x1))/(d*k2)-sin(k*x1)/k;
    br2=-(cos(k*(d+x1))-cos(k*x1))/(d*k2)+sin(k*(x1+d))/k;
    
    bm1=(sin(k*(d+x1))-sin(k*x1))/(d*k2)+cos(k*x1)/k;
    bm2=-(sin(k*(d+x1))-sin(k*x1))/(d*k2)-cos(k*(x1+d))/k;
    
     index1=indices(n1);
     index2=indices(n2);
     
     if(index1>0)
     b(index1)= b(index1)-br1*(1-k2)-1i*bm1*(1-k2);
     end;
     if(index2>0)
     b(index2)= b(index2)-br2*(1-k2)-1i*bm2*(1-k2);
     end;     

end

  
end