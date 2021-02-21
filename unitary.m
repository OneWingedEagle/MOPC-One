

y0=-9.36+1.17/4;
for i=0:16
y=y0+i*1.17  
 end
unitar=1;
if(unitar==1)
eps=5.5;
gamab=.01;

e=[eps 0 -1i*gamab;0 eps 0;1i*gamab 0 eps ];
U=1./sqrt(2)*[1 0 1i;0 sqrt(2) 0;1 0 -1i ];

Uc=inv(U)

U
Uc



D=U*e*Uc

e2=Uc*D*U

nn=sqrt(D)
end


