
unitar=1;
if(unitar==1)
epsbx=2.72;
epsby=2.72;
epsbz=2.72;
gamab=.03;

e=[epsbx 0 -1i*gamab;0 epsby 0;1i*gamab 0 epsbz ];
e1=[epsbx -1i*gamab 0;0 i*gamab 0; 0 0 epsbz ];
U=1./sqrt(2)*[1 0 1i;0 sqrt(2) 0;1 0 -1i ];

Uc=inv(U)

U
Uc



D=U*e*Uc

e2=Uc*D*U

nn=sqrt(D)
end


