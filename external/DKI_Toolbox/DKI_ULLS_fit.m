function [DT,KT, S0]=...
DKI_ULLS_fit(y,bval,bvec)

Nvol=length(y);

Z15=zeros(Nvol,15);% old Z
Z6=zeros(Nvol,6);% old Z2
Ad=Z6;
Ak=Z15;

for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
    
    Ak(v,1:15)=[b*b*bvec(1,v)^4,...
        b*b*bvec(2,v)^4, ...
        b*b*bvec(3,v)^4, ...
        4*b*b*bvec(1,v)^3*bvec(2,v),...
        4*b*b*bvec(1,v)^3*bvec(3,v),...
        4*b*b*bvec(2,v)^3*bvec(1,v),...
        4*b*b*bvec(2,v)^3*bvec(3,v),...
        4*b*b*bvec(3,v)^3*bvec(1,v),...
        4*b*b*bvec(3,v)^3*bvec(2,v),...
        6*b*b*bvec(1,v)^2*bvec(2,v)^2,...
        6*b*b*bvec(1,v)^2*bvec(3,v)^2,...
        6*b*b*bvec(2,v)^2*bvec(3,v)^2,...
        12*b*b*bvec(1,v)^2*bvec(2,v)*bvec(3,v),...
        12*b*b*bvec(2,v)^2*bvec(1,v)*bvec(3,v),...
        12*b*b*bvec(3,v)^2*bvec(1,v)*bvec(2,v)];
    
end

%A
A=[ -Ad 1/6*Ak ones(Nvol, 1)]; %SO is now the last column
B=log(y);
piA=pinv(A); 
X=piA*B;
MD=mean(X(1:3));
DT=X(1:6);
KT=X(7:21)/(MD^2);
S0=exp(X(22));
