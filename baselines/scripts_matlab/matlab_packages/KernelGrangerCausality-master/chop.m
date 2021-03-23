function c=chop(v,val);
if nargin <2
    val=1.e-10;
end
ind=find(abs(imag(v))<val);
v(ind)=real(v(ind));
c=v;
end