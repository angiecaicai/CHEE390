function [ r,A,B ] = solvroot(am,bm,P,T )

A=(am*P)/(8.314^2*T^2);
B=(bm*P)/(8.314*T);

fun=[1,-(1-B),(A-3*(B^2)-2*B),-(A*B-B^2-B^3)];
rtemp=roots(fun);
%r should be a 3x1 array

r=zeros(1,3);
j=1;
for i=1:3
    if isreal(rtemp(i))==1 && rtemp(i)>=0
        r(j)=rtemp(i);
        j=j+1;
    end
end

r(j:end)=[];
    
    

end

