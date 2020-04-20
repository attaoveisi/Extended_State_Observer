function y=PBH(A,B)
%PBH test
n=size(B,1);
Lambda=eig(A);
j=0;
for i=1:length(Lambda)
    if rank([Lambda(i)*eye(n)-A B])<n
        j=j+1;
        break
    end
end
if j==0
    y='Controllable';
else
    y='Uncontrollable';
end