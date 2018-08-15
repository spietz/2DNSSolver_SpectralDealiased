function h=Lagrange_poly(i,x,nodes)

h=1;
N=length(nodes)-1;
for j=1:N+1
 if(j-1~=i)   
  term=(x-nodes(j))/(nodes(i+1)-nodes(j));   
  h=h.*term;
 end
end
