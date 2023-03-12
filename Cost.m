function sol=Cost(sol,data)
   
sol.x=CB(sol.x,data.lb,data.ub);

sol.fit=sum(sol.x.^2);

end

function x=CB(x,lb,ub)

x=max(x,lb);
x=min(x,ub);

end
