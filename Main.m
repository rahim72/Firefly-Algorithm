%%
%% The FireFly Algorithm
%%
clc
clear
close all
 format shortG

%% Insert Data
%%
data=[];

%% Parameters Definiterion
%%
nvar=5;                           % Number of  Variables
lb=-10*ones(1,nvar);             %  Variables Lower Bound
ub= 10*ones(1,nvar);             %  Variables Upper Bound
%% Menu
%%
prompt={ ' Maximum Number of Iterations','Number of FireFly', 'Attraction Coefficient Base Value', ' Mutation Coefficient1'...
    ,'Radius Reduction Factor' };        
title='Enter a ITER';
nline=([1 40;1 40;1 40;1 40;1 40]);
default={ '200','70','2','0.02','0.96'};
ev=inputdlg(prompt,title,nline,default);
ITER=ev(1,:);ITER=str2num(ITER{:});
cc=ev(2,:);cc=str2num(cc{:});
aa=ev(3,:);aa=str2num(aa{:});
bb=ev(4,:);bb=str2num(bb{:});
dd=ev(5,:);dd=str2num(dd{:});
%%
%%
maxiter=ITER; % Maximum Number of iterations
n=cc;              % Number of Fireflies  n
L=1;
gamma=1./sqrt(L);            % Light Absorption Coefficient
beta0=bb;                     % Attraction Coefficient Base Value
alpha=aa;                   % Mutation Coefficient
alpha_RF=dd;                %Radius Reduction Factor 
data.lb=lb;
data.ub=ub;
%% Create Random Pop
tic
emp.x=[];
emp.fit=[];
emp.info=[];
pop=repmat(emp,n,1);
for i=1:n
   pop(i).x=unifrnd(lb,ub);
   pop(i)=Cost(pop(i),data);
end
%% Best Solution
%%
[~,ind]=min([pop.fit]);
gpop=pop(ind);
%% Main Loop
%%
BEST=zeros(maxiter,1);
MEAN=zeros(maxiter,1);
for iter=1:maxiter
    newpop=pop;
    k=n;  
    for i=1:n
        for j=1:n
            if pop(j).fit<=pop(i).fit
                k=k+1;
                                     % move firefly i towards j in all dimensions
                rij=norm(pop(i).x-pop(j).x,2); % The distance between fireflies i and j (rij) is evaluated as  the Euclidean distance presented in 
                                                                
                beta=beta0*exp(-gamma*rij^2);  % Beta 
                E=alpha*(unifrnd(-1,1,1,nvar).*(ub-lb));  
                                                           
                newpop(i).x=pop(i).x... 
                            +beta*(pop(j).x-pop(i).x)...
                            +E;         
                newpop(i)=Cost(newpop(i),data);
                newpop(k)=newpop(i);
            end
        end
    end
    % Merge
    [pop]=[pop;newpop;gpop];
    % Sort and Select
    [~, ind]=sort([pop.fit]);
    pop=pop(ind);
    pop=pop(1:n);
   % Select Best Sol 
    gpop=pop(1);     
    BEST(iter)=gpop.fit;
    MEAN(iter)=mean([pop.fit]);
    disp('----------------------------------------------------');
    disp(['iter ' num2str(iter) ' Best= ' num2str(BEST(iter))]); 
    % Reduction Mutation Coefficient
    alpha=alpha*alpha_RF;
end

%% Results
%%
disp(' ')
disp('=================================')
disp([ ' BEST solution = '  num2str(gpop.x)]);
disp([ ' BEST fitness = '  num2str(gpop.fit)]);
disp([ ' Time = '  num2str(toc)]);
disp('=================================')
figure()
plot(BEST,'r','LineWidth',1)
hold on
plot(MEAN,'b','LineWidth',1)
xlabel(' Iteration ')
ylabel(' Fitness')
legend( 'BEST','MEAN')
% title('Firefly Algorithm');

hold off