clc
clear all
% First Window: Read cost function, max/min, number of constraints
options.Resize = 'on';
title = 'Enter your LP data';
prompt = {'Cost function (c vector)',...
    'Max (enter 1) Min (enter 2)',...
    'Number of constraints'};
DefaultText = {'[ ]','',''};
a = inputdlg(prompt,title,1,DefaultText,options);

c    = eval(a{1});
type = eval(a{2});
nConstraints = eval(a{3});
struct1 = struct('vari',{},'Type',{});
%% Second window: Sense of each inequality constraint
title = '(In)equality Constraints Data';

for i = 1:nConstraints
    prompt = {['Specify (<=,>=,=) for constraint no. ', num2str(i)]};
    DefaultText = {''};
    a = inputdlg(prompt,title,1,DefaultText,options);
    struct1(1,i).Type = a{1};
end

%% Third window: matrix A
title = 'Inequality constraint coefficients';
prompt = {'GIVE ME MATLAB MATRIX A'};
DefaultText = {'[ ]'}; % MATLAB MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
A = eval(a{1});
%% Fourth window: vector b
title = 'Inequality (upper/lower) Bounds';
prompt = {'GIVE ME VECTOR b'};
DefaultText = {'[ ]'}; % VECTOR MATRIX
a = inputdlg(prompt,title,1,DefaultText,options);
b = eval(a{1});
%% Slack, Surplus, Non-basic (activity)
v_e = []; %slack
A1 = []; %matrix for slacks
v_ari = []; %artificial or surplus
A2 = []; %matrix for surplus
v_b = []; % vector of bounds

for i = 1:nConstraints
    n = struct1(1,i).Type;
    
    switch n
        case '<='
            v_e = [v_e b(i)];
            A1(i,length(v_e)) = 1;
            
        case '>='
            v_e = [v_e 0];
            A1(i,length(v_e)) = -1;
            v_ari =[v_ari b(i)];
            A2(i,length(v_ari)) = 1;
            
        case '='
            v_ari =[v_ari b(i)];
            A2(i,length(v_ari)) = 1;
    end
    v_b=[v_b b(i)];
end
% try
% A =[A A1 A2];
% catch
%    if size(A1,1) < size(A2,1)
%       A = [A [A1;zeros(size(A2,1)-size(A1,1),size(A1,2))] A2]; 
%    end
% end
A = [A A1 A2]; 
vari = [];
vari_ar = []; %artificials (surplus)
vari_e = []; %slacks
vari_a = []; % non-basic (acivity)

for i = 1:size(A,2)
    struct1(1,i).vari = ['x',num2str(i)];
    vari = [vari, struct1(1,i).vari,' '];
    if i <= length(c)
        vari_a=[vari_a, struct1(1,i).vari,' ']; % activity variables x1 , x2
    elseif i <= length(c) + length(v_e)
        vari_e =[vari_e, struct1(1,i).vari,' ']; % slack
    else
        vari_ar=[vari_ar, struct1(1,i).vari,' ']; %artificial
    end
end
v_a = 0*c;

x = [v_a,v_e,v_ari];
if(isempty(v_ari))
    v_ar = [];
else
    if type == 1 %max
        v_ar = -9999*ones(1,length(v_ari));
    else
        v_ar = +9999*ones(1,length(v_ari));
    end
end
%% Initialize Cj,Zj


Cj = [c,0.*v_e,v_ar];
Ci = [];
tabl = [];
Vb = [];
Q = v_b;

for i = 1:length(Q)
    tabl = [tabl; ' | '];
    struct2(1,i).value = Q(i);
    ind = find(x == Q(i));
    struct2(1,i).var_base = struct1(1,ind).vari;
    Vb = [Vb,struct2(1,i).var_base,' '];
    Ci = [Ci, Cj(ind)];
end
Z = sum(Ci.*Q);
for i = 1:length(Cj)
    Zj(i) = sum(Ci'.*A(:,i));
end


vars_at_moment=[];
for i = 1:nConstraints
    if(length(struct2(1,i).var_base) == 2)
        vars_at_moment=[vars_at_moment;struct2(1,i).var_base,' '];
    else
        vars_at_moment=[vars_at_moment;struct2(1,i).var_base];
    end
    
end


fprintf('\n')

disp('==== LP in Standard Form =====')
disp(['Variables: ', vari])
disp(['     -Activity (non-basic) Variables: ', vari_a])
disp(['     -Slack Variables: ', vari_e])
disp(['     -Artificial (surplus) Variables: ', vari_ar])

disp('==== Iteration 0 =====')
disp([' Initializing my variables: ', vari])
disp(['     -Activity (non-basic) Variables: ', num2str(v_a)])
disp(['     -Slack Variables: ', num2str(v_e)])
disp(['     -Artificial (surplus) Variables: ', num2str(v_ar)])
disp('===================')
disp(['Cj: ',num2str(Cj)])
disp([tabl,num2str(Ci'),tabl,vars_at_moment,tabl,num2str(Q'),tabl,num2str(A),tabl])
disp('===================')
disp(['Zj: ',num2str(Zj)])
disp(['Cj-Zj: ',num2str(Cj-Zj)])
disp(['Z: ',num2str(Z)])

iterNum = 1;

while iterNum <= 10
    % entering variable
    %max
    if type == 1
        num = max(Cj - Zj);
        num = num(1);
        num1 = find(Cj - Zj == num);
        num1 = num1(1);
        V_enter = struct1(1,num1);
    else
        %min
        num = min(Cj - Zj);
        num = num(1);
        num1 = find(Cj - Zj == num);
        num1 = num1(1);
        V_enter = struct1(1,num1);
    end
    
    b = A(:,num1);
    k=0;
    d=1e4;
    for i = 1:length(Q)
       if(b(i)>0)
          div = Q(i)/b(i);
          if(d > div)
             d = div;
             k = i;
          end
       end
    end
    if k==0
        disp('Solution is infinity');
        break;
    else
        num2 = k;
    end
    
    
    V_leave = struct2(1,num2).var_base;
    struct2(1,num2).var_base = struct1(1,num1).vari;
    pivot = A(num2,num1);
    Ci(num2) = Cj(num1);
    A(num2,:) = A(num2,:)./pivot;
    Q(num2) = Q(num2)./pivot;
    h = size(A,1);
    for i = 1:h
        if i ~= num2
           Q(i) = Q(i) - A(i,num1)*Q(num2);
           A(i,:) = A(i,:) - A(i,num1)*A(num2,:);
        end
    end
    
    Z = sum(Ci.*Q);
    for i = 1:size(A,2)
       Zj(i) = sum(Ci'.*A(:,i)); 
    end
    
    vars_at_moment=[];
    for i = 1:nConstraints
        if(length(struct2(1,i).var_base) == 2)
            vars_at_moment=[vars_at_moment;struct2(1,i).var_base,' '];
        else
            vars_at_moment=[vars_at_moment;struct2(1,i).var_base];
        end

    end
    
    disp(['==== Iteration', num2str(iterNum) ,'=====']);
    disp(['ENTER: ', V_enter.vari]);
    disp(['LEAVE: ', V_leave]);
    disp(['PIVOT: ', num2str(pivot)]);
    
    disp('===================')
    disp(['Cj: ',num2str(Cj)])
    disp([tabl,num2str(Ci'),tabl,vars_at_moment,tabl,num2str(Q'),tabl,num2str(A),tabl])
    disp('===================')
    disp(['Zj: ',num2str(Zj)])
    disp(['Cj-Zj: ',num2str(Cj-Zj)])
    disp(['Z: ',num2str(Z)])
    
    % stopping criterion
    if type == 1
       temp = max(Cj - Zj);
       temp = temp(1);
       if (temp <= 0)
          break; 
       end
    else
        temp = min(Cj - Zj);
        temp = temp(1);
       if (temp >= 0)
          break; 
       end 
    end
    
    iterNum=iterNum+1;
end

