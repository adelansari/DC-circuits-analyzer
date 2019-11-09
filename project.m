%%   Engineering Computation & Linear Algebra   %%
%*   Project                                    *%
%*   Adel Ali Ansari                            *%
%*   U00038673                                  *%
%*   Output file: results                       *%
%************************************************%
clear all
clc

%%  The input code

file= fopen('results','w');
disp('****Engineering Computation & Linear Algabra Project');
disp(sprintf('****Finding the mesh current in purely resistive independent DC sources circuit\n'));
disp(sprintf('\t\tName:\t\tAdel Ali Ansari\n\t\tID:\t\t\tU00038673\n\t\tInstructor:\tDr.Ismail Shahin'));

fprintf(file,'**\tEngineering Computation & Linear Algabra Project\n');
fprintf(file,'**\tFinding the mesh current in purely resistive independent DC sources circuit\n');
fprintf(file,'\t\tName:\t\tAdel Ali Ansari\n\t\tID:\t\t\tU00038673\n\t\tInstructor:\tDr.Ismail Shahin');
fprintf(file,'\n===========================================================================\n');

fprintf(file,'\t\t\t\t\t\t\t\t{INPUTS}\n');

%   The number of meshes
m=input('\n\nEnter the number of mesh currents in the circuit (max. = 10): ');

    %  Checking the number of meshes
while m>10 || m<=0

    error1= sprintf('\nError!!! The number of mesh currents is more than 10');
    error2= sprintf('\nError!!! The number of mesh currents can''t be zero');
    error3= sprintf('\nError!!! The number of mesh currents can''t be negative');
    if m>10
        disp(error1);
        m=input('\n\nPlease enter the number of mesh currents in the circuit again: ');
    else if m==0
            disp(error2);
            m=input('\n\nPlease enter the number of mesh currents in the circuit again: ');
        else
            disp(error3);
            m=input('\n\nPlease enter the number of mesh currents in the circuit again: ');
        end
    end
end

%   Sending the meshes number to output file
fprintf(file,'\n\nThe number of mesh currents = %d\n',m);

%   Entering the resistance values
for i=1:m
    R(i,i)= input(sprintf('\nEnter the equivalent resistance of the mesh (%d) in Kohm = ',i));
    %   Checking the values:
    while R(i,i)<0
        disp(sprintf('\nError!!! The equivalent resistance can''t be negative'));
        R(i,i)= input(sprintf('Please enter the equivalent resistance of the mesh (%d) in Kohm again = ',i));
    end
end

for i=1:m
    for j=i+1:m
        R(i,j)= input(sprintf('\nEnter the shared resistance between the meshes (%d) and (%d) in Kohm = ',i,j));
        %   Checking the values:
        while R(i,j)<0
            disp(sprintf('\nError!!! The shared resistance can''t be negative'));
            R(i,j)= input(sprintf('Please enter the shared resistance of the meshes (%d) and (%d) in Kohm again = ',i,j));
        end
        R(i,j)=-R(i,j);
        R(j,i)=R(i,j);
    end
end


%   Sending the resistance matrix to the output file
fprintf(file,'\nR(Kohm)=\n');
fprintf(file, [repmat('\t\t%10.6f', 1, m) '\n'], R');

% Explaination:
% A = (1:3)'
% A =
% 
%      1
%      2
%      3
% Create a horizontal stack of four copies of A.
% 
% B = repmat(A,1,4)
% B =
% 
%      1     1     1     1
%      2     2     2     2
%      3     3     3     3
% B contains 1 copy of A in the first dimension and 4 copies in the second dimension

%*********************************************************
% % Alternative Method:                                 **
% % Sending the resistance matrix to the output file    **
%   fprintf(file,'\nR(Kohm)=\n');                       **
%   for i=1:m                                           **
%       for j=1:m                                       **
%           if j==1                                     ** 
%               fprintf(file,'\t\t%10.6f',R(i,j));      **
%           else                                        **
%               fprintf(file,'\t%10.6f',R(i,j));        **
%           end                                         **     
%       end                                             **
%   fprintf(file,'\n');                                 **
%   end                                                 **
%*********************************************************


%   Entering the voltage values
for i=1:m
    V(i,1)=input(sprintf('\nEnter the algebraic sum of the voltage in mesh (%d) in volt [Please take care of the sign]= ',i));
end

%   Sending the voltage vector to the output file
fprintf(file,'\nV(volt)=\n');
fprintf(file,'\t\t%10.6f\n',V);

%   Entering the value of the tolerance and sending it to the output file
tol=input(sprintf('\n\nEnter the tolerance for the iterative method in mA = '));
    %   Checking The values
while tol<=0 || tol>=1
    disp(sprintf('\nError!!! The tolerance is either very small, negative or very large'));
    tol=input(sprintf('Please enter the tolerance in mA again = '));
end
fprintf(file,'\nThe tolerance in (mA) = %8.6f\n',tol);

%   Entering the initial values of the currents and seding it to the output
for i=1:m
    I_0(i,1)=input(sprintf('\nEnter the initial value for the mesh current (%d) in mA = ',i));
end
fprintf(file,'\nThe initial values for the current Io(mA)=\n');
fprintf(file,'\t\t%10.6f\n',I_0);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Triangular Factorization Method       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(file,'\n===========================================================================\n');
fprintf(file,'\t\t\t\t\t{Triangular Factorization Method}\n');


%   Test

C=abs(R);
solution=1;
for j=1:m
    b = max(C(j:m,j));
    p = find(C(:,j)==b);
    if C(p,j)==0
        disp(sprintf('\n\t**Error!!! The maximum element in column (%d) is zero.\n',j));
        disp(sprintf('\t  The matrix is singular. There''s no exact solution.\n'));
        fprintf(file,'\n\n\t\t***Error!!! There''s no exact solution for this circuit');
        solution=0;
        break
    end
end

if solution==1
    fprintf(file,'\n\n**There''s an exact solution.\n');
    disp(sprintf('\n\n**There''s an exact solution since matrix R is a nonsingular matrix.\n'));
end


%   Initialize L to the identity matrix:
L = eye(m);

%   Initialize U to the R matrix:
U = R;

for j=1:m-1
    for i=j+1:m
        L(i,j)= U(i,j)/U(j,j);
        U(i,:)=U(i,:)-U(i,j)/U(j,j)*U(j,:);
        
    end
end


%   Finding the matrix Y using forward substitution
Y= zeros(m,1);   %init vector to zeros
Y(1)= V(1);
for i=2:m
    Y(i)=V(i)-L(i,:)*Y;
end

%   Finding the Mesh Currents I using backward substitution
I= zeros(m,1);
I(m)=Y(m)/U(m,m);
for i=m-1:-1:1
    I(i)=(Y(i)-U(i,:)*I)/U(i,i);
end


%   Sending the Lower Triangular matrix (L) to the output file
fprintf(file,'\nThe Lower Triangular matrix (L):\n');
fprintf(file,'L=\n');
fprintf(file, [repmat('\t\t%10.6f', 1, m) '\n'], L');

%   Sending the Upper Triangular matrix (U) to the output file
fprintf(file,'\nThe Upper Triangular matrix (U):\n');
fprintf(file,'U=\n');
fprintf(file, [repmat('\t\t%10.6f', 1, m) '\n'], U');

%   Sending the Y Matrix to the output file
fprintf(file,'\nThe Y matrix :\n');
fprintf(file,'Y =\n');
for i=1:m
    fprintf(file,'\t\t%10.6f\n',Y(i));
end


%   Printing the Mesh currents, I
fprintf(file,'\nThe mesh currents, I :\n');
fprintf(file,'I =\n');
for i=1:m
    fprintf(file,'\t\t%10.6f\n',I(i));
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Jacobi Iterative method     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf(file,'\n===========================================================================\n');
fprintf(file,'\t\t\t\t\t\t\t{Jacobi Iterative Method}\n');


%  Testing if there is iterative solution
AbsR= abs(R);
solution2=1;
for i=1:m
    if AbsR(i,i) < (sum(AbsR(i,:))-AbsR(i,i))
        disp(sprintf('\n\t**Error!The diagonal magnitude in row(%d) is less the sumation of all magnitude of the row.\n\t  There''s no iterative solution.\n\t  End of the programme.',r));
        fprintf(file,'\n\n\t\t***Error!There''s no iterative solution for this circuit\n\n');
        solution2=0;
        return      %   end of the programme
    end
end

if solution2==1
    fprintf(file,'\n\n\t\t**There''s an iterative solution for this circuit.\n\n');
    disp(sprintf('\n**There''s an iterative solution for this circuit.\n'));
end


I_Jac(:,1)=I_0;

% Initializing I(k+1) to the initial guess
I = I_0;
I_old = I_0;        
% Initializing the Norm with a value greater than the tolerance (avoiding the loop to be aborted from the first iteration)
n=100;
        
% Counting the number of iterations (k)
k = 0;


while n >= tol

    k = k + 1;
    for i=1:m
        segma = 0;
        for j=1:m
            if j ~= i
                segma = segma + R(i,j)*I_old(j);
            end
        end
        I(i) = (1/R(i,i))*(V(i)-segma);
    end

    % Calculating the infinite norm between I(k+1) & I(k)
    n = max(abs(I-I_old));    
    I_Jac(:,k+1)=I(:,1);
    I_old = I;

end

%  Sending the number of iterations to the output file
fprintf(file,'\nThe number of iterations done (K): %d\n\n',k);
fprintf(file,'The iterations results are in the table.(All the currents and the Norm are in (mA))\n\n');

for i=1:k
    Norm(i+1)=norm(I_Jac(:,i+1)-I_Jac(:,i),inf);
end

fprintf(file,'K\t\t\t');
for i=1:m
    fprintf(file,'I %d\t\t\t\t',i);
end
fprintf(file,'Norm\n');

for i=1:m+2
    if i==1
        fprintf(file,'_______');
    else
        fprintf(file,'________________');
    end
end

for i=0:k
    fprintf(file,'\n%d\t\t',i);
    for j=1:m
        fprintf(file,'%10.6f\t\t',I_Jac(j,i+1));
    end
    fprintf(file,'%10.6f\t\t',Norm(i+1));
end
fprintf(file,'\n\n');

%  Ploting the currents & the norm
%  The legend of the plot
for i=1:m
    name(i,:)=['I ' num2str(i) ' ' ];
end


name(m+1,:) = 'Norm' ;

figure
plot(0:1:k,[I_Jac ; Norm]);
legend(name,'Location','NorthEastOutside');
title('Jacobi Iterative Method');
xlabel('Iterations');
ylabel('Currents & Norm (mA)');
grid;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Gauss-Seidel Iterative Method      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf(file,'\n===========================================================================\n');
fprintf(file,'\t\t\t\t\t\t{Gauss-Seidel Iterative Method}\n');



I_GS(:,1)=I_0;

% Initializing I(k+1) to the initial guess
I = I_0;
I_old = I_0;        
% Initializing the Norm with a value greater than the tolerance (avoiding the loop to be aborted from the first iteration)
n2=100;
        
% Counting the number of iterations (k)
L = 0;


while n2 >= tol

    L = L + 1;
    
    for i=1:m
        segma_1 = 0;
        for j=1:(i-1)
            segma_1 = segma_1 + R(i,j)*I(j);
        end
        
        segma_2 = 0;
        for j=(i+1):m
            segma_2 = segma_2 + R(i,j)*I_old(j);
        end
        
        I(i) = (1/R(i,i))*(V(i)-segma_1-segma_2);
        
    end

    % Calculating the infinite norm between I(k+1) & I(k)
    n2 = max(abs(I-I_old));    
    I_GS(:,L+1)=I(:,1);
    I_old = I;

end

%  Sending the number of iterations to the output file
fprintf(file,'\nThe number of iterations done (K): %d\n\n',L);
fprintf(file,'The iterations results are in the table.(All the currents and the Norm are in (mA))\n\n');

for i=1:L
    Norm2(i+1)=norm(I_GS(:,i+1)-I_GS(:,i),inf);
end

fprintf(file,'K\t\t\t');
for i=1:m
    fprintf(file,'I %d\t\t\t\t',i);
end
fprintf(file,'Norm\n');

for i=1:m+2
    if i==1
        fprintf(file,'_______');
    else
        fprintf(file,'________________');
    end
end

for i=0:L
    fprintf(file,'\n%d\t\t',i);
    for j=1:m
        fprintf(file,'%10.6f\t\t',I_GS(j,i+1));
    end
    fprintf(file,'%10.6f\t\t',Norm2(i+1));
end
fprintf(file,'\n\n');

%  Ploting the currents & the norm
%  The legend of the plot
for i=1:m
    name2(i,:)=['I ' num2str(i) ' ' ];
end


name2(m+1,:) = 'Norm' ;

figure
plot(0:1:L,[I_GS ; Norm2]);
legend(name2,'Location','NorthEastOutside');
title('Gauss-Seidel Iterative Method');
xlabel('Iterations');
ylabel('Currents & Norm (mA)');
grid;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Message for the user
disp(sprintf('\n***The results are in an output file called "results"***\n\t\t\t\t\tEND OF THE PROGRAMME'));

fclose(file);
%  End of the programme
