% Function to solve Diagonally dominant linear algebraic system using Guass Sedial iterative algorithm:
% Here in GSM-> A: Coefficient Matrix.
%               b: Constants' vector
% NOTE: In this method to ensure convergence the coefficient matrix must be Daigonally dominant.
% Author: Kausadikar Varad Prashant


function Y = Gauss_sedial_fun(A,b)
  %% Argument validations begins
           % Test if A is sqaure matrix
            [numRows,numCols]=size(A);
            if (numRows~=numCols)
                error('Matrix A is not a square matrix')
            end            
            [numRows_b,numCols_b]=size(b);
         % Test if A and b both have n rows
            if ~isequal(numRows,numRows_b)
                error('Matrix A and vector b must have same number of rows')
            end 
        % Test if A is daigonally dominant sqaure matrix
            flag=0;
          for i=1:numRows
              x=0;
              for j=1:numRows
                  if i==j
                      continue
                  end
                  x= x + abs(A(i,j));
              end
              if ( abs(A(i,i))> x)
                  flag=flag+2;
              end
              if ( abs(A(i,i))== x)
                  flag=flag+1;
              end
          end
          if(flag <(numRows+1))
               error('Matrix A is not a Daigonaly dominant and thus this method cannot be applied');
          end
 %% Argument validations Ends and the function code begins:
    
    n=numRows;

    % Intial guess: (n x 1) vector:
    x_io=zeros(n,1);  

    %% Intializing x_f vector:
    x_f=zeros(n,1);

    x_i=x_io;
    % First Iteration:
    for i=1:n
        x_f(i,1)=b(i,1);
        for j=1:n
            if (i==j)
                continue
            end
            x_f(i,1)=x_f(i,1)-A(i,j)*x_i(j,1);
        end
        x_i(i,1)=x_f(i,1)/A(i,i);
    end

    x_f=x_i;
    % Convergence limit:
    con_limit=abs(x_f-x_io);

    %% Iteration in loop:
    flag=1;

    while (con_limit>= 10^-6)
        x_io=x_f;
        for i=1:n
            x_f(i,1)=b(i,1);
            for j=1:n
                if (i==j)
                    continue
                end
                x_f(i,1)=x_f(i,1)-A(i,j)*x_i(j,1);
            end
            x_i(i,1)=x_f(i,1)/A(i,i);
        end
    
        con_limit = abs(x_f-x_io); % Convergence limit
    
        flag=flag+1;
    end

    %% Final answer:
    Y=x_i;
end

