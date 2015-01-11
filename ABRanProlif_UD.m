% Implementing spatial Moran 2D algorithm: Normal cells=1's
% and PM cells = 0's.
% Input: 21 X 21 matrix with one 0 and rest all 1's
% Generate a randomproliferation matrix for A and B-
% <rB> = 1.5 and <rA> = 1.0
% rB: D = [1.4,1.6] and rA: D = [0.9,1.1]

clc; clear all;

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); % Random seed

% create a grid
xmax = 21;
ymax = 21;
x = 0:1:xmax;
y = 0:1:ymax;

A_won = 0;
B_won = 0;

lowerlimitA = 0.9;
upperlimitA = 1.1;
rA = lowerlimitA + (upperlimitA-lowerlimitA).*rand(21,21);

lowerlimitB = 0.1;
upperlimitB = 2.9;
rB = lowerlimitB + (upperlimitB-lowerlimitB).*rand(21,21);

fid = fopen('RanProlifB-A[0.9,1.1]B[0.1,2.9]1.txt', 'wt');

%-------------Randomly choose (i,j) to place PM cell by 0--------------

sets_max = 1;
iter_max = 50000;
for sets=1:sets_max   % Number of sets
    sets;
    for iterations = 1:iter_max
        iterations;
        fprintf(fid, 'Iterations: %g\n', iterations);
        
        X = ones(xmax,ymax);
        M = sum(X);
        M1 = sum(M);
        M2 = 0;
        M1init = xmax*ymax-1;
        %i1 = randint(1,1,[1,xmax]);
        %j1 = randint(1,1,[1,ymax]);
        
        i1 = randi(xmax,1,1);
        j1 = randi(ymax,1,1);
        
        X(i1,j1) = 0;
        while ((M1 ~= sum(sum(X))) && (M2 ~= sum(sum(X))))
            X;
            %i2 = randint(1,1,[1,xmax]);
            %j2 = randint(1,1,[1,ymax]);
            
            i2 = randi(xmax,1,1);
            j2 = randi(ymax,1,1);
            
            X(i2,j2) = -1;
            
            X_withcelldeath = X;
            
            while (any(any(X==-1))==1)  % two times any since it checks both rows and columns
                X_minusonepresent = X;
                [k l] = ind2sub(size(X),find(X==-1)); % check the position i,j for -1
                i2 = k;
                j2 = l;
                % Interior part (2,2)...(2,4) and (4,2)...(4,4)
                if ((i2>1) && (i2<xmax) && (j2>1) && (j2<ymax))
                    neighbours(1) = X(i2,j2-1); % left side
                    neighbours(2) = X(i2,j2+1); % right side
                    neighbours(3) = X(i2-1,j2); % Top
                    neighbours(4) = X(i2+1,j2) ;% Bottom
                end
                
                % Top part (1,2)...(1,4)
                if (i2==1)
                    
                    if ((j2>1) && (j2<ymax))
                        k = 1;
                        neighbours(1) = X(i2,j2-1);
                        neighbours(2) = X(i2,j2+1);
                        neighbours(3) = 10;
                        neighbours(4) = X(i2+1,j2);
                    end
                    
                end
                
                % Bottom part (5,2)...(5,4)
                if (i2==xmax)
                    if ((j2>1) && (j2<xmax))
                        k = 1;
                        neighbours(1) = X(i2,j2-1);
                        neighbours(2) = X(i2,j2+1);
                        neighbours(3) = X(i2-1,j2);
                        neighbours(4) = 10;
                    end
                    
                end
                
                % Left Side part (2,1)...(4,1)
                if (j2==1)
                    
                    if ((i2>1) && (i2<xmax))
                        k = 1;
                        neighbours(1) = 10;
                        neighbours(2) = X(i2,j2+1);
                        neighbours(3) = X(i2-1,j2);
                        neighbours(4) = X(i2+1,j2);
                    end
                end
                
                % Right Side part (2,5)...(4,5)
                if (j2==ymax)
                    
                    if ((i2>1) && (i2<xmax))
                        k = 1;
                        neighbours(1) = X(i2,j2-1);
                        neighbours(2) = 10;
                        neighbours(3) = X(i2-1,j2);
                        neighbours(4) = X(i2+1,j2);
                    end
                end
                
                if ((i2==1) && (j2== 1))
                    k = 1;
                    neighbours(1) = 10;
                    neighbours(2) = X(i2,j2+1);
                    neighbours(3) = 10;
                    neighbours(4) = X(i2+1,j2);
                end
                
                if ((i2==1) && (j2== ymax))
                    k = 1;
                    neighbours(1) = X(i2,j2-1);
                    neighbours(2) = 10;
                    neighbours(3) = 10;
                    neighbours(4) = X(i2+1,j2);
                end
                
                if ((i2==xmax) && (j2== 1))
                    k = 1;
                    neighbours(1) = 10;
                    neighbours(2) = X(i2,j2+1);
                    neighbours(3) = X(i2-1,j2);
                    neighbours(4) = 10 ;
                end
                
                if ((i2==xmax) && (j2== ymax))
                    k = 1;
                    neighbours(1) = X(i2,j2-1);
                    neighbours(2) = 10;
                    neighbours(3) = X(i2-1,j2);
                    neighbours(4) = 10;
                end
                neighbours;
                
                
                
                % Count the number of normal and pre-malignant cells surrounding
                % the empty spot where cell death happened
                
                a = sum(neighbours);
                b = floor(a/10);
                NCells = a - b*10;
                PMCells = 4 - b - NCells;
                
                %---------Calculate the probabilities----------
                
                nA = NCells;
                nB = PMCells;
                [i2 j2];
                
                %-----Calculate the average r------
                
                if (length(find(neighbours==1) ~=0 ))
                    for r = 1:length(neighbours)
                        if (neighbours(1)==10)
                            rAtemp(1) = 0;
                        else
                            linearindex1 = sub2ind(size(X),i2,j2-1);
                            rAtemp(1) = rA(linearindex1);
                        end
                        if (neighbours(2)==10)
                            rAtemp(2) = 0;
                        else
                            linearindex2 = sub2ind(size(X),i2,j2+1);
                            rAtemp(2) = rA(linearindex2);
                        end
                        if (neighbours(3)==10)
                            rAtemp(3) = 0;
                        else
                            linearindex3 = sub2ind(size(X),i2-1,j2);
                            rAtemp(3) = rA(linearindex3);
                        end
                        if (neighbours(4)==10)
                            rAtemp(4) = 0;
                        else
                            linearindex4 = sub2ind(size(X),i2+1,j2);
                            rAtemp(4) = rA(linearindex4);
                        end
                    end
                    X;
                    rA;
                    rAtemp;
                    %pause
                    rAtemp1 = find(neighbours==1);
                    for t1 = 1:length(rAtemp1)
                        rAtemp2 = rAtemp(rAtemp1);
                        temp3 = rAtemp2;
                    end
                    temp3;
                    rAave = mean(temp3);
                    %pause
                    
                else
                    rAave = 0;
                end
                %------end of average rA code------
                
                rAave;
                
                %-----Calculate the average rB------
                
                if (length(find(neighbours==0) ~=0 ))
                    for r = 1:length(neighbours)
                        if (neighbours(1)==10)
                            rBtemp(1) = 0;
                        else
                            linearindex1 = sub2ind(size(X),i2,j2-1);
                            rBtemp(1) = rB(linearindex1);
                        end
                        if (neighbours(2)==10)
                            rBtemp(2) = 0;
                        else
                            linearindex2 = sub2ind(size(X),i2,j2+1);
                            rBtemp(2) = rB(linearindex2);
                        end
                        if (neighbours(3)==10)
                            rBtemp(3) = 0;
                        else
                            linearindex3 = sub2ind(size(X),i2-1,j2);
                            rBtemp(3) = rB(linearindex3);
                        end
                        if (neighbours(4)==10)
                            rBtemp(4) = 0;
                        else
                            linearindex4 = sub2ind(size(X),i2+1,j2);
                            rBtemp(4) = rB(linearindex4);
                        end
                    end
                    X;
                    rB;
                    rBtemp;
                    %pause
                    rBtemp1 = find(neighbours==0);
                    for t1 = 1:length(rBtemp1)
                        rBtemp2 = rBtemp(rBtemp1);
                        temp3 = rBtemp2;
                    end
                    temp3;
                    rBave = mean(temp3);
                    %pause
                    
                else
                    rBave = 0;
                end
                %------end of average rB code------
                
                rBave;
                mA = 0;
                mB = 0;
                rBave;
                Ktilde = nA*(rAave+mA) + nB*(rBave+mB);
                
                ProbAdiv     = (nA*rAave)/Ktilde;
                ProbAmig     = (nA*mA)/Ktilde;
                ProbBdiv     = (nB*rBave)/Ktilde;
                ProbBmig     = (nB*mB)/Ktilde;
                
                Prob = [ProbAdiv ProbBdiv ProbAmig ProbBmig];
                
                % Generate a random number to see which of the above events is
                % likely to occur
                
                r = rand;
                
                if (r<=Prob(1))
                    %display('Event: A divides and places a cell in empty slot')
                    X(i2,j2) = 1 ;
                    M1updated = sum(sum(X));
                    M1init = M1updated;
                    X_Aafterbirth=X;
                    %pause
                    %return;
                end
                if ((Prob(1)< r) && (r<= Prob(1)+Prob(2)))
                    %display('Event: B divides and places a cell in empty slot')
                    [Q1 Q2] = ind2sub([xmax ymax],find(X==-1)); % check the position i,j for -1
                    Q = [Q1 Q2];
                    X(Q(1),Q(2)) = 0 ;
                    X;
                    Mupdated = sum(X);
                    M1updated = sum(Mupdated);
                    M1init = M1updated;
                    X_Bafterbirth = X;
                    %pause
                    %return;
                end
                if ((Prob(1)+Prob(2)< r) && (r<= Prob(1)+Prob(2)+Prob(3)))
                    %display('Event: A migrates into the empty spot')
                    %return;
                end
                if (Prob(1)+Prob(2)+Prob(3)< r)
                    %display('Event: B migrates into the empty spot')
                    Y = find(neighbours==0);
                    interval = Prob(4)/nB;
                    
                    a1 = r - (Prob(1)+Prob(2)+Prob(3));
                    a2 = ceil(a1/interval);
                    if (Y(a2)==1)
                        X(i2,j2-1) = -1;
                        X(i2,j2) = 0;
                    end
                    if (Y(a2)==2)
                        X(i2,j2+1) = -1;
                        X(i2,j2) = 0;
                    end
                    if (Y(a2)==3)
                        X(i2-1,j2) = -1;
                        X(i2,j2) = 0;
                    end
                    if (Y(a2)==4)
                        X(i2+1,j2) = -1;
                        X(i2,j2) = 0;
                    end
                    X_aftermig = X;
                    %pause
                    
                end
                
                Mupdated = sum(X);
                M1updated = sum(Mupdated);
                M1init = M1updated;
                
            end
            X;
        end
        if (sum(sum(X)) == xmax*ymax)
            A_won = A_won + 1;
        end
        
        if (sum(sum(X)) == 0)
            B_won = B_won + 1;
            
        end
        
    end
    A(sets) = A_won;
    B(sets)  = B_won;
    
end

AverageA = A(sets_max)/(sets_max*iter_max);
AverageB = B(sets_max)/(sets_max*iter_max);

fprintf(fid, 'AverageA: %g\n', AverageA);
fprintf(fid, 'AverageB: %g ', AverageB);
fclose(fid);


