      % Input the model name
      % ^^^^^^^^^^^^^^^^^^^^
      inputstring = ['Enter the name of the ',prefix,' model : '];
      b1 = 1; b2 = 1; b3 = 1; b4 = 1;
    while b4        % check for parameter-redundant models
      while b3      % check that any covariates exist and are of the right length
        while b2,   % check for reasonable model name
          while b1, % check for 3 slashes in the model name
            b1 = 0; b2 = 0; b3 = 0; b4 = 0;
            modelname = input(inputstring,'s');
            slash = findstr(modelname,'/');
            if length(slash) ~= 3,
              b1 = 1;
              fprintf('Model name must consist of 4 submodels separated by slashes\n\n')
            end 
          end  % while b1
         
          b1 = 1;

          model1_orig_name = modelname(1:slash(1)-1);
          modela_orig_name = modelname(slash(1)+1:slash(2)-1);
          modellam_orig_name = modelname(slash(2)+1:slash(3)-1);
          modelp_orig_name = modelname(slash(3)+1:length(modelname));
        
        model1name = cnvrtnm(model1_orig_name);
        modelaname = cnvrtnm(modela_orig_name);
        modellamname = cnvrtnm(modellam_orig_name);
        modelpname = cnvrtnm(modelp_orig_name);
        modelname = [model1name,'_',modelaname,'_',modellamname,'_',modelpname];

          
          s = model1name;
          if length(s)==1,
            if s~='C' & s~='T',
              b2 = 1;
            end
          elseif s(1)~='V',
              b2 = 1;
          end
          if b2 == 1,            
            fprintf('The model name for 1st-year survival must be C or T\n')
            fprintf('or be of the form V(covariate pair list)\n\n')
          end
        
          s = modelaname;
          phiatype = 'age';
          if length(s)==1,
            if s~='C' & s~='T' & s~='A',
              b2 = 1;
            end
          elseif s(1)~='V' & s(1)~='A',
              b2 = 1;
          end
          if s == 'T' | s(1) == 'V', 
            phiatype = 'time';
          elseif s == 'A',
            amax = ncols-1;
          elseif s(1) == 'A',
            x = abs(s(2:length(s))) - 48;
            amax = sum(x.*10.^(length(x)-1:-1:0));
            if amax > ncols-1,
              b2 = 1;
            end
          end
          if b2 == 1,            
            fprintf('The model name for adult survival must be C or T or A or\n')
            fprintf('An for some integer n, or be of the form V(covariate pair list)\n\n')
          end
        
          s = modellamname;
          if length(s)==1,
            if s~='C' & s~='T',
              b2 = 1;
            end
          elseif length(s)==2,
            if ~strcmp(s,'T1'),
              b2 = 1;
            end
          elseif s(1)~='V',
              b2 = 1;
          end
          if b2 == 1,            
            fprintf('The model name for recoveries must be C or T or T1\n')
            fprintf('or be of the form V(covariate list)\n\n')
          end

          s = modelpname;
          if length(s)==1,
            if s~='C'
              b2 = 1;
            end
          elseif s(1)~='V',
              b2 = 1;
          end
          if b2 == 1,            
            fprintf('The model name for productivities must be C\n')
            fprintf('or be of the form V(covariate list)\n\n')
          end
        
        end % while b2
    
        b2 = 1;
        
        % Calculate the X matrices for the model
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        % Calculate the X1's
        % ^^^^^^^^^^^^^^^^^
        Xrr1name = ['Xrr1_',model1name];
        eval([Xrr1name,'=[];'])
        Xc1name = ['Xc1_',model1name];
        eval([Xc1name,'=[];'])

        Xrr1_C = ones(nrows,1);
        Xrr1_T = eye(nrows);
        Xc1_C = ones(ncensus,1);
        Xc1_T = eye(ncensus);

        if model1name == 'T',
          if nrows ~= ncensus,
            fprintf('For this model the ring-recovery and census data must contain the same\n')
            fprintf('1st year survival parameters\n')
            fprintf('You may be able to fit this model to a subset of the data\n')
            b3 = 1;
          end
        elseif  model1name(1)=='V',
          if length(findstr(model1_orig_name,'(')) ~= length(findstr(model1_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression\n')
            b3 = 1;
          else
            splits1 = findstr(model1_orig_name,',');
            ncovs1 = length(splits1)+1;
            if rem(ncovs1,2) ~= 0,
              b3 = 1;
              fprintf(['Incomplete covariate pair(s) for 1st year survival. Add covariate(s)\n'])
            else
              splits1 = [2, splits1, length(model1_orig_name)];
              % Calculate ring-recovery X1
              j=0;
              for i=1:2:ncovs1,
                j=j+1;
                s = model1_orig_name(splits1(i)+1:splits1(i+1)-1);
                paren_op = findstr(s,'(');
                if isempty(paren_op),
                  paren_op = length(s) + 1;
                end
                cov_name = s(1:paren_op-1);
                if ~exist(cov_name),
                  b3 = 1;
                  fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
                elseif length(eval(s)) ~= nrows,
                  b3 =1;
                  fprintf('covariate %2.0f for 1st year survival should be of length %2.0f\n',i,nrows)
                else
                  eval([Xrr1name,'(:,j+1) = ',s,';']);
                end
              end
              % Calculate census X1
              j=0;
              for i=2:2:ncovs1,
                j=j+1;
                s = model1_orig_name(splits1(i)+1:splits1(i+1)-1);
                paren_op = findstr(s,'(');
                if isempty(paren_op),
                  paren_op = length(s) + 1;
                end
                cov_name = s(1:paren_op-1);
                if ~exist(cov_name),
                  b3 = 1;
                  fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
                elseif length(eval(s)) ~= ncensus,
                  b3 = 1;
                  fprintf('covariate %2.0f for 1st year survival should be of length %2.0f\n',i,ncensus)
                else
                  eval([Xc1name,'(:,j+1) = ',s,';']);
                end
              end
            end  %if ~covariate pairs
          end
          if b3 == 0, 
            eval([Xc1name,' = ',Xc1name,'-ones(ncensus,1)*mean(',Xrr1name,');']);
            eval([Xc1name,'(:,1) = ones(ncensus,1);']); 
            eval([Xrr1name,' = ',Xrr1name,'-ones(nrows,1)*mean(',Xrr1name,');']);
            eval([Xrr1name,'(:,1) = ones(nrows,1);']); 
          end
        end
        eval(['Xrr1 = ',Xrr1name,';'])
        eval(['Xc1 = ',Xc1name,';'])
  
        % Calculate the Xa's
        % ^^^^^^^^^^^^^^^^^^
        Xrraname = ['Xrra_',modelaname];
        eval([Xrraname,'=[];'])
        Xcaname = ['Xca_',modelaname];
        eval([Xcaname,'=[];'])

        Xrra_C = ones(ncols-1,1);
        Xrra_T = eye(ncols-1);
        Xca_C = ones(ncensus,1);
        Xca_T = eye(ncensus);

        if modelaname(1)=='A',
          eval([Xrraname,'= [eye(amax); zeros(ncols-1-amax,amax-1), ones(ncols-1-amax,1)];'])
          eval([Xcaname,'= [eye(amax); zeros(ncensus-amax,amax-1), ones(ncensus-amax,1)];'])
        elseif  modelaname(1)=='T',
          if ncols ~= ncensus, 
            fprintf('For this model the recovery and census periods must coincide\n')
            fprintf('You may be able to fit this model to a subset of the census data\n')
            b3 = 1;
          end
        elseif  modelaname(1)=='V',
          if length(findstr(modela_orig_name,'(')) ~= length(findstr(modela_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression\n')
            b3 = 1;
          else
            splitsa = findstr(modela_orig_name,',');
            ncovsa = length(splitsa)+1;
            if rem(ncovsa,2) ~= 0,
              b3 = 1;
              fprintf('Incomplete covariate pair(s) for adult survival. Add covariate(s)\n')
            else
              splitsa = [2, splitsa, length(modela_orig_name)];
              % Calculate ring-recovery Xa
              j=0;
              for i=1:2:ncovsa,
                j=j+1;
                s = modela_orig_name(splitsa(i)+1:splitsa(i+1)-1);
                paren_op = findstr(s,'(');
                if isempty(paren_op),
                  paren_op = length(s) + 1;
                end
                cov_name = s(1:paren_op-1);
                if ~exist(cov_name),
                  b3 = 1;
                  fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
                elseif length(eval(s)) ~= ncols-1,
                  b3 =1;
                  fprintf('Covariate %2.0f for adult survival should be of length %2.0f\n',i,ncols-1)
                  if length(eval(s)) == ncols,
                    fprintf('Don''t forget that there is no adult survival parameter from ring-recovery for the first year\n')
                    fprintf(['You can use ',s,'(2:',ncols,') as a covariate\n'])
                  end
                else
                  eval([Xrraname,'(:,j+1) = ',s,';']);
                end
              end
              % Calculate census Xa
              j=0;
              for i=2:2:ncovsa,
                j=j+1;
                s = modela_orig_name(splitsa(i)+1:splitsa(i+1)-1);
                paren_op = findstr(s,'(');
                if isempty(paren_op),
                  paren_op = length(s) + 1;
                end
                cov_name = s(1:paren_op-1);
                if ~exist(cov_name),
                  b3 = 1;
                  fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
                elseif length(eval(s)) ~= ncensus,
                  b3 =1;
                  fprintf('Covariate %2.0f for adult survival should be of length %2.0f\n',i,ncensus)
                else
                  eval([Xcaname,'(:,j+1) = ',s,';']);
                end
              end
            end %if ~covariate pair
          end
          if b3 == 0,
            eval([Xcaname,' = ',Xcaname,'-ones(ncensus,1)*mean(',Xrraname,');']);
            eval([Xcaname,'(:,1) = ones(ncensus,1);']); 
            eval([Xrraname,' = ',Xrraname,'-ones(ncols-1,1)*mean(',Xrraname,');']);
            eval([Xrraname,'(:,1) = ones(ncols-1,1);']); 
          end
        end
        eval(['Xrra = ',Xrraname,';'])
        eval(['Xca = ',Xcaname,';'])
        
        % Calculate Xlam
        % ^^^^^^^^^^^^^^
        Xlamname = ['Xlam_',modellamname];
        eval([Xlamname,'=[];'])

        Xlam_C = ones(ncols,1);
        Xlam_T = eye(ncols);
        Xlam_T1 = [eye(ncols-1);ones(1,ncols-1)/(ncols-1)];
        
        if  modellamname(1)=='V',
          if length(findstr(modellam_orig_name,'(')) ~= length(findstr(modellam_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression \n')
            b3 = 1;
          else
            splitslam = findstr(modellam_orig_name,',');
            ncovslam = length(splitslam)+1;
            splitslam = [2, splitslam, length(modellam_orig_name)];        
            eval([Xlamname,'(:,1) = ones(ncols,1);'])
            for i=1:ncovslam,
              s = modellam_orig_name(splitslam(i)+1:splitslam(i+1)-1);
              paren_op = findstr(s,'(');
              if isempty(paren_op),
                paren_op = length(s) + 1;
              end
              cov_name = s(1:paren_op-1);
              if ~exist(cov_name),
                b3 = 1;
                fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
              elseif length(eval(s)) ~= ncols,
                b3 =1;
                fprintf('covariate %2.0f for recovery probabilities should be of length %2.0f\n',i,ncols)
              else
                eval([Xlamname,'(:,i+1) = ',s,'-mean(',s,');']);
              end
            end
          end
        end
        eval(['Xlam = ',Xlamname,';'])


        % Calculate Xp
        % ^^^^^^^^^^^^^^
        Xpname = ['Xp_',modelpname];
        eval([Xpname,'=[];'])

        Xp_C = ones(ncensus,1);
        %Xp_T = eye(ncols);
        %Xp_T1 = [eye(ncols)-1);ones(1,ncols-1)/(ncols-1)];
        
        if  modelpname(1)=='V',
          if length(findstr(modelp_orig_name,'(')) ~= length(findstr(modelp_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression \n')
            b3 = 1;
          else
            splitsp = findstr(modelp_orig_name,',');
            ncovsp = length(splitsp)+1;
            splitsp = [2, splitsp, length(modelp_orig_name)];        
            eval([Xpname,'(:,1) = ones(ncensus,1);'])
            for i=1:ncovsp,
              s = modelp_orig_name(splitsp(i)+1:splitsp(i+1)-1);
              paren_op = findstr(s,'(');
              if isempty(paren_op),
                paren_op = length(s) + 1;
              end
              cov_name = s(1:paren_op-1);
              if ~exist(cov_name),
                b3 = 1;
                fprintf(['There is no such covariate as ',s,' Check spelling.\n'])
              elseif length(eval(s)) ~= ncensus,
                b3 =1;
                fprintf('covariate %2.0f for productivities should be of length %2.0f\n',i,ncensus)
              else
                eval([Xpname,'(:,i+1) = ',s,'-mean(',s,');']);
              end
            end
          end
        end
        eval(['Xp = ',Xpname,';'])

      end    % while b3

      b3 = 1;    

      if strcmp(modelname,'C_A_C') | (strcmp(modelaname,'T') & strcmp(modellamname,'T')),
        fprintf(['Sorry, model ',modelname,' is not identifiable\n\n'])
        b4 = 1;
      end
    end  % while b4

    b4 = 1;
