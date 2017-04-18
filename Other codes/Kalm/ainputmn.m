      % Input the model name
      % ^^^^^^^^^^^^^^^^^^^^
      inputstring = ['Enter the name of the ',prefix,' model : '];
      b1 = 1; b2 = 1; b3 = 1; b4 = 1;
    while b4        % check that model has "nmu" survival parameters
      while b3      % check that any covariates exist and are of the right length
        while b2,   % check for reasonable model name
          while b1, % check for 2 slashes in the model name
            b1 = 0; b2 = 0; b3 = 0; b4 = 0;
            modelname = input(inputstring,'s');
            slash = findstr(modelname,'/');
            if length(slash) ~= 2,
              b1 = 1;
              fprintf('Model name must consist of 3 submodels separated by slashes\n\n')
            end 
          end  % while b1
         
          b1 = 1;

          model1_orig_name = modelname(1:slash(1)-1);
          modela_orig_name = modelname(slash(1)+1:slash(2)-1);
          modelp_orig_name = modelname(slash(2)+1:length(modelname));
        
        model1name = cnvrtnm(model1_orig_name);
        modelaname = cnvrtnm(modela_orig_name);
        modelpname = cnvrtnm(modelp_orig_name);
        modelname = [model1name,'_',modelaname,'_',modelpname];

          
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
            if s~='C' & s~='T',
              b2 = 1;
            end
          elseif s(1)~='V' & s(1)~='A',
              b2 = 1;
          end
          if s == 'T' | s(1) == 'V', 
            phiatype = 'time';
          elseif s(1) == 'A',
            x = abs(s(2:length(s))) - 48;
            amax = sum(x.*10.^(length(x)-1:-1:0));
          end
          if b2 == 1,            
            fprintf('The model name for adult survival must be C or T or\n')
            fprintf('An for some integer n, or be of the form V(covariate pair list)\n\n')
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
        X1name = ['X1_',model1name];
        eval([X1name,'=[];'])

        X1_C = ones(ncensus,1);
        X1_T = eye(ncensus);

        if  model1name(1)=='V',
          if length(findstr(model1_orig_name,'(')) ~= length(findstr(model1_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression\n')
            b3 = 1;
          else
            splits1 = findstr(model1_orig_name,',');
            ncovs1 = length(splits1)+1;
            splits1 = [2, splits1, length(model1_orig_name)];
            eval([X1name,'(:,1) = ones(ncensus,1);'])
            for i=1:ncovs1,
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
                b3 =1;
                fprintf('covariate %2.0f for 1st year survival should be of length %2.0f\n',i,ncensus)
              else
                  %eval([X1name,'(:,i+1) = ',s,'-mean(',s,');']);
                  eval([X1name,'(:,i+1) = ',s,';']);
              end
            end
          end
        end
        eval(['X1 = ',X1name,';'])
  
        % Calculate the Xa's
        % ^^^^^^^^^^^^^^^^^^
        Xaname = ['Xa_',modelaname];
        eval([Xaname,'=[];'])

        Xa_C = ones(ncensus,1);
        Xa_T = eye(ncensus);

        if modelaname(1)=='A',
          eval([Xaname,'= [eye(amax); zeros(ncensus-amax,amax-1), ones(ncensus-amax,1)];'])
        elseif  modelaname(1)=='V',
          if length(findstr(modela_orig_name,'(')) ~= length(findstr(modela_orig_name,')')),
            fprintf('Unmatched parentheses in covariate expression\n')
            b3 = 1;
          else
            splitsa = findstr(modela_orig_name,',');
            ncovsa = length(splitsa)+1;
            splitsa = [2, splitsa, length(modela_orig_name)];
            eval([Xaname,'(:,1) = ones(ncensus,1);'])
            for i=1:ncovsa,
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
                %eval([Xaname,'(:,i+1) = ',s,'-mean(',s,');']);
                eval([Xaname,'(:,i+1) = ',s,';']);
              end
            end
          end
        end
        eval(['Xa = ',Xaname,';'])
        

        % Calculate Xp
        % ^^^^^^^^^^^^^^
        Xpname = ['Xp_',modelpname];
        eval([Xpname,'=[];'])

        Xp_C = ones(ncensus,1);
        
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

      npar1 = size(X1,2);
      npara = size(Xa,2); 
      if npar1+npara-nmu ~=0, 
        fprintf('Model should have %2.0f survival parameters\n',nmu)
        b4 = 1;
      end
    end  % while b4

    b4 = 1;
