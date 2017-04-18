% Copyright (C) 2003 by P.T. Besbeas
 % All rights reserved.
 %
 % This software may be freely copied, modified and redistributed without
 % fee for non-commerical purposes provided that this copyright notice is
 % preserved intact on all copies and modified copies.
 % 
 % There is no warranty or other guarantee of fitness of this software.
 % It is provided solely "as is". The author disclaims all
 % responsibility and liability with respect to this software's usage
 % or its effect upon hardware or computer systems.
 
% required files:
% calcsv.m
% chisq.m
% cnvrtnm.m
% grad.m
% ilogit.m
% infoT.m
% inputmn.m
% lfact.m
% likT.m
% lt.m
% mexpT.m
% pchisq.m
% preproc.m
% probT.m
% prtupsgn.m
% prtupmat.m
% scoreJT.m
% stretch.m

global Ringed Mobs unrec nrows ncols Census ncensus const Toe ellmax sopts uopts Mu Sigma nmu;

% more on
matlab_varlength = 19;   % the maximum length of a variable name in Matlab
maxerr = 0.001;

fprintf('          -----------------------------------------------------------\n')
fprintf('         | Kalm version 1.0 . . by Besbeas & Morgan . . 2003/07/01   |\n')
fprintf('         |                                                           |\n')
fprintf('         |    A program for fitting ring-recovery and census data    |\n')
fprintf('         |                                                           |\n')
fprintf('         |    Bug reports, praise, offers of large grants etc, to    |\n')
fprintf('         |    p.t.besbeas@kent.ac.uk      b.j.t.morgan@kent.ac.uk    |\n')
fprintf('          ----------------------------------------------------------- \n\n')

% Introduction
% ^^^^^^^^^^^^
type = input('Is this (1) an exact analysis (2) an approximate analysis? [1] :','s');
if ( isempty(type) | type ~= '2' ), 
  type = 1;
  newsession = input('Is this (1) a new session (2) a continuation session? [1] : ','s');
  if ( isempty(newsession) | newsession ~= '2' ),
    newsession = 1;
    fprintf('\n\n')
    fprintf('The data should be stored in a matlab file in the form of a vector\n')
    fprintf('"Ringed" of ringing numbers, a matrix "Mobs" of recovery numbers \n')
    fprintf('and a vector "Census" of census numbers. \n')
    fprintf('For example the matlab file example.m might contain the lines\n\n')
    fprintf('Ringed = [1000;1000;1000;1000];\n\n')
    fprintf('Mobs = [35  7  0  1\n')
    fprintf('         0 23 11  3\n')
    fprintf('         0  0 45  13\n')
    fprintf('         0  0  0  33 ];\n\n')
    fprintf('Census = [500;700;800;600;900];\n\n')
  
    % Read the data
    % ^^^^^^^^^^^^^
    data_file_exists = 0;
    while  data_file_exists ~= 2,
      datafilename = input('Enter the name of the data file, without the ".m" extension : ','s');
      data_file_exists = exist([datafilename,'.m']);
  
      if data_file_exists ~= 2,
        fprintf(['I could not find file ',datafilename,'.m']);
      end
    end;  % while data_file_exists
    fprintf('\n')
  
    diarydefaultname = [datafilename,'.res'];
    eval(datafilename);
    preproc;

  else % if newsession=='2'

    % Load the old results
    % ^^^^^^^^^^^^^^^^^^^^
    load_file_exists = 0;
    while  load_file_exists ~= 2,
      fprintf('\nEnter name of file of stored results, without the ".mat" extension\n')
      loadfilename = input(['File name = : '],'s');
      load_file_exists = exist([loadfilename,'.mat']);
      if load_file_exists ~= 2,
        fprintf(['I could not find file ',loadfilename,'.m\n']);
      end
    end;  % while load_file_exists
    fprintf('\n')

    eval(['load ',loadfilename])
    fprintf(['Previous results and data from ',loadfilename,'.mat now loaded\n'])

    diarydefaultname = [loadfilename,'.res'];

  end  % if newsession

a = input(['Enter name of file to store session transcript [',diarydefaultname,'] : '],'s');
if isempty(a),
  diaryname = diarydefaultname;
else
  diaryname = a;
end
diary(diaryname);
diarybox = '--';
for i = 1:length(diaryname),
  diarybox = ['-',diarybox];
end
fprintf([' ',diarybox,'\n'])
fprintf(['| ',diaryname,' |\n'])
fprintf([' ',diarybox,'\n\n'])

if newsession==1,
  fprintf('\nModel names use the notation x/y/z/w\n')
  fprintf('where x = C      if first-year survival probabilities are constant\n')
  fprintf('          T      if they vary with time (year of ringing)\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('      y = C      if adult survival probabilities are constant\n')
  fprintf('          A      if they are fully age dependent\n')
  fprintf('          An     if they vary with age up to age n\n')
  fprintf('          T      if they vary with time (year of recovery)\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('      z = C      if recovery probabilities are constant\n')
  fprintf('          T      if they vary with time (year of recovery)\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('      w = C      if productivity rates are constant\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('The covariate names should be specified inside the parentheses,\nseparated by commas.\n')
  fprintf('The survival probabilities should depend on pairs of variables per covariate,\n')
  fprintf('one variable for the ring-recovery years and the other for the census years.\n')
  fprintf('For example, if the 1st-year survival probabilities depend on covariates\n')
  fprintf('"January" and "February", the adult survival probabilities vary with age\n')
  fprintf('up to age 1 (i.e. have the same value in the 2nd, 3rd,... years of life),\n')
  fprintf('and the recovery and productivity rates are constant,\n')
  fprintf('then the model should be specified as V(Janrr,Janc,Febrr,Febc)/A2/C/C\n\n')
  fprintf('Do not use covariate names that include spaces.\n\n')
  fprintf('Press any key when ready . . .\n')
  pause
end

% sopts=optimset('TolX',2e-10,'TolFun',2e-10,'MaxFunEval',5000000, 'MaxIter',50000000);
% uopts=optimset('TolX',2e-4,'TolFun',2e-4,'MaxFunEval',5000000, 'MaxIter',50000000);
uopts=optimset('LargeScale','off','TolFun',2e-6, 'MaxFunEvals',1000000000000, 'MaxIter',1000000000);

while 1,
    fprintf('--------------------------------------------------------------------\n')
    fprintf('\nDo you want to (1) fit a model\n')
    fprintf('or             (0) quit from this program ?\n\n')
    dofit = input('Please type 0 or 1 : ');
    
    if dofit == 0,
      diary off;
      a = input('Do you want to store the results of this session? [Y] : ','s');
      if isempty(a),
        b = 1;
      elseif a(1)=='y' | a(1)=='Y',
        b = 1;
      else b = 0;
      end
      if b, 
          if newsession==1,
            matfilename = datafilename;
          else
            matfilename = loadfilename;
          end
          fprintf('Enter name of file to store results in, without the ".mat" extension\n')
          a = input(['File name = [',matfilename,'] : '],'s');
          if isempty(a),
            savefilename = matfilename;
          else
            savefilename = a;
          end
          clear loadfilename datafilename newsession
          eval(['save ',savefilename])
          fprintf(['Results saved in ',savefilename,'.mat\n'])
        end
      fprintf('\nThank you for using Kalm.\n')
      return

    elseif dofit == 1,
      % Input the model name
      % ^^^^^^^^^^^^^^^^^^^^
      prefix = '';
      inputmn;
      
      % Calculating starting values
      % ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      calcsv;

      % Fit the model
      % ^^^^^^^^^^^^^
      fprintf('Estimating the parameters, please wait . . .\n\n')
      betaname = ['b_',modelname];
      pname = ['p_',modelname];
      temptime = clock;
%       if exist('fminunc')==2,
        eval([betaname,' = fminunc(''likj'',betastart,uopts,',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
%       else
%         eval([betaname,' = fminsearch(''likj'',betastart,sopts,',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
%       end
      runtime = etime(clock,temptime);
      fprintf('time taken =%8.2f secs\n\n',runtime)
      fprintf('Checking for convergence\n')
      fprintf('and calculating standard errors, please wait . . .\n\n')
      temptime = clock;
      eval([pname,' = [ilogit(',betaname,'(1:npar1+npara+nparlam));','exp(',betaname,'(npar1+npara+nparlam+1:npars))];'])
      eval(['g = grad(''likj'',',betaname,',Xrr1,Xc1,Xrra,Xca,Xlam,Xp,phiatype);'])
      eval(['jay = hessian(''likj'',',betaname,',',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
      runtime = etime(clock,temptime);
      fprintf('time taken =%8.2f secs\n\n',runtime)
      fprintf('gradient = %8.4g\n\n',norm(g))
      errb = jay\g;
      eval(['beta = ',betaname,';'])
      errp = ilogit(beta) - ilogit(beta-errb);

      newsopts = sopts;
      newuopts = uopts;
      while max(errp) > maxerr,
        fprintf('This gradient is too large. The likelihood maximisation has not converged\n')
        fprintf('fully. This introduces an extra error into the parameter estimates,\n')
        fprintf('of the order of %6.4f\n\n',max(errp))
        another_run = input('Would you like a new run with a smaller tolerance? [No] ','s');
        if isempty(another_run),
          break
        elseif another_run(1)=='N'| another_run(1)=='n',
          break
        else 
          newsTol = optimget(newsopts,'TolX')/sqrt(10);
          newsopts = optimset(newsopts,'TolX',newsTol,'TolFun',newsTol);
          newuTol = optimget(newuopts,'TolX')/sqrt(10);
          newuopts = optimset(newuopts,'TolX',newuTol,'TolFun',newuTol);
%           if exist('fminunc')==2,
            eval([betaname,' = fminunc(''likj'',betastart,newuopts,',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
%           else
%             eval([betaname,' = fminsearch(''likj'',betastart,newsopts,',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
%           end
          runtime = etime(clock,temptime);
          fprintf('time taken =%8.2f secs\n\n',runtime)
          fprintf('Checking for convergence\n')
          fprintf('and calculating standard errors, please wait . . .\n\n')
          temptime = clock;
          eval([pname,' = [ilogit(',betaname,'(1:npar1+npara+nparlam));','exp(',betaname,'(npar1+npara+nparlam+1:npars))];'])
          eval(['g = grad(''likj'',',betaname,',Xrr1,Xc1,Xrra,Xca,Xlam,Xp,phiatype);'])
          eval(['jay = hessian(''likj'',',betaname,',',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
          runtime = etime(clock,temptime);
          fprintf('time taken =%8.2f secs\n\n',runtime)
          fprintf('gradient = %8.4g\n\n',norm(g))
          errb = jay\g;
          eval(['beta = ',betaname,';'])
          errp = ilogit(beta) - ilogit(beta-errb);
        end % if another_run
      end % while max(errp) > maxerr

      seb = sqrt(diag(inv(jay)));
      eval(['sep = (ilogit(',betaname,'+0.5*seb)-ilogit(',betaname,'-0.5*seb));'])
      indxparc = npars-nparp:npars;
      eval(['sep(indxparc) = (exp(',betaname,'(indxparc)).*seb(indxparc));'])

      % Correlation matrix
      corname = ['cor_',modelname];
      V = inv(jay);
      D = diag(seb.^(-1));
      C = D*V*D;
      eval([corname,' = C;']);

%     Output the results
%     ^^^^^^^^^^^^^^^^^^
      if model1name(1)=='V',

        fprintf('The logistic regression parameter estimates for 1st-year survival are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(1)']),seb(1))
        j=0;
        for i=1:2:ncovs1,
          j=j+1;
          s = model1_orig_name(splits1(i)+1:splits1(i+2)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(j+1)']),seb(j+1))
        end
        fprintf('\n')
        eval(['beta1 = ',betaname,'(1:npar1);'])
        V1 = V(1:npar1,1:npar1);
        yhat = Xrr1*beta1;
        phi1hat = ilogit(yhat);
        sey = sqrt(diag(Xrr1*V1*Xrr1'));
        sephi1 = ilogit(yhat + 0.5*sey) - ilogit(yhat - 0.5*sey);
        eval(['phi1_',modelname,'=phi1hat;'])
        fprintf('The corresponding 1st-year survival probabilities will be saved')
        a = input(['as phi1_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:nrows,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phi1hat(i),sephi1(i))
           end
        end
      else
        fprintf('The estimated 1st-year survival probabilities are\n')
        for i = 1:npar1,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      if modelaname(1)=='V',
        fprintf('The logistic regression parameter estimates for adult survival are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(npar1+1)']),seb(npar1+1))
        j=0;
        for i=1:2:ncovsa,
          j=j+1;
          s = modela_orig_name(splitsa(i)+1:splitsa(i+2)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(npar1+j+1)']),seb(npar1+j+1))
        end
        fprintf('\n')
        eval(['betaa = ',betaname,'(npar1+1:npar1+npara);'])
        Va = V(npar1+1:npar1+npara,npar1+1:npar1+npara);
        yhat = Xrra*betaa;
        phiahat = ilogit(yhat);
        sey = sqrt(diag(Xrra*Va*Xrra'));
        sephia = ilogit(yhat + 0.5*sey) - ilogit(yhat - 0.5*sey);
        eval(['phia_',modelname,'=phiahat;'])
        fprintf('The corresponding adult survival probabilities will be saved')
        a = input(['as phia_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncols-1,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phiahat(i),sephia(i))
           end
        end
      else
        fprintf('The estimated adult survival probabilities are\n')
        if modelaname(1) == 'T', 
          fprintf('(Note the 1st probability below is estimated from census data alone)\n')
        end
        for i = npar1+1:npar1+npara,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      if modellamname(1)=='V',
        fprintf('The logistic regression parameter estimates for recovery rates are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(npar1+npara+1)']),seb(npar1+npara+1))
        for i=1:ncovslam,
          s = modellam_orig_name(splitslam(i)+1:splitslam(i+1)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(npar1+npara+i+1)']),seb(npar1+npara+i+1))
        end
        fprintf('\n')
        eval(['betalam = ',betaname,'(npar1+npara+1:npar1+npara+nparlam);'])
        Vlam = V(npar1+npara+1:npar1+npara+nparlam,npar1+npara+1:npar1+npara+nparlam);
        yhat = Xlam*betalam;
        lamhat = ilogit(yhat);
        sey = sqrt(diag(Xlam*Vlam*Xlam'));
        selam = ilogit(yhat + 0.5*sey) - ilogit(yhat - 0.5*sey);
        eval(['lam_',modelname,'=lamhat;'])
        fprintf('The corresponding recovery probabilities will be saved')
        a = input(['as lam_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncols,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,lamhat(i),selam(i))
           end
        end
      else
        fprintf('The estimated recovery probabilities are\n')
        for i = npar1+npara+1:npar1+npara+nparlam,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      if modelpname(1)=='V',
        fprintf('The logarithmic regression parameter estimates for productivity rates are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(npar1+npara+nparlam+1)']),seb(npar1+npara+nparlam+1))
        for i=1:ncovsp,
          s = modelp_orig_name(splitsp(i)+1:splitsp(i+1)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(npar1+npara+nparlam+i+1)']),seb(npar1+npara+nparlam+i+1))
        end
        fprintf('\n')
        eval(['betap = ',betaname,'(npar1+npara+nparlam+1:npars-1);'])
        Vp = V(npar1+npara+nparlam+1:npars-1,npar1+npara+nparlam+1:npars-1);
        yhat = Xp*betap;
        phat = exp(yhat);
        sey = sqrt(diag(Xp*Vp*Xp'));
        seprod = exp(yhat).*sey;
        eval(['prod_',modelname,'=phat;'])
        fprintf('The corresponding productivity rates will be saved')
        a = input(['as prod_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncensus,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phat(i),seprod(i))
           end
        end
      else
        fprintf('The estimated productivity rates are\n')
        for i = npar1+npara+nparlam+1:npars-1,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      fprintf('The estimated measurement error st. deviation is\n')
      fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(npars)']),sep(npars));
      fprintf('\n')

%     Print correlation matrix
%     ^^^^^^^^^^^^^^^^^^^^^^^^
      fprintf('\n')
      fprintf(['The estimated parameter correlations are stored in the matrix ',corname,'\n'])
      a = input('Do you want to view them now [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y'
          fprintf(['The parameters appear in ',corname,' in the order\n'])
          fprintf('first-year survival, adult survival, recovery probabilities,\n')
          fprintf('productivity rate, measurement error st. deviation\n')
          prtupmat(C,2);
        end
      end
      fprintf('\n')

%     Print max loglik and AIC
%     ^^^^^^^^^^^^^^^^^^^^^^^^
      likname = ['lik_',modelname];
      eval([likname,' = -likj(',betaname,',',Xrr1name,',',Xc1name,',',Xrraname,',',Xcaname,',',Xlamname,',',Xpname,',''',phiatype,''');'])
      fprintf('The maximized loglikelihood is %-10.4f\n',eval(likname))
        
      AIC = 2*(npars-eval(likname));
      fprintf('The Akaike information criterion is %-10.3f\n\n',AIC)
    
      eval(['Mexp = mexpT(',betaname,'(1:npar1+npara+nparlam),',Xrr1name,',',Xrraname,',',Xlamname,',''',phiatype,''');'])
      X2 = chisq(Mobs,Mexp);
      df = nrows*(2*ncols-nrows+1)/2 - npars;
      pval = pchisq(X2,df);
      fprintf('The naive Pearson chi-square goodness-of-fit statistic is \n')
      fprintf('           X2 = %-8.4g on %3d degrees of freedom;  p-value = %6.3g\n',X2,df,pval)
      %dev = 2*(lik_max - eval(likname));
      %pval = pchisq(dev,df);
      %fprintf('The deviance is %-8.4g on %3d degrees of freedom;  p-value = %6.3g\n\n',dev,df,pval)

%     Print fitted values and residual pattern
%     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      a = input('Do you want to see the matrix of fitted recoveries [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y',
          prtupmat(Mexp,2);
        end
      end

      resid = Mobs - Mexp;
      signres = (resid>1) - (resid<-1);
      a = input('Do you want to see the pattern of residuals [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y',
          fprintf(' + denotes a residual >  1,\n')
          fprintf(' - denotes a residual < -1,\n')
          fprintf(' . denotes a residual with absolute value < 1\n\n')
          prtupsgn(signres);
        end
      end
      fprintf('\n')
     
%     Print smoothed estimates and fitted values
%     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      a0=zeros(2,1); P0=1e6*eye(2);
      smoothname = ['smooth_',modelname];
      eval([smoothname,' = ksc(',betaname,'([1:npar1+npara npars-nparp:npars]),',Xc1name,',',Xcaname,',',Xpname,',''',phiatype,''');'])
      fprintf(['The smoothed estimates of the population numbers are stored in the matrix ',smoothname,'\n'])
      a = input('Do you want to view them now [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y'
          fprintf(['The estimates appear in ',smoothname,', in columns, in the order\n'])
          if modelaname(1) == 'C' | modelaname(1) == 'T',
            fprintf('[1-year olds | adults]\n')
          else
            fprintf('[1-year olds | 2-year olds | ... | adults]\n')
          end
          eval(['prtupmat(',smoothname,',1);'])
        end
      end
      fprintf('\n')

      a = input('Do you want to plot the observed vs the fitted census [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y',
          plot(Census,':')
          hold on
          eval(['plot(sum(',smoothname,'(:,2:end),2));'])
          legend('observed','fitted')
          hold off
        end
      end
      fprintf('\n')

  end %  if dofit == 1
  
end % endless while loop


else % if type == 2
  newsession = input('Is this (1) a new session (2) a continuation session? [1] : ','s');
  if ( isempty(newsession) | newsession ~= '2' ),
    newsession = 1;
    fprintf('\n\n')
    fprintf('The data should be stored in a matlab file in the form of a vector\n "Mu"')
    fprintf('of the ring-recovery parameter estimates for survival, in logistic space, a matrix\n')
    fprintf('"Sigma" of their variances-covariances and a vector "Census" of census numbers.\n')
    fprintf('For example the matlab file aexample.m might contain the lines\n\n')
    fprintf('Mu = [-0.982;-1.364];\n\n')
    fprintf('Sigma = [0.038  0.028\n')
    fprintf('         0.028  0.309];\n\n')
    fprintf('Census = [500;700;800;600;900];\n\n')
  
    % Read the data
    % ^^^^^^^^^^^^^
    data_file_exists = 0;
    while  data_file_exists ~= 2,
      datafilename = input('Enter the name of the data file, without the ".m" extension : ','s');
      data_file_exists = exist([datafilename,'.m']);

      if data_file_exists ~= 2,
        fprintf(['I could not find file ',datafilename,'.m\n\n']);
      end
    end;  % while data_file_exists
    fprintf('\n')
  
    diarydefaultname = [datafilename,'.res'];
    eval(datafilename);
    apreproc; 

  else % if newsession=='2'

    % Load the old results
    % ^^^^^^^^^^^^^^^^^^^^
    load_file_exists = 0;
    while  load_file_exists ~= 2,
      fprintf('\nEnter name of file of stored results, without the ".mat" extension\n')
      loadfilename = input(['File name = : '],'s');
      load_file_exists = exist([loadfilename,'.mat']);
      if load_file_exists ~= 2,
        fprintf(['I could not find file ',loadfilename,'.m\n']);
      end
    end;  % while load_file_exists
    fprintf('\n')

    eval(['load ',loadfilename])
    fprintf(['Previous results and data from ',loadfilename,'.mat now loaded\n'])

    diarydefaultname = [loadfilename,'.res'];

  end  % if newsession

a = input(['Enter name of file to store session transcript [',diarydefaultname,'] : '],'s');
if isempty(a),
  diaryname = diarydefaultname;
else
  diaryname = a;
end
diary(diaryname);
diarybox = '--';
for i = 1:length(diaryname),
  diarybox = ['-',diarybox];
end
fprintf([' ',diarybox,'\n'])
fprintf(['| ',diaryname,' |\n'])
fprintf([' ',diarybox,'\n\n'])

if newsession==1,
  fprintf('\nModel names use the notation x/y/w\n')
  fprintf('where x = C      if first-year survival probabilities are constant\n')
  fprintf('          T      if they vary with time (year of ringing)\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('      y = C      if adult survival probabilities are constant\n')
  fprintf('          An     if they vary with age up to age n\n')
  fprintf('          T      if they vary with time (year of recovery)\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('      w = C      if productivity rates are constant\n')
  fprintf('          V(...) if they depend on one or more time-dependent covariates\n\n')
  fprintf('The covariate names should be specified inside the parentheses,\nseparated by commas.\n')
  fprintf('For example, if the 1st-year survival probabilities depend on covariates\n')
  fprintf('"January" and "February", the adult survival probabilities vary with age\n')
  fprintf('up to age 3 (i.e. have the same value in the 4th, 5th,... years of life),\n')
  fprintf('and the productivity rate is constant,\n')
  fprintf('then the model should be specified as V(January,February)/A3/C\n\n')
  fprintf('Do not use covariate names that include spaces and centre and scale the\n')
  fprintf('covariates for survival as in the ring-recovery exercise.\n\n')
  fprintf('Press any key when ready . . .\n')
  pause
end

sopts=optimset('TolX',2e-10,'TolFun',2e-10,'MaxFunEval',5000000, 'MaxIter',50000000);
uopts=optimset('TolX',2e-4,'TolFun',2e-4,'MaxFunEval',5000000, 'MaxIter',50000000);
%uopts=optimset('LargeScale','off','TolX',2e-4,'TolFun',2e-4,'MaxFunEval',50000);

while 1,
    fprintf('--------------------------------------------------------------------\n')
    fprintf('\nDo you want to (1) fit a model\n')
    fprintf('or             (0) quit from this program ?\n\n')
    dofit = input('Please type 0 or 1 : ');
    
    if dofit == 0,
      diary off;
      a = input('Do you want to store the results of this session? [Y] : ','s');
      if isempty(a),
        b = 1;
      elseif a(1)=='y' | a(1)=='Y',
        b = 1;
      else b = 0;
      end
      if b, 
          if newsession==1,
            matfilename = datafilename;
          else
            matfilename = loadfilename;
          end
          fprintf('Enter name of file to store results in, without the ".mat" extension\n')
          a = input(['File name = [',matfilename,'] : '],'s');
          if isempty(a),
            savefilename = matfilename;
          else
            savefilename = a;
          end
          clear loadfilename datafilename newsession
          eval(['save ',savefilename])
          fprintf(['Results saved in ',savefilename,'.mat\n'])
        end
      fprintf('\nThank you for using Kalm.\n')
      return

    elseif dofit == 1,
      % Input the model name
      % ^^^^^^^^^^^^^^^^^^^^
      prefix = '';
      ainputmn;
      
      % Calculating starting values
      % ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      acalcsv;

      % Fit the model
      % ^^^^^^^^^^^^^
      fprintf('Estimating the parameters, please wait . . .\n\n')
      betaname = ['b_',modelname];
      pname = ['p_',modelname];
      temptime = clock;
%       if exist('fminunc')==2,
        eval([betaname,' = fminunc(''alikj'',betastart,uopts,',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
%       else
%         eval([betaname,' = fminsearch(''alikj'',betastart,sopts,',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
%       end
      runtime = etime(clock,temptime);
      fprintf('time taken =%8.2f secs\n\n',runtime)
      fprintf('Checking for convergence\n')
      fprintf('and calculating standard errors, please wait . . .\n\n')
      temptime = clock;
      eval([pname,' = [ilogit(',betaname,'(1:npar1+npara));','exp(',betaname,'(npar1+npara+1:npars))];'])
      eval(['g = grad(''alikj'',',betaname,',X1,Xa,Xp,phiatype);'])
      eval(['jay = hessian(''alikj'',',betaname,',',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
      runtime = etime(clock,temptime);
      fprintf('time taken =%8.2f secs\n\n',runtime)
      fprintf('gradient = %8.4g\n\n',norm(g))
      errb = jay\g;
      eval(['beta = ',betaname,';'])
      errp = ilogit(beta) - ilogit(beta-errb);

      newsopts = sopts;
      newuopts = uopts;
      while max(errp) > maxerr,
        fprintf('This gradient is too large. The likelihood maximisation has not converged\n')
        fprintf('fully. This introduces an extra error into the parameter estimates,\n')
        fprintf('of the order of %6.4f\n\n',max(errp))
        another_run = input('Would you like a new run with a smaller tolerance? [No] ','s');
        if isempty(another_run),
          break
        elseif another_run(1)=='N'| another_run(1)=='n',
          break
        else 
          newsTol = optimget(newsopts,'TolX')/sqrt(10);
          newsopts = optimset(newsopts,'TolX',newsTol,'TolFun',newsTol);
          newuTol = optimget(newuopts,'TolX')/sqrt(10);
          newuopts = optimset(newuopts,'TolX',newuTol,'TolFun',newuTol);
%           if exist('fminunc')==2,
            eval([betaname,' = fminunc(''alikj'',betastart,newuopts,',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
%           else
%             eval([betaname,' = fminsearch(''alikj'',betastart,newsopts,',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
%           end
          runtime = etime(clock,temptime);
          fprintf('time taken =%8.2f secs\n\n',runtime)
          fprintf('Checking for convergence\n')
          fprintf('and calculating standard errors, please wait . . .\n\n')
          temptime = clock;
          eval([pname,' = [ilogit(',betaname,'(1:npar1+npara));','exp(',betaname,'(npar1+npara+1:npars))];'])
          eval(['g = grad(''alikj'',',betaname,',X1,Xa,Xp,phiatype);'])
          eval(['jay = hessian(''alikj'',',betaname,',',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
          runtime = etime(clock,temptime);
          fprintf('time taken =%8.2f secs\n\n',runtime)
          fprintf('gradient = %8.4g\n\n',norm(g))
          errb = jay\g;
          eval(['beta = ',betaname,';'])
          errp = ilogit(beta) - ilogit(beta-errb);
        end % if another_run
      end % while max(errp) > maxerr

      seb = sqrt(diag(inv(jay)));
      eval(['sep = (ilogit(',betaname,'+0.5*seb)-ilogit(',betaname,'-0.5*seb));'])
      indxparc = npars-nparp:npars;
      eval(['sep(indxparc) = (exp(',betaname,'(indxparc)).*seb(indxparc));'])

      % Correlation matrix
      corname = ['cor_',modelname];
      V = inv(jay);
      D = diag(seb.^(-1));
      C = D*V*D;
      eval([corname,' = C;']);

%     Output the results
%     ^^^^^^^^^^^^^^^^^^
      if model1name(1)=='V',

        fprintf('The logistic regression parameter estimates for 1st-year survival are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(1)']),seb(1))
        for i=1:ncovs1,
          s = model1_orig_name(splits1(i)+1:splits1(i+1)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(i+1)']),seb(i+1))
        end
        fprintf('\n')
        eval(['beta1 = ',betaname,'(1:npar1);'])
        V1 = V(1:npar1,1:npar1);
        yhat = X1*beta1;
        phi1hat = ilogit(yhat);
        sey = sqrt(diag(X1*V1*X1'));
        sephi1 = ilogit(yhat + 0.5*sey) - ilogit(yhat - 0.5*sey);
        eval(['phi1_',modelname,'=phi1hat;'])
        fprintf('The corresponding 1st-year survival probabilities will be saved')
        a = input(['as phi1_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncensus,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phi1hat(i),sephi1(i))
           end
        end
      else
        fprintf('The estimated 1st-year survival probabilities are\n')
        for i = 1:npar1,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      if modelaname(1)=='V',
        fprintf('The logistic regression parameter estimates for adult survival are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(npar1+1)']),seb(npar1+1))
        for i=1:ncovsa,
          j=j+1;
          s = modela_orig_name(splitsa(i)+1:splitsa(i+1)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(npar1+i+1)']),seb(npar1+i+1))
        end
        fprintf('\n')
        eval(['betaa = ',betaname,'(npar1+1:npar1+npara);'])
        Va = V(npar1+1:npar1+npara,npar1+1:npar1+npara);
        yhat = Xa*betaa;
        phiahat = ilogit(yhat);
        sey = sqrt(diag(Xa*Va*Xa'));
        sephia = ilogit(yhat + 0.5*sey) - ilogit(yhat - 0.5*sey);
        eval(['phia_',modelname,'=phiahat;'])
        fprintf('The corresponding adult survival probabilities will be saved')
        a = input(['as phia_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncensus,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phiahat(i),sephia(i))
           end
        end
      else
        fprintf('The estimated adult survival probabilities are\n')
        for i = npar1+1:npar1+npara,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      if modelpname(1)=='V',
        fprintf('The logarithmic regression parameter estimates for productivity rates are\n')
        fprintf('Constant        %8.4f\t\tst.err = %7.4f\n', eval([betaname,'(npar1+npara+1)']),seb(npar1+npara+1))
        for i=1:ncovsp,
          s = modelp_orig_name(splitsp(i)+1:splitsp(i+1)-1);
          blanks = [];
          for blanksi = 1:16-length(s)
            blanks = [blanks,' '];
          end
          fprintf([s,blanks,'%8.4f\t\tst.err = %7.4f\n'],eval([betaname,'(npar1+npara+i+1)']),seb(npar1+npara+i+1))
        end
        fprintf('\n')
        eval(['betap = ',betaname,'(npar1+npara+1:npars-1);'])
        Vp = V(npar1+npara+1:npars-1,npar1+npara+1:npars-1);
        yhat = Xp*betap;
        phat = exp(yhat);
        sey = sqrt(diag(Xp*Vp*Xp'));
        seprod = exp(yhat).*sey;
        eval(['prod_',modelname,'=phat;'])
        fprintf('The corresponding productivity rates will be saved')
        a = input(['as prod_',modelname,'. Want to see them now? [N] :'],'s');
        if isempty(a),
          a = 'N';
        elseif a(1)=='y' | a(1)=='Y',
          a = 'Y';;
        end
        if a == 'Y',
           for i = 1:ncensus,
             fprintf('%2.0f:\t\t%8.4f\t\tst.err = %7.4f\n',i,phat(i),seprod(i))
           end
        end
      else
        fprintf('The estimated productivity rates are\n')
        for i = npar1+npara+1:npars-1,
          fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(i)']),sep(i));
        end
      end
      fprintf('\n')

      fprintf('The estimated measurement error st. deviation is\n')
      fprintf('\t\t%8.4f\t\tst.err = %7.4f\n',eval([pname,'(npars)']),sep(npars));
      fprintf('\n')

%     Print correlation matrix
%     ^^^^^^^^^^^^^^^^^^^^^^^^
      fprintf('\n')
      fprintf(['The estimated parameter correlations are stored in the matrix ',corname,'\n'])
      a = input('Do you want to view them now [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y'
          fprintf(['The parameters appear in ',corname,' in the order\n'])
          fprintf('first-year survival, adult survival, productivity rate,\n')
          fprintf('measurement error st. deviation\n')
          prtupmat(C,2);
        end
      end
      fprintf('\n')

%     Print max loglik and AIC
%     ^^^^^^^^^^^^^^^^^^^^^^^^
      likname = ['lik_',modelname];
      eval([likname,' = -alikj(',betaname,',',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
      fprintf('The maximized loglikelihood is %-10.4f\n',eval(likname))
        
      AIC = 2*(npars-eval(likname));
      fprintf('The Akaike information criterion is %-10.3f\n\n',AIC)
      fprintf('\n')
     
%     Print smoothed estimates and fitted values
%     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      smoothname = ['smooth_',modelname];
      eval([smoothname,' = ksc(',betaname,',',X1name,',',Xaname,',',Xpname,',''',phiatype,''');'])
      fprintf(['The smoothed estimates of the population numbers are stored in the matrix ',smoothname,'\n'])
      a = input('Do you want to view them now [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y'
          fprintf(['The estimates appear in ',smoothname,', in columns, in the order\n'])
          if modelaname(1) == 'C' | modelaname(1) == 'T',
            fprintf('[1-year olds | adults]\n')
          else
            fprintf('[1-year olds | 2-year olds | ... | adults]\n')
          end
          eval(['prtupmat(',smoothname,',1);'])
        end
      end
      fprintf('\n')

      a = input('Do you want to plot the observed vs the fitted census [N] ? ','s');
      if ~isempty(a),
        if a(1)=='y' | a(1)=='Y',
          plot(Census,':')
          hold on
          eval(['plot(sum(',smoothname,'(:,2:end),2));'])
          legend('observed','fitted')
          hold off
        end
      end
      fprintf('\n')
     

%    elseif dofit == 2,
%      % Start of score test
%      % ^^^^^^^^^^^^^^^^^^^
%      % Input the null model
%      % ^^^^^^^^^^^^^^^^^^^^
%      prefix = 'null';
%      inputmn;
%      nm = modelname;
%      nm1 = model1name;
%      nma = modelaname;
%      nmlam = modellamname;
%      parname = ['b_',nm];
%      pnl = min(length(parname),matlab_varlength);
%      if exist(parname(1:pnl)),
%        nparnull1 = length(eval(['X1_',nm1,'(1,:)']));
%        nparnulla = length(eval(['Xa_',nma,'(1,:)']));
%        nparnulllam = length(eval(['Xlam_',nmlam,'(1,:)']));
%        nparnull = nparnull1 + nparnulla + nparnulllam;
%    
%        % Input the alternative model
%        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^
%        prefix = 'alt';
%        inputmn;
%        am = modelname;
%        am1 = model1name;
%        ama = modelaname;
%        amlam = modellamname;
%        nparalt1 = length(eval(['X1_',am1,'(1,:)']));
%        nparalta = length(eval(['Xa_',ama,'(1,:)']));
%        nparaltlam = length(eval(['Xlam_',amlam,'(1,:)']));
%        nparalt = nparalt1 + nparalta + nparaltlam;
%      
%        % Do the score test
%        % ^^^^^^^^^^^^^^^^^
%        starttime = clock;
%        eval(['scorevalue=scoreJT(b_',nm,',X1_',nm1,',Xa_',nma,',Xlam_',nmlam,',X1_',am1,',Xa_',ama,',Xlam_',amlam,',''',phiatype,''');'])
%        runtime = etime(clock,starttime);
%        fprintf('\ntime taken = %5.1f secs\n',runtime)
%        df = nparalt - nparnull;
%        p_val = pchisq(scorevalue,df);
%        fprintf(['\nThe score statistic for testing the alternative model ',am,'\n'])
%        fprintf(['against the null model ',nm,'\n'])
%        fprintf('is %8.4f on %2.0f degrees of freedom     (p-val = %8.5g)\n\n',scorevalue,df,p_val)
%
%    else % if not exist beta_nullmodelname
%
%      fprintf('You must fit the null model before you can use it in a score test.\n\n')    
%
%    end % if not exist beta_nullmodelname
%
  end %  if dofit == 2
  
end % endless while loop
 
end % type == 2
