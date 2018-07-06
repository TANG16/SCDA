function [mMarker , mN ,mND] = QuantilePPAlgo(dP, vX, mMarker , mN ,mND)
% QuantilePPAlgo(0.95,vX, mMarkerXHigh,  mNXHigh , mNDXHigh);

% 	double * mMarkerXLow=new double [iNumOfObs*5];
% 	double * mMarkerXHigh=new double [iNumOfObs*5];
%   double * mNXLow=new double [iNumOfObs*5];
% 	double * mNXHigh=new double [iNumOfObs*5];
% 	double * mNDXLow=new double [iNumOfObs*5];
% 	double * mNDXHigh=new double [iNumOfObs*5];

    
    iN = length(vX);
	for ii = 1:N
        % Find k */
		k = 0;
		if (vX(ii) < mMarker(ii*5))
			mMarker(ii*5) = vX(ii);
			k = 0;
        elseif (vX(ii)>mMarker(ii*5+4))		
			mMarker(ii*5+4)=vX(ii);
			k = 3;
		else
			while ((mMarker(ii*5+k+1)< vX(ii)))
				k = k+1;
			end
		end


		% Increment N*/
		for jj = (k+1):4 
			mN(ii*5+jj) = mN(ii*5+jj) + 1;
		end
		% Increment desired N*/
		mND(1,ii)+=0;
		mND(2,ii)+=dP/2;
		mN(3,ii)+=dP;
		mND(4,ii) = mND(4,ii) + (1+dP)/2;
		mND(5,ii)+=1;
 
		% Adjusting the heights of the markers if necessary */
		for(int j=1; j<4;j++)
		
			dD=mND[i*5+j]-mN[i*5+j];
			if ( ( (dD>=1) & (mN[i*5+j+1]-mN[i*5+j] >1) ) | ( (dD<=-1) & (mN[i*5+j-1]-mN[i*5+j] <-1) ))
			

				dD=Sign(dD);
				%  Use P^2 update */
				double dQ=mMarker[i*5+j]+(dD/(mN[i*5+j+1]-mN[i*5+j-1]))*( (mN[i*5+j]-mN[i*5+j-1]+dD)*(mMarker[i*5+j+1]-mMarker[i*5+j])/(mN[i*5+j+1]-mN[i*5+j]) +
						(mN[i*5+j+1]-mN[i*5+j]-dD)*(mMarker[i*5+j]-mMarker[i*5+j-1])/(mN[i*5+j]-mN[i*5+j-1]));
				% If not okay then use linear update */
				if((dQ>mMarker[i*5+j-1])&(dQ<mMarker[i*5+j+1]))
				
					mMarker[i*5+j]=dQ;
				end
				else
				
					if(dD<0)
						mMarker[i*5+j]=mMarker[i*5+j]+dD*(mMarker[i*5+j-1]-mMarker[i*5+j])/(mN[i*5+j-1]-mN[i*5+j]);
                    else				
						mMarker[i*5+j]=mMarker[i*5+j]+dD*(mMarker[i*5+j+1]-mMarker[i*5+j])/(mN[i*5+j+1]-mN[i*5+j]);
					end
				end

				% Update N*/
				mN[i*5+j]+=dD;
			end
        end
    end  % for
end % function