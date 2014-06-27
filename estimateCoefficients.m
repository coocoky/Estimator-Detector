
function  estimateCoefficients(sample, n,sampBforWindow)
                                            
global DensityVars;
                                        
% For the first iteration of the online density estimator, the coefficients 
% are usually set to zeros. 



% Check to see if the sample is outside the domain. If the sample is 
% outside the domain, don't count it.
inRange = checkRange(sample, DensityVars.DensityDomain);
if( ~inRange )
    return
end % if( ~checkRange(sample(k), DensityVars.DensityDomain) )


% Upper and lower values of wavelet basis support.
lowerSupp = min(DensityVars.WaveSupport);
upperSupp = max(DensityVars.WaveSupport);

% Scaling coefficients only.
scalCoeff = DensityVars.Coefficients{1,1};    


% Corresponding coefficients and translates for the scaling function with 
% current sample.
[ correspCoeff,...
  relevantKs,...
  lowerKIndex,...
  upperKIndex] = findRelevantCoefficients(sample,...
                                          lowerSupp,...
                                          upperSupp,...
                                          DensityVars.StartLevel,...
                                          scalCoeff,...
                                          DensityVars.ScalTranslates);

% Perform interpolation across each translate and sample.
x = 2^DensityVars.StartLevel*sample - relevantKs;
phiHere = getPhiAt(x);

% Updating relevant coefficients.
% If agingFlag == 2 then updating using the window method is performed.
if (DensityVars.AgingFlag == 2)
    % Update corresponding coefficients with addition as default.
    updCoeffAdd = correspCoeff + 2^(DensityVars.StartLevel / 2)*(phiHere / DensityVars.WindowSize);
    scalCoeff(lowerKIndex:upperKIndex) = updCoeffAdd;
    
        % Check to see if n > windowSize. If n > windowSize then begin
        % discounting the sampleBeforeWindow's corresponding
        % coefficients.
        if (n > DensityVars.WindowSize)
            % Corresponding coefficients and translates for the scaling function 
            % the sampleBeforeWindow.
            [corrCoefSampBfWin,...
             relKsSampBfWin,...
             lowKIndSampBfWin,...
             uppKIndSampBfWin] = findRelevantCoefficients(sampBforWindow,...
                                                          lowerSupp,...
                                                          upperSupp,...
                                                          DensityVars.StartLevel,...
                                                          scalCoeff,...
                                                          DensityVars.ScalTranslates); 

            % Perform interpolation across each translate and sample for sampleBforWindow.
            xSampBfWin = 2^DensityVars.StartLevel*sampBforWindow - relKsSampBfWin;
            phiHereSampBfWin = getPhiAt(xSampBfWin);
            updCoeffSubtr = corrCoefSampBfWin - 2^(DensityVars.StartLevel / 2)*(phiHereSampBfWin / DensityVars.WindowSize);
            
            scalCoeff(lowKIndSampBfWin : uppKIndSampBfWin) = updCoeffSubtr;
            

            DensityVars.Coefficients{1,1} = scalCoeff;                
        else
            scalCoeff(lowerKIndex : upperKIndex) = updCoeffAdd;
            DensityVars.Coefficients{1,1} = scalCoeff;
        end % if (n > windowSize)
else        
    updCoeff = (n / (n + 1))*correspCoeff + (1 / (n + 1))*2^(DensityVars.StartLevel / 2)*phiHere;
    allCoeffUpd = (n / (n + 1))*scalCoeff;
    allCoeffUpd(lowerKIndex : upperKIndex) = updCoeff;
    DensityVars.Coefficients{1,1} = allCoeffUpd;   
end % (agingFlag == 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVELETS.
if (DensityVars.WaveletFlag == 1)
    for j = DensityVars.StartLevel : DensityVars.StopLevel   
        % Translates for resolution j.
        transOfResJ = DensityVars.WaveTranslates{((j - DensityVars.StartLevel) + 1),1};

        % Current j coefficients.
        coeffIndVal = ((j - DensityVars.StartLevel) + 2) ; % Only wavelet coefficients
        currWaveCoeff = DensityVars.Coefficients{coeffIndVal,1};

        % Corresponding coefficients and translates for the wavelet function with 
        % current sample.
        [corresCoeffWavCurrSamp,...
              relKsWavCurrSamp,...
              lowKIndexWavCurrSamp,...
              uppKIndexWavCurrSamp] = findRelevantCoefficients(sample,...
                                                      lowerSupp,...
                                                      upperSupp,...
                                                      j,...
                                                      currWaveCoeff,...
                                                      transOfResJ);

        % Perform interpolation with the psi wavelet for the current sample.
        xWavCurrSamp = 2^j*sample - relKsWavCurrSamp;
        psiHereWavCurrSamp = getPsiAt( xWavCurrSamp );
        
        % If agingFlag == 2 then updating using the window method is performed.
        if (DensityVars.AgingFlag == 2)
            % Update corresponding wavelet coefficients with addition as default.
            updWavCoeffAddCurrSamp = corresCoeffWavCurrSamp + 2^(j / 2)*(psiHereWavCurrSamp / DensityVars.WindowSize);
            currWaveCoeff(lowKIndexWavCurrSamp:uppKIndexWavCurrSamp) = updWavCoeffAddCurrSamp;
            % Check to see if n > windowSize. If n > windowSize then begin
            % discounting the sampleBeforeWindow's corresponding
            % coefficients.
            if (n > DensityVars.WindowSize)
                % Corresponding coefficients and translates for the wavelet function 
                % the sampleBeforeWindow.
                [corrCoefWavSampBfWin,...
                 relKsWavSampBfWin,...
                 lowKIndWavSampBfWin,...
                 uppKIndWavSampBfWin] = findRelevantCoefficients(sampBforWindow,...
                                                              lowerSupp,...
                                                              upperSupp,...
                                                              j,...
                                                              currWaveCoeff,...
                                                              transOfResJ); 

                % Perform interpolation across each translate and sample for sampleBforWindow.
                xWavSampBfWin = 2^j*sampBforWindow - relKsWavSampBfWin;
                psiHereWavSampBfWin = getPsiAt(xWavSampBfWin);
                
                % Subtraction portion for the sample before window.                    
                updWavCoeffSubtr = corrCoefWavSampBfWin - 2^(j / 2)*(psiHereWavSampBfWin / DensityVars.WindowSize);
                currWaveCoeff(lowKIndWavSampBfWin:uppKIndWavSampBfWin) = updWavCoeffSubtr;
                

                  DensityVars.Coefficients{coeffIndVal,1} = currWaveCoeff;         
            else
                currWaveCoeff(lowKIndexWavCurrSamp : uppKIndexWavCurrSamp) = updWavCoeffAddCurrSamp;
                DensityVars.Coefficients{coeffIndVal,1} = currWaveCoeff;
            end % if (n > windowSize)

        else
            % Update the relevant coefficients.
            updWaveCoeff = (n / (n + 1))*corresCoeffWavCurrSamp + (1 / (n + 1))*2^(j / 2)*psiHereWavCurrSamp;
            currWaveCoeff = (n / (n + 1))*currWaveCoeff;
            currWaveCoeff(lowKIndexWavCurrSamp : uppKIndexWavCurrSamp) = updWaveCoeff;
            DensityVars.Coefficients{coeffIndVal,1} = currWaveCoeff;
        end % if (agingFlag == 2)
    end  % for j = startLevel : stopLevel 
end % if (waveletFlag == 1)     
       

end % estimateCoefficients()











