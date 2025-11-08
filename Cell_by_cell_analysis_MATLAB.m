%%%%%%%%%%%%%%%%%%%%%%INPUT!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create 2 variables in Matlab and input the data
% These variables MUST be named "pe", and "none" for this code to work

% The first two rows of "pe" and "none" have to be: slope and interceptor
% of the PE beads calibration curve

% The i-th column of "pe" and the i-th column of "none" should have the
% same slope and intercept

%%% You can change the number of bins here for the exported histogram
numhistbins = 120; % how many bins to represent the histogram? 
% Check exported figure and see if the curve is intact from left to right,
% if not then increase numhistbins
% If the curve is not very smooth then decrease numhistbins

%%% EXPORT!! %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% median, c.v., kurtosis and skewness of the cell-by-cell distribution
%%% EXPORT(x,y) is the histogram of the distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%CAUTION!!%%%%%%%%%%%%%%%%
% type "clear all" and enter in the command window 
% before you add another pair of "none" and "pe"

%After each run, type in "save FILENAME" inthe command window
%For example, save HDF_R1_CONTROL
%or, save HUVEC_N1_PDGFAA5
%Then type in "clear all" to release the memory space before you start
%another analysis

%If you need to revisit your analysis, i.e. EXPORT, none, pe, just type in
%"load FILENAME"

%% Don't alter anything below this line 


slope_b=none(1,:);
intercept_b=none(2,:);

for j=1:length(none(1,:))
    temp = none(3:end,j);
    temp = temp(temp>0);
    S2(j)=sum(temp)./length(temp);
    clear temp;
end

slope=pe(1,:);
intercept=pe(2,:);

for j=1:length(pe(1,:))
    temp = pe(3:end,j);
    temp = temp(temp>0);
    O1(j)=sum(temp)./length(temp);
    clear temp;
end

K=zeros();
K=S2./O1;
       
for j=1:length(pe(1,:))
    for i=3: length(pe(:,j))
     
        pe(i,j)=pe(i,j)*(1-K(j));
    
        PE(i-2,j)=10.^((log10(pe(i,j))-intercept(j))./slope(j)); 
     
    end
    
    PE5000 (:,j)=datasample(PE(PE~=0),5000);%make the sample size equal before pooling.
end


PE2 = reshape(PE5000,numel(PE5000),1);
PE2=PE2(PE2>0);
Median=median(PE2);
Coeff_variance=std(PE2)./mean(PE2);
Kurtosis=kurtosis(PE2);
Skewness=skewness(PE2);
CellNum=length(PE2);
IQR=iqr(PE2);
mu=mean(PE2);
gmu=geomean(PE2);
sigma=std(PE2);


PE3 = sort(log(PE2));

[A,S] = hist(PE3,numhistbins); % S returns position of the bin centers in X.
    A=A';S=S';
    for k = 1:length(A)
        if A(k) < 0
        A(k) = 0;
        end
    end
    x = A;
    x(1) = (x(1) + x(2) + x(3))/3;
    x(2) = (x(1) + x(2) + x(3) + x(4))/4;
    for v = 3:(length(x)-2)
        x(v) = (x(v-2) + x(v-1) + x(v) + x(v+1) + x(v+2))/5;
    end
    x(length(x)-1) = (x(length(x)-3) + x(length(x)-2) + x(length(x)-1) + x(length(x)))/4;
    x(length(x)) = (x(length(x)-2) + x(length(x)-1) + x(length(x)))/3;
    x = x./sum(x);

    S1 = exp(S);
    semilogx(S1,x); 
    EXPORT(:,1) = S1; EXPORT(:,2) = x; 
     
    g1=exp(-(S-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    g1 = g1./sum(g1);
    [h1,p1]=kstest2(x,g1,'Alpha',0.05);
    p1
    
 
    qe=0;
    
    for i=1:numhistbins-1
        for j=(i+1):numhistbins
            qe=qe+x(i)*x(j)*(S(j)-S(i));
        end
    end
    
   qe;        
 T=table({'CellNum';'Median';'GeoMean'; 'CV';'Skewness'; 'QE';'Kurtosis';'IQR';'Mean';'SD';},[CellNum;Median;gmu;Coeff_variance;Skewness;qe;Kurtosis;IQR;mu;sigma])  
   %% Determining the cutoff of the data
PE2 = sort(PE2);
PE2 = log(PE2);
CutoffThreshold = 0.01;
[bincountsPE2,bincentersPE2]=hist(PE2,numhistbins);
maxhit = max(bincountsPE2);
count = 1;
for i = 2:(length(bincentersPE2)-1)
  if (bincountsPE2(i) <= CutoffThreshold*maxhit) && (bincountsPE2(i-1) <= CutoffThreshold*maxhit) && (bincountsPE2(i+1) <= CutoffThreshold*maxhit)
      PE2_Cutoff(count) = bincountsPE2(i);
      count = count +1;
  end
end        
PercentOutlier = 100*sum(PE2_Cutoff)/sum(bincountsPE2);

save ('/OUTPUT_FOLDER/FILE_NAME.mat');

%% funtion: write with headers
% This function functions like the build in MATLAB function csvwrite but
% allows a row of headers to be easily inserted
%
% known limitations
% 	The same limitation that apply to the data structure that exist with 
%   csvwrite apply in this function, notably:
%       m must not be a cell array
%
% Inputs
%   
%   filename    - Output filename
%   m           - array of data
%   headers     - a cell array of strings containing the column headers. 
%                 The length must be the same as the number of columns in m.
%   r           - row offset of the data (optional parameter)
%   c           - column offset of the data (optional parameter)
%
%
% Outputs
%   None
function csvwrite_with_headers(filename,m,headers,r,c)
%% initial checks on the inputs
if ~ischar(filename)
    error('FILENAME must be a string');
end
% the r and c inputs are optional and need to be filled in if they are
% missing
if nargin < 4
    r = 0;
end
if nargin < 5
    c = 0;FD_F
end
if ~iscellstr(headers)
    error('Header must be cell array of strings')
end
 
if length(headers) ~= size(m,2)
    error('number of header entries must match the number of columns in the data')
end
%% write the header string to the file 
%turn the headers into a single comma seperated string if it is a cell
%array, 
header_string = headers{1};
for i = 2:length(headers)
    header_string = [header_string,',',headers{i}];
end
%if the data has an offset shifting it right then blank commas must
%be inserted to match
if r>0
    for i=1:r
        header_string = [',',header_string];
    end
end
%write the string to a file
fid = fopen(filename,'w');
fprintf(fid,'%s\r\n',header_string);
fclose(fid);
%% write the append the data to the file
%
% Call dlmwrite with a comma as the delimiter
%
dlmwrite(filename, m,'-append','delimiter',',','roffset', r,'coffset',c);

%% Output save 
headers_1 = {'X','Y'};
csvwrite_with_headers('/OUTPUT_FOLDER/FILE_NAME_Plot.csv',EXPORT,headers_1)
headers_2 = {'1'};
csvwrite_with_headers('/OUTPUT_FOLDER/FILE_NAME_PE5000.csv',PE5000,headers_2)
writetable(T,'/OUTPUT_FOLDER/FILE_NAME_stats.csv','Delimiter',',','QuoteStrings',true)
end
