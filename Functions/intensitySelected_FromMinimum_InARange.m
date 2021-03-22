
function [mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph, mRNA_mRNA,n_minima,mRNA_Ser5ph,mRNA_Unph] = intensitySelected_FromMinimum_InARange (mRNAFluc_Rave,UnphFluc_Rave,Ser5phFluc_Rave,vector_minima_mRNA, withRespectToMinValue, PlusMinusRange)
%% Section that returns the intensity values in a range +/- N seconds with respect to the minimal values of

% INPUT:
%      - mRNAFluc_Rave,UnphFluc_Rave,Ser5phFluc_Rave are the complete temporal intensities
%      for mRNA, Ser5ph and Unph. Matrices (200 x 7), 200 time points, 7
%      experimental cells.
%      - selectedPairs: is a cell (No_Experimental_cells , maxNumberOfPairs ) that contains
%      the minima value of two signals.

%      - withRespectToMinValue: is a string that represents the the
%      signal that is selected as the central element of the selection. the
%      accepted values mRNA, Ser5ph and Unph.
%      PlusMinusRange: is the range that will be selected to report.

% OUTPUT:
%       - Intensities for mRNA, Ser5ph, Unph in the selected local minima within a PlusMinusRange
%       - int_mRNA, int_Ser5ph, int_Unph. Each one is a cell (No_Experimental_cells , maxNumberOfPairs) and each element in the cell containts a vector with the intenisty


switch withRespectToMinValue
    case 'mRNA'
        indexWithRespect =1;
    case 'Unph'
        indexWithRespect =2;
    case 'Ser5ph'
        indexWithRespect =3;
end
% 
% for i =1: size(selectedPairs,1)
%     if isempty( selectedPairs{i}) ==0
%         for k =1:size( selectedPairs,2)
%             if isempty(selectedPairs{i,k})==0
%                 int_mRNA{i,k} = mRNAFluc_Rave( selectedPairs{i,k}(indexWithRespect)-PlusMinusRange: selectedPairs{i,k}(indexWithRespect)+PlusMinusRange ,i);
%                 int_Ser5ph{i,k} = Ser5phFluc_Rave( selectedPairs{i,k}(indexWithRespect)-PlusMinusRange: selectedPairs{i,k}(indexWithRespect)+PlusMinusRange ,i);
%                 int_Unph{i,k} = UnphFluc_Rave( selectedPairs{i,k}(indexWithRespect)-PlusMinusRange: selectedPairs{i,k}(indexWithRespect)+PlusMinusRange ,i);
%             end
%         end
%     else
%         int_mRNA{i,k} =[];
%         int_Ser5ph{i,k} =[];
%         int_Unph{i,k} =[];
%     end
% end
mRNA_Unph =[]; mRNA_mRNA=[]; mRNA_Ser5ph =[];
NoSelectedCells = 0;
for i =1:size(vector_minima_mRNA,2)
    
    if isempty(vector_minima_mRNA{i})==0
        NoSelectedCells = NoSelectedCells+1;
        for k =1:length( vector_minima_mRNA{i})
                index = vector_minima_mRNA{i}(k);
%                 int_mRNA{i,k} = mRNAFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i);
%                 int_Ser5ph{i,k} = Ser5phFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i);
%                 int_Unph{i,k} = UnphFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i);
                mRNA_mRNA =   [mRNA_mRNA,mRNAFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i)];
                mRNA_Unph =   [mRNA_Unph,UnphFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i)];
                mRNA_Ser5ph = [mRNA_Ser5ph,Ser5phFluc_Rave( index-PlusMinusRange: index+PlusMinusRange ,i)];
        
        end
                
    else
        continue
%         int_mRNA{i,1} =[];
%         int_Ser5ph{i,1} =[];
%         int_Unph{i,1} =[];
    end
end
   
    


%% creates a table with the intensities for each signal.
% this table has 200 rows (time points) and N columns that represent the
% toal number of minimal values.


% NoSelectedCells=0;
% for i =1:size(vector_minima_mRNA,2)
%     
%     if isempty( selectedPairs{i}) ==0
%         NoSelectedCells= NoSelectedCells+1;
%         for k =1:size( selectedPairs,2)
%             % For mRNA
%             mRNA_mRNA =   [mRNA_mRNA,int_mRNA{i,k}];
%             mRNA_Unph =   [mRNA_Unph,int_Unph{i,k}];
%             mRNA_Ser5ph = [mRNA_Ser5ph,int_Ser5ph{i,k}];
%         end
%     end
% end

%% Calculating Averafge quantities

% Calculating mean values
mean_mRNA = mean (mRNA_mRNA,2)';
mean_Ser5ph = mean (mRNA_Ser5ph,2)';
mean_Unph = mean (mRNA_Unph,2)';

% calculating SEM
n_minima = size(mRNA_mRNA,2);
sem_mRNA = std (mRNA_mRNA,0,2)'/sqrt(n_minima) ;
sem_Ser5ph = std (mRNA_Ser5ph,0,2)'/sqrt(n_minima);
sem_Unph = std (mRNA_Unph,0,2)'/sqrt(n_minima);

n_minima = size(mRNA_mRNA,2);

end
