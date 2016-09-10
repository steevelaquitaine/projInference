
 %author: Steeve Laquitaine
   %date: 140304
%purpose: simulate a vector of voxels BOLD response from neurons firing
%usage: BOLD=makeVoxels(1:1:15,rand(15,1))
        %e.g., "whichNeurInaVoxel" is a matrix (m voxels, 1000 neurons) 
        %indicating which of 15 neurons composed each voxel.
        
        %BOLD is a column vector of Bold response for 500 voxels.       
        
%Descriptions: we average firing with a voxel instead of summing them 
%because too many neurons leads to very high BOLD and this affects forward 
%encoding modeling for example. It's just scaling the data.
        
function BOLD=makeVoxels(whichNeurInaVoxel,NeuronsFiring)

%check that NeuronsFiring is a column vector
if size(NeuronsFiring,1)<size(NeuronsFiring,2)
    NeuronsFiring=NeuronsFiring';
end

%make Bold responses of each voxels by summing neurons' firing
%FiringInaVoxel is a matrix (N neurons,m voxels)
NeuronsFiring=NeuronsFiring';
FiringInaVoxel=NeuronsFiring(whichNeurInaVoxel);
BOLD=mean(FiringInaVoxel,1);
BOLD=BOLD';





