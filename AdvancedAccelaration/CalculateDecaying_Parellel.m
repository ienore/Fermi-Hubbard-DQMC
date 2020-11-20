core_number = 8;
T_hop = 1.0;
NumInEdge = 8;
NumOfWarm = 1;
NumOfWarm_inside = 0;
OutputNumber = 8;
NumOfCores = core_number;
NumOfEpochInBin = 1;
N_wrap = 10;
Uene = 0.0;
Beta = 8;
D_Tau = 0.2;
NumOfVertexs = NumInEdge^2;
TempSlice = Beta/D_Tau;
propa_green_group = zeros([NumOfVertexs,NumOfVertexs,TempSlice,OutputNumber,NumOfCores]);
propa_green_sample = zeros([NumOfVertexs,NumOfVertexs,TempSlice,OutputNumber*NumOfCores]);
parfor zjy_index = 1:1:NumOfCores
    propa_green_group(:,:,:,:,zjy_index) = CalculateDecaying_ParellelKernel(zjy_index,Beta,D_Tau,NumOfWarm,NumOfWarm_inside,N_wrap,Uene,T_hop,NumInEdge,OutputNumber,NumOfEpochInBin);
end
% for zjy_index = 1:1:NumOfCores
%     for pa_index = 1:1:OutputNumber
%         group_index = (zjy_index-1) * OutputNumber + pa_index;
%         propa_green_sample(:,:,:,group_index) = propa_green_group(:,:,:,pa_index,zjy_index);
%     end
% end
