function [ sourcePos ] = computeSourcePosDeformationFRS( srcTfrs, VecIEC )
%COMPUTESOURCEPOSDEFORMATIONFRS computes position of X-ray source in FRS CS
%from transformation matrix and deformation vector (X,Y,Z) given in local
%IEC CS
%
% srcTfrs (4x4 matrix)
% VecIEC ( 3-element vector: x,y,z)

sourcePos = srcTfrs*trans(VecIEC(1),VecIEC(2),VecIEC(3))*[0,0,0,1]'; 
sourcePos = sourcePos(1:3)';

end
