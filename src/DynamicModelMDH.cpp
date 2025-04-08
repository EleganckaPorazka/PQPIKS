#include "PQPIKS/DynamicModelMDH.h"

namespace PQPIKS
{
DynamicModelMDH::DynamicModelMDH( int DOF, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ )
: DynamicModel( DOF, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity ), MDHangles(MDHangles), MDHlengthsX(MDHlengthsX), MDHlengthsZ(MDHlengthsZ)
{
}

DynamicModelMDH::DynamicModelMDH( const DynamicModelMDH &other )
: DynamicModel(other), MDHangles(other.MDHangles), MDHlengthsX(other.MDHlengthsX), MDHlengthsZ(other.MDHlengthsZ)
{
}

DynamicModelMDH& DynamicModelMDH::operator= ( const DynamicModelMDH &other )
{
	if (&other == this)
	{
		return *this;
	}
	
	DynamicModel::operator= (other);
	
	MDHangles = other.MDHangles;
	MDHlengthsX = other.MDHlengthsX;
	MDHlengthsZ = other.MDHlengthsZ;
	
	return *this;
}

DynamicModelMDH::~DynamicModelMDH()
{
}

void DynamicModelMDH::TransformationMatrix( MatrixD &T_ij, const Eigen::VectorXd &JointPos, int j )
{
//	Computes the transformation matrix from the previous link frame to the current link frame.

	T_ij << cos( JointPos(j) ), -sin( JointPos(j) ), 0.0, MDHlengthsX(j),
			sin( JointPos(j) ) * cos( MDHangles(j) ), cos( JointPos(j) ) * cos( MDHangles(j) ), -sin( MDHangles(j) ), -MDHlengthsZ(j) * sin( MDHangles(j) ),
			sin( JointPos(j) ) * sin( MDHangles(j) ), cos( JointPos(j) ) * sin( MDHangles(j) ), cos( MDHangles(j) ), MDHlengthsZ(j) * cos( MDHangles(j) ),
			0.0, 0.0, 0.0, 1.0;
}

} //namespace PQPIKS
