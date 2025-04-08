#include "PQPIKS/DynamicModel.h"

namespace PQPIKS
{

DynamicModel::DynamicModel( int DOF, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity ) : DOF(DOF), AllRotAxes(AllRotAxes), Masses(Masses), AllFirstMomentsOfMass(AllFirstMomentsOfMass), AllInertias(AllInertias), Gravity(Gravity)
{
}

DynamicModel::DynamicModel( const DynamicModel &other )
: DOF(other.DOF), AllRotAxes(other.AllRotAxes), Masses(other.Masses), AllFirstMomentsOfMass(other.AllFirstMomentsOfMass), AllInertias(other.AllInertias), Gravity(other.Gravity)
{
}

DynamicModel& DynamicModel::operator= ( const DynamicModel &other )
{
	if (&other == this) // self-assignment
	{
		return *this;
	}
	
	DOF = other.DOF;
	AllRotAxes = other.AllRotAxes;
	Masses = other.Masses;
	AllFirstMomentsOfMass = other.AllFirstMomentsOfMass;
	AllInertias = other.AllInertias;
	Gravity = other.Gravity;
	
	return *this;
}

DynamicModel::~DynamicModel()
{
}

void DynamicModel::InertiaMatrix( MatrixD &M, const MatrixD &JointPos )
{
	//~ M = 2.0 * MatrixD::Identity(DOF, DOF);
	
	MatrixD PhiH = MatrixD::Zero(6 * DOF, DOF);
	MatrixD MassMatrix = MatrixD::Zero(6 * DOF, 6 * DOF);
	
	MatrixD T_i = MatrixD::Identity(4, 4);	// transformation matrix from the base frame to the i-th link frame (currently, before the loop, it's the transformation from the base to the base, i.e. identity)
	MatrixD T_j(4, 4);	// transformation matrix from the base frame to the j-th link frame
	MatrixD T_ij(4, 4);	// transformation matrix from the i-th link frame to j-th link frame
	Eigen::Vector3d RotAxis;	// unit vector along rotational axis of the j-th joint expressed in the j-th link frame
	
	
	for (int j = 0; j < DOF; j++)
	{
		TransformationMatrix( T_ij, JointPos, j );
		T_j = T_i * T_ij;
		
		RotAxis = T_j.block<3,3>(0,0) * AllRotAxes.col(j);
		
		MatrixD H_j = MatrixD::Zero(6, 1);
		H_j.block<3, 1>(3, 0) = RotAxis;
		
		PhiH.block(6 * j, j, 6, 1) = H_j;
		
		if ( j > 0 )
		{
			int i = j - 1;
			
			MatrixD r_ij = T_j.block<3,1>(0,3) - T_i.block<3,1>(0,3);	// vector linking the j-th and i-th link frame origins, expressed in the base frame
			
			MatrixD B_ij = MatrixD::Identity(6, 6);
			B_ij.block<3,3>(0,3) = -skew(r_ij);
			
			PhiH.block(6 * j, 0, 6, j) = B_ij * PhiH.block(6 * i, 0, 6, j);
		}
		
		MatrixD Inertia = T_j.block<3,3>(0,0) * AllInertias.block<3,3>(0, 3 * j) * (T_j.block<3,3>(0,0)).transpose();
		MatrixD FirstMomentOfMass = T_j.block<3,3>(0,0) * AllFirstMomentsOfMass.col(j);
		
		MassMatrix.block<3,3>(6 * j, 6 * j) = Masses(j) * MatrixD::Identity(3, 3);
		MassMatrix.block<3,3>(6 * j, 6 * j + 3) = -skew(FirstMomentOfMass);
		MassMatrix.block<3,3>(6 * j + 3, 6 * j) = skew(FirstMomentOfMass);
		MassMatrix.block<3,3>(6 * j + 3, 6 * j + 3) = Inertia;
		
		T_i = T_j;
	}
	MatrixD PhiHtranspose = PhiH.transpose();
	
	//~ Manipulator inertia matrix
	M = PhiHtranspose * MassMatrix * PhiH;
}

MatrixD DynamicModel::skew( const Eigen::Vector3d &Vector )
{
//	Returns a skew-symmetric matrix from a vector.
	MatrixD SkewSymmetricMatrix(3, 3);

	SkewSymmetricMatrix << 0.0, -Vector(2), Vector(1),
			Vector(2), 0.0, -Vector(0),
			-Vector(1), Vector(0), 0.0;

	return SkewSymmetricMatrix;
}

} //namespace PQPIKS
