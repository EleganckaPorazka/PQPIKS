#ifndef DYNAMICMODEL_H_
#define DYNAMICMODEL_H_

#include <iostream>
#include <Eigen/Dense>
#include "typedefs.h"

namespace PQPIKS
{
	
class DynamicModel
{	
public:
	DynamicModel( int DOF, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity );
	
	DynamicModel( const DynamicModel &other );
	
	DynamicModel& operator= ( const DynamicModel &other );
	
	virtual ~DynamicModel();
	
	virtual void TransformationMatrix( MatrixD &T_ij, const Eigen::VectorXd &JointPos, int j ) = 0;
	
	void InertiaMatrix( MatrixD &M, const MatrixD &JointPos );
	
	MatrixD skew( const Eigen::Vector3d &Vector );

protected:
	int DOF; // number of degrees of freedom
	MatrixD AllRotAxes;
	MatrixD Masses;
	MatrixD AllFirstMomentsOfMass;
	MatrixD AllInertias;
	MatrixD Gravity;
};

} //namespace PQPIKS

#endif /* DYNAMICMODEL_H_ */
