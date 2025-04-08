#ifndef DYNAMICMODELMDH_H_
#define DYNAMICMODELMDH_H_

#include "DynamicModel.h"

namespace PQPIKS
{

class DynamicModelMDH : public DynamicModel
{
public:
	DynamicModelMDH( int DOF, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ );
	
	DynamicModelMDH( const DynamicModelMDH &other );
	
	DynamicModelMDH& operator= ( const DynamicModelMDH &other );
	
	virtual ~DynamicModelMDH();
	
	virtual void TransformationMatrix( MatrixD &T_ij, const Eigen::VectorXd &JointPos, int j );

protected:
	Eigen::VectorXd MDHangles;
	Eigen::VectorXd MDHlengthsX;
	Eigen::VectorXd MDHlengthsZ;
};

} //namespace PQPIKS

#endif /* DYNAMICMODELMDH_H_ */
