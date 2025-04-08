#ifndef KINEMATICMODELMDH_H_
#define KINEMATICMODELMDH_H_

// Provides a forward kinematics solution and a jacobian, using modified Denavit-Hartenberg (MDH) parameters
// Repurposed from the RedundantRobot project

#include <iostream>
#include <Eigen/Dense>
#include "typedefs.h"

//~ typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixD;

namespace PQPIKS
{

class KinematicModelMDH
{
public:
//	Constructor without an EE
	KinematicModelMDH( int DOF, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, int *ObstacleAvoidanceDOF, double ObstacleAvoidanceNominalVel, double ObstacleSphereOfInfluence, double ObstacleUnityGain, const MatrixD &ObstacleAvoidancePoints );

//	Full parameterized constructor
	KinematicModelMDH( int DOF, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE, int *ObstacleAvoidanceDOF, double ObstacleAvoidanceNominalVel, double ObstacleSphereOfInfluence, double ObstacleUnityGain, const MatrixD &ObstacleAvoidancePoints );

	KinematicModelMDH( const KinematicModelMDH &other );
	
	KinematicModelMDH& operator= ( const KinematicModelMDH &other );
	
	virtual ~KinematicModelMDH();

//	Forward kinematics
	void ForwardKinematics( MatrixD &T, const MatrixD &JointPos );

	void GetJacobian( MatrixD &J, const MatrixD &JointPos ); // computes the Jacobian of the manipulator
	//~ virtual Eigen::VectorXd GetEndEffectorVel( const Eigen::VectorXd &JointPos, const Eigen::VectorXd &NextEndEffectorPos, double StepSize ); // computes the desired end effector velocity
	//~ virtual Eigen::VectorXd GetEndEffectorError( const Eigen::VectorXd &JointPos, const Eigen::VectorXd &EndEffectorPos ); // computes the end effector error
	void GetObstacleJacobian( MatrixD &Jo, int CriticalPointNumber, const MatrixD &JointPos ); // computes the Jacobian of the part of the manipulator doing the obstacle avoidance task
	void ObstacleKinematics( MatrixD &CriticalPoint, int CriticalPointNumber, const MatrixD &JointPos ); // computes the location of the critical point (point avoiding the obstacle) on the robot
	void GetObstacleAvoidance( MatrixD &ObstacleJacobian, double &ObstacleAvoidanceVel, double &ObstacleSmoothingFactor, const MatrixD &JointPos, const MatrixD &obstacles );
	double GetObstacleDistance( MatrixD &ObstacleUnitVector, int CriticalPointNumber, const MatrixD &ObstacleLocation, const MatrixD &JointPos );

//	Calculates the transformation matrix from the last link frame to the end effector frame
	MatrixD EEFrame( const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE );

//	Allows to change the location and orientation of the end effector frame relative to the last link frame
	void SetNewTool( const Eigen::Vector3d &New_angles_EE, const Eigen::Vector3d &New_r_EE );

//	Eigen::Matrix3d Skew( const Eigen::Vector3d &Vector );

protected:
	int DOF;
	MatrixD MDHangles;
	MatrixD MDHlengthsX;
	MatrixD MDHlengthsZ;
	int *ObstacleAvoidanceDOF;
	double ObstacleAvoidanceNominalVel;
	double ObstacleSphereOfInfluence;
	double ObstacleUnityGain;
	MatrixD ObstacleAvoidancePoints;
	MatrixD A_EE;
};

} //namespace PQPIKS

#endif /* KINEMATICMODELMDH_H_ */
