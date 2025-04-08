#include "PQPIKS/KinematicModelMDH.h"

namespace PQPIKS
{

KinematicModelMDH::KinematicModelMDH( int DOF, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, int *ObstacleAvoidanceDOF, double ObstacleAvoidanceNominalVel, double ObstacleSphereOfInfluence, double ObstacleUnityGain, const MatrixD &ObstacleAvoidancePoints )
: DOF(DOF), MDHangles(MDHangles), MDHlengthsX(MDHlengthsX), MDHlengthsZ(MDHlengthsZ), ObstacleAvoidanceDOF(ObstacleAvoidanceDOF), ObstacleAvoidanceNominalVel(ObstacleAvoidanceNominalVel), ObstacleSphereOfInfluence(ObstacleSphereOfInfluence), ObstacleUnityGain(ObstacleUnityGain), ObstacleAvoidancePoints(ObstacleAvoidancePoints)
{
	A_EE = MatrixD::Identity(4,4);
}

KinematicModelMDH::KinematicModelMDH( int DOF, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE, int *ObstacleAvoidanceDOF, double ObstacleAvoidanceNominalVel, double ObstacleSphereOfInfluence, double ObstacleUnityGain, const MatrixD &ObstacleAvoidancePoints )
: DOF(DOF), MDHangles(MDHangles), MDHlengthsX(MDHlengthsX), MDHlengthsZ(MDHlengthsZ), ObstacleAvoidanceDOF(ObstacleAvoidanceDOF), ObstacleAvoidanceNominalVel(ObstacleAvoidanceNominalVel), ObstacleSphereOfInfluence(ObstacleSphereOfInfluence), ObstacleUnityGain(ObstacleUnityGain), ObstacleAvoidancePoints(ObstacleAvoidancePoints)
{
	A_EE = EEFrame(angles_EE, r_EE);
}

KinematicModelMDH::KinematicModelMDH( const KinematicModelMDH &other )
: DOF(other.DOF), MDHangles(other.MDHangles), MDHlengthsX(other.MDHlengthsX), MDHlengthsZ(other.MDHlengthsZ), ObstacleAvoidanceDOF(other.ObstacleAvoidanceDOF), ObstacleAvoidanceNominalVel(other.ObstacleAvoidanceNominalVel), ObstacleSphereOfInfluence(other.ObstacleSphereOfInfluence), ObstacleUnityGain(other.ObstacleUnityGain), ObstacleAvoidancePoints(other.ObstacleAvoidancePoints), A_EE(other.A_EE)
{
}

KinematicModelMDH& KinematicModelMDH::operator= ( const KinematicModelMDH &other )
{
	if (&other == this) // self-assignment
	{
		return *this;
	}
	
	DOF = other.DOF;	
	MDHangles = other.MDHangles;
	MDHlengthsX = other.MDHlengthsX;
	MDHlengthsZ = other.MDHlengthsZ;
	ObstacleAvoidanceDOF = other.ObstacleAvoidanceDOF;
	ObstacleAvoidanceNominalVel = other.ObstacleAvoidanceNominalVel;
	ObstacleSphereOfInfluence = other.ObstacleSphereOfInfluence;
	ObstacleUnityGain = other.ObstacleUnityGain;
	ObstacleAvoidancePoints = other.ObstacleAvoidancePoints;
	A_EE = other.A_EE;
	
	return *this;
}

KinematicModelMDH::~KinematicModelMDH()
{
}

void KinematicModelMDH::ForwardKinematics( MatrixD &T, const MatrixD &JointPos )
{
//	Solves the forward kinematics problem as a transformation matrix from the base frame to the end effector frame.
//	A_j after the loop is the transformation matrix from the base frame to the last link frame,
//	A_EE is the transformation matrix from the last link to the end effector frame,
//	T = A_j * A_EE is the transformation matrix from the base frame to the end effector frame.

	MatrixD A_i = MatrixD::Identity(4, 4);

	MatrixD A_ij(4, 4), A_j(4, 4);

	for (int j = 0; j < DOF; j++)
	{
		A_ij << cos( JointPos(j) ), -sin( JointPos(j) ), 0.0, MDHlengthsX(j),
				sin( JointPos(j) ) * cos( MDHangles(j) ), cos( JointPos(j) ) * cos( MDHangles(j) ), -sin( MDHangles(j) ), -MDHlengthsZ(j) * sin( MDHangles(j) ),
				sin( JointPos(j) ) * sin( MDHangles(j) ), cos( JointPos(j) ) * sin( MDHangles(j) ), cos( MDHangles(j) ), MDHlengthsZ(j) * cos( MDHangles(j) ),
				0.0, 0.0, 0.0, 1.0;

		A_j = A_i * A_ij;
		A_i = A_j;
	}

	T = A_j * A_EE;
}

void KinematicModelMDH::GetJacobian( MatrixD &J, const MatrixD &JointPos )
{
//	Computes the manipulator's Jacobian J (size 3 * DOF).

	MatrixD A_i = MatrixD::Identity(4, 4);

	MatrixD A_ij(4, 4), A_j(4, 4);

	MatrixD z(3, DOF), p(3, DOF);

	for (int j = 0; j < DOF; j++)
	{
		p.col(j) = A_i.block<3,1>(0,3);

		A_ij << cos( JointPos(j) ), -sin( JointPos(j) ), 0.0, MDHlengthsX(j),
				sin( JointPos(j) ) * cos( MDHangles(j) ), cos( JointPos(j) ) * cos( MDHangles(j) ), -sin( MDHangles(j) ), -MDHlengthsZ(j) * sin( MDHangles(j) ),
				sin( JointPos(j) ) * sin( MDHangles(j) ), cos( JointPos(j) ) * sin( MDHangles(j) ), cos( MDHangles(j) ), MDHlengthsZ(j) * cos( MDHangles(j) ),
				0.0, 0.0, 0.0, 1.0;

		A_j = A_i * A_ij;
		A_i = A_j;
		MatrixD one(3, 1);
		one << 0.0, 0.0, 1.0;
		z.col(j) = A_i.block<3,3>(0,0) * one;
		//~ z.col(j) = A_i.block<3, 1>(0, 2);
	}

	MatrixD p_EE(3, 1);
	p_EE = (A_j * A_EE).block<3,1>(0,3);

	for (int j = 0; j < DOF; j++)
	{
		Eigen::Vector3d pTemp = p_EE - p.col(j);
		J.col(j) = Eigen::Vector3d(z.col(j)).cross( pTemp );
	}
}

void KinematicModelMDH::GetObstacleJacobian( MatrixD &Jo, int CriticalPointNumber, const MatrixD &JointPos )
{
//	Computes the Jacobian of the part of the manipulator doing the obstacle avoidance task (size 3 * DOF).

	MatrixD A_i = MatrixD::Identity(4, 4);

	MatrixD A_ij(4, 4), A_j(4, 4);
	
	MatrixD z(3, DOF), p(3, DOF);

	// we will use the number of joints equal to the ObstacleAvoidanceDOF to avoid the obstacle
	for (int j = 0; j < ObstacleAvoidanceDOF[CriticalPointNumber]; j++)
	{
		p.col(j) = A_i.block<3,1>(0,3);

		A_ij << cos( JointPos(j) ), -sin( JointPos(j) ), 0.0, MDHlengthsX(j),
				sin( JointPos(j) ) * cos( MDHangles(j) ), cos( JointPos(j) ) * cos( MDHangles(j) ), -sin( MDHangles(j) ), -MDHlengthsZ(j) * sin( MDHangles(j) ),
				sin( JointPos(j) ) * sin( MDHangles(j) ), cos( JointPos(j) ) * sin( MDHangles(j) ), cos( MDHangles(j) ), MDHlengthsZ(j) * cos( MDHangles(j) ),
				0.0, 0.0, 0.0, 1.0;

		A_j = A_i * A_ij;
		A_i = A_j;
		MatrixD one(3, 1);
		one << 0.0, 0.0, 1.0;
		z.col(j) = A_i.block<3,3>(0,0) * one;
		//~ z.col(j) = A_i.block<3, 1>(0, 2);
	}
	
	//~ The vector from the beginning of the ObstacleAvoidanceDOF-th link frame locating the critical point
	MatrixD A_ObstacleAvoidance = MatrixD::Identity(4, 4);
	//~ MatrixD ObstacleAvoidancePoint = ObstacleAvoidancePoints.col(0); //(0.0, MDHlengthsZ(2), 0.0);
	A_ObstacleAvoidance.block<3,1>(0,3) = ObstacleAvoidancePoints.col(CriticalPointNumber);
	
	MatrixD p_EE = (A_j * A_ObstacleAvoidance).block<3,1>(0,3);

	for (int j = 0; j < ObstacleAvoidanceDOF[CriticalPointNumber]; j++)
	{
		Eigen::Vector3d pTemp = p_EE - p.col(j);
		Jo.col(j) = Eigen::Vector3d(z.col(j)).cross( pTemp );
	}
}

void KinematicModelMDH::ObstacleKinematics( MatrixD &CriticalPoint, int CriticalPointNumber, const MatrixD &JointPos )
{
//~ Computes the location of the critical point (point avoiding the obstacle) on the robot

	MatrixD A_i = MatrixD::Identity(4, 4);

	MatrixD A_ij(4, 4), A_j(4, 4);
	
	// we will use the number of joints equal to the ObstacleAvoidanceDOF to avoid the obstacle
	for (int j = 0; j < ObstacleAvoidanceDOF[CriticalPointNumber]; j++)
	{
		A_ij << cos( JointPos(j) ), -sin( JointPos(j) ), 0.0, MDHlengthsX(j),
				sin( JointPos(j) ) * cos( MDHangles(j) ), cos( JointPos(j) ) * cos( MDHangles(j) ), -sin( MDHangles(j) ), -MDHlengthsZ(j) * sin( MDHangles(j) ),
				sin( JointPos(j) ) * sin( MDHangles(j) ), cos( JointPos(j) ) * sin( MDHangles(j) ), cos( MDHangles(j) ), MDHlengthsZ(j) * cos( MDHangles(j) ),
				0.0, 0.0, 0.0, 1.0;

		A_j = A_i * A_ij;
		A_i = A_j;
	}

	//~ The vector from the beginning of the ObstacleAvoidanceDOF-th link frame locating the critical point
	MatrixD A_ObstacleAvoidance = MatrixD::Identity(4, 4);
	A_ObstacleAvoidance.block<3,1>(0,3) = ObstacleAvoidancePoints.col(CriticalPointNumber);
	
	CriticalPoint = (A_j * A_ObstacleAvoidance).block<3,1>(0,3);
}

void KinematicModelMDH::GetObstacleAvoidance( MatrixD &ObstacleJacobian, double &ObstacleAvoidanceVel, double &ObstacleSmoothingFactor, const MatrixD &JointPos, const MatrixD &obstacles )
{
	//~ T. Petrič et al. “Obstacle Avoidance with Industrial Robots”. In: Motion and Operation Planning of Robotic Systems: Background and Practical Approaches. Ed. by Giuseppe Carbone and Fernando Gomez-Bravo. Cham: Springer International Publishing, 2015, pp. 113–145. isbn: 978-3-319-14705-5. doi: 10.1007/978-3-319-14705-5_5 . url: https://doi.org/10.1007/978-3-319-14705-5_5
	
	//~ Anthony A. Maciejewski and Charles A. Klein. “Obstacle Avoidance for Kinematically Redundant Manipulators in Dynamically Varying Environments”. In: The International Journal of Robotics Research 4.3 (1985), pp. 109–117. doi: 10.1177/027836498500400308 . eprint: https://doi.org/10.1177/027836498500400308 . url: https://doi.org/10.1177/027836498500400308 .
	
	int CriticalPointNumber = 0;

	//~ Get the location of the obstacle
	MatrixD ObstacleLocation = obstacles.col(CriticalPointNumber);
	
	//~ Compute the unit vector in the direction between the critical point and the obstacle
	MatrixD ObstacleUnitVector(3, 1);
	double ObstacleDistance = GetObstacleDistance( ObstacleUnitVector, CriticalPointNumber, ObstacleLocation, JointPos );
	
	if (ObstacleDistance < ObstacleSphereOfInfluence)
	{
		//~ Compute the ObstacleSmoothingFactor
		if (ObstacleDistance <= ObstacleUnityGain)
		{
			ObstacleSmoothingFactor = 1.0;
		}
		else
		{
			ObstacleSmoothingFactor = 0.5 * ( 1.0 - cos( M_PI * (ObstacleDistance - ObstacleSphereOfInfluence) / (ObstacleUnityGain - ObstacleSphereOfInfluence) ) );
		}
		
		//~ Compute the ObstacleAvoidanceVel
		double ObstacleAvoidanceVelGain;
		if (ObstacleDistance < ObstacleUnityGain)
		{
			ObstacleAvoidanceVelGain = (ObstacleUnityGain / ObstacleDistance) * (ObstacleUnityGain / ObstacleDistance) - 1.0;
		}
		else
		{
			ObstacleAvoidanceVelGain = 0.0;
		}
		ObstacleAvoidanceVel = ObstacleAvoidanceVelGain * ObstacleAvoidanceNominalVel;
		
		//~ Compute the obstacle avoidance task Jacobian
		MatrixD Jo = MatrixD::Zero(3, DOF);
		GetObstacleJacobian(Jo, CriticalPointNumber, JointPos);
		ObstacleJacobian = ObstacleUnitVector.transpose() * Jo;
	}
	else
	{
		//~ If outside the SphereOfInfluence:
		ObstacleJacobian = MatrixD::Zero(1, DOF);
		ObstacleAvoidanceVel = 0.0;
		ObstacleSmoothingFactor = 0.0;
	}
}

double KinematicModelMDH::GetObstacleDistance( MatrixD &ObstacleUnitVector, int CriticalPointNumber, const MatrixD &ObstacleLocation, const MatrixD &JointPos )
{
	//~ Get the location of the critical point on the robot
	MatrixD CriticalPoint(3, 1);
	ObstacleKinematics( CriticalPoint, CriticalPointNumber, JointPos );
	
	//~ Compute the vector connecting the critical point on the robot with the point on the obstacle
	MatrixD ObstacleVector = CriticalPoint - ObstacleLocation;
	
	//~ Compute the unit vector in the direction between the critical point and the obstacle
	double ObstacleDistance = ObstacleVector.norm();
	ObstacleUnitVector = ObstacleVector / ObstacleDistance;
	return ObstacleDistance;
}

MatrixD KinematicModelMDH::EEFrame( const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE )
{
//	Computes the transformation matrix from the last link frame to the end effector frame.
	MatrixD A_EE(4, 4);

	A_EE << cos( angles_EE(2) ) * cos( angles_EE(0) ) - sin( angles_EE(2) ) * cos( angles_EE(1) ) * sin( angles_EE(0) ), -sin( angles_EE(2) ) * cos( angles_EE(0) ) - cos( angles_EE(2) ) * cos( angles_EE(1) ) * sin( angles_EE(0) ), sin( angles_EE(1) ) * sin( angles_EE(0) ), r_EE(0),
			sin( angles_EE(2) ) * cos( angles_EE(1) ) * cos( angles_EE(0) ) + cos( angles_EE(2) ) * sin( angles_EE(0) ), cos( angles_EE(2) ) * cos( angles_EE(1) ) * cos( angles_EE(0) ) - sin( angles_EE(2) ) * sin( angles_EE(0) ), -sin( angles_EE(1) ) * cos( angles_EE(0) ), r_EE(1),
			sin( angles_EE(2) ) * sin( angles_EE(1) ), cos( angles_EE(2) ) * sin( angles_EE(1) ), cos( angles_EE(1) ), r_EE(2),
			0.0, 0.0, 0.0, 1.0;

	return A_EE;
}

void KinematicModelMDH::SetNewTool( const Eigen::Vector3d &New_angles_EE, const Eigen::Vector3d &New_r_EE )
{
//	Allows to change the location and orientation of the end effector frame relative to the last link frame.
	A_EE = EEFrame( New_angles_EE, New_r_EE );
}

/*Eigen::Matrix3d Skew( const Eigen::Vector3d &Vector )
{
//	Returns a skew-symmetric matrix from a vector.
	Eigen::Matrix3d SkewSymmetricMatrix;

	SkewSymmetricMatrix << 0.0, -Vector(2), Vector(1),
			Vector(2), 0.0, -Vector(0),
			-Vector(1), Vector(0), 0.0;

	return SkewSymmetricMatrix;
}*/

} //namespace PQPIKS
