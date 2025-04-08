#include "PQPIKS/PQPIKS_MDH.h"

namespace PQPIKS
{

PQPIKS_MDH::PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData )
: PQPIKS(inputData)
{
	int n = 7;
	//	Modified Denavit-Hartenberg parameters:
	MatrixD MDHangles(n, 1); // angles
	MDHangles << 0.0, pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2;
	MatrixD MDHlengthsX(n, 1); // lengths on X axis
	MDHlengthsX << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	MatrixD MDHlengthsZ(n, 1); // lengths on Z axis
	MDHlengthsZ << 0.31, 0.0, 0.4, 0.0, 0.39, 0.0, 0.0;

	Kin = new KinematicModelMDH( n, MDHangles, MDHlengthsX, MDHlengthsZ, inputData.obstacleAvoidanceDOF, inputData.obstacleAvoidanceNominalVel, inputData.obstacleSphereOfInfluence, inputData.obstacleUnityGain, inputData.obstacleAvoidancePoints );
	
	//~ Inertial parameters (source: A. Jubien, M. Gautier, A. Janot. Dynamic identification of the KUKA lightweight robot: Comparison between actual and confidential KUKA’s parameters. In: 2014 IEEE/ASME International Conference on Advanced Intelligent Mechatronics. Proceedings, July, 2014, pp. 483–488)
	MatrixD AllRotAxes = MatrixD::Zero(3, n);
	AllRotAxes.row(2) = MatrixD::Ones(1, n);
	
	MatrixD Masses = MatrixD::Zero(n, 1);
	
	MatrixD AllFirstMomentsOfMass(3, n);
	AllFirstMomentsOfMass.col(0) = MatrixD::Zero(3, 1);
	AllFirstMomentsOfMass.col(1) = Eigen::Vector3d( 1.35 * pow(10.0, -3.0), 3.46, 0.0 );
	AllFirstMomentsOfMass.col(2) = Eigen::Vector3d( 9.45 * pow(10.0, -4.0), -4.72 * pow(10.0, -4.0), 0.0 );
	AllFirstMomentsOfMass.col(3) = Eigen::Vector3d( -3.5 * pow(10.0, -3.0), -1.33, 0.0 );
	AllFirstMomentsOfMass.col(4) = Eigen::Vector3d( -6.56 * pow(10.0, -4.0), 4.07 * pow(10.0, -2.0), 0.0 );
	AllFirstMomentsOfMass.col(5) = Eigen::Vector3d( 8.35 * pow(10.0, -4.0), 2.86 * pow(10.0, -2.0), 0.0 );
	AllFirstMomentsOfMass.col(6) = Eigen::Vector3d( -2.98 * pow(10.0, -4.0), 9.54 * pow(10.0, -4.0), 0.0 );
	
	MatrixD AllInertias(3, 3 * n);
	AllInertias.block<3,3>(0,0) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15 * pow(10.0, -2.0);
	AllInertias.block<3,3>(0,3) << 1.25, 0.0, 0.0, 0.0, 0.0, -5.43 * pow(10.0, -4.0), 0.0, 5.43 * pow(10.0, -4.0), 1.25;
	AllInertias.block<3,3>(0,6) << 6.36 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 7.26 * pow(10.0, -4.0), 0.0, -7.26 * pow(10.0, -4.0), 1.08 * pow(10.0, -2.0);
	AllInertias.block<3,3>(0,9) << 0.413, 0.0, 0.0, 0.0, 0.0, 5.32 * pow(10.0, -4.0), 0.0, -5.32 * pow(10.0, -4.0), 0.418;
	AllInertias.block<3,3>(0,12) << 3.95 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 4.22 * pow(10.0, -4.0), 0.0, -4.22 * pow(10.0, -4.0), 6.33 * pow(10.0, -3.0);
	AllInertias.block<3,3>(0,15) << 1.17 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 1.02 * pow(10.0, -5.0), 0.0, -1.02 * pow(10.0, -5.0), 3.77 * pow(10.0, -3.0);
	AllInertias.block<3,3>(0,18) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2 * pow(10.0, -4.0);
	
	MatrixD Gravity(3, 1);
	Gravity << 0.0, 0.0, -9.80665;
	
	Dyn = new DynamicModelMDH( n, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity, MDHangles, MDHlengthsX, MDHlengthsZ );
	
	initTrajectory(trajectoryData);
	Traj = new Trajectory(trajectoryData);
}

PQPIKS_MDH::PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE )
: PQPIKS(inputData)
{
	int n = inputData.n;
	
	// TODO: add some error handling: if the sizes of input matrices are wrong, then fail gracefully
	
	Kin = new KinematicModelMDH( n, MDHangles, MDHlengthsX, MDHlengthsZ, angles_EE, r_EE, inputData.obstacleAvoidanceDOF, inputData.obstacleAvoidanceNominalVel, inputData.obstacleSphereOfInfluence, inputData.obstacleUnityGain, inputData.obstacleAvoidancePoints );
	
	//~ Inertial parameters (source: A. Jubien, M. Gautier, A. Janot. Dynamic identification of the KUKA lightweight robot: Comparison between actual and confidential KUKA’s parameters. In: 2014 IEEE/ASME International Conference on Advanced Intelligent Mechatronics. Proceedings, July, 2014, pp. 483–488)
	MatrixD AllRotAxes = MatrixD::Zero(3, n);
	AllRotAxes.row(2) = MatrixD::Ones(1, n);
	
	MatrixD Masses = MatrixD::Zero(n, 1);
	
	MatrixD AllFirstMomentsOfMass(3, n);
	AllFirstMomentsOfMass.col(0) = MatrixD::Zero(3, 1);
	AllFirstMomentsOfMass.col(1) = Eigen::Vector3d( 1.35 * pow(10.0, -3.0), 3.46, 0.0 );
	AllFirstMomentsOfMass.col(2) = Eigen::Vector3d( 9.45 * pow(10.0, -4.0), -4.72 * pow(10.0, -4.0), 0.0 );
	AllFirstMomentsOfMass.col(3) = Eigen::Vector3d( -3.5 * pow(10.0, -3.0), -1.33, 0.0 );
	AllFirstMomentsOfMass.col(4) = Eigen::Vector3d( -6.56 * pow(10.0, -4.0), 4.07 * pow(10.0, -2.0), 0.0 );
	AllFirstMomentsOfMass.col(5) = Eigen::Vector3d( 8.35 * pow(10.0, -4.0), 2.86 * pow(10.0, -2.0), 0.0 );
	AllFirstMomentsOfMass.col(6) = Eigen::Vector3d( -2.98 * pow(10.0, -4.0), 9.54 * pow(10.0, -4.0), 0.0 );
	
	MatrixD AllInertias(3, 3 * n);
	AllInertias.block<3,3>(0,0) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15 * pow(10.0, -2.0);
	AllInertias.block<3,3>(0,3) << 1.25, 0.0, 0.0, 0.0, 0.0, -5.43 * pow(10.0, -4.0), 0.0, 5.43 * pow(10.0, -4.0), 1.25;
	AllInertias.block<3,3>(0,6) << 6.36 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 7.26 * pow(10.0, -4.0), 0.0, -7.26 * pow(10.0, -4.0), 1.08 * pow(10.0, -2.0);
	AllInertias.block<3,3>(0,9) << 0.413, 0.0, 0.0, 0.0, 0.0, 5.32 * pow(10.0, -4.0), 0.0, -5.32 * pow(10.0, -4.0), 0.418;
	AllInertias.block<3,3>(0,12) << 3.95 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 4.22 * pow(10.0, -4.0), 0.0, -4.22 * pow(10.0, -4.0), 6.33 * pow(10.0, -3.0);
	AllInertias.block<3,3>(0,15) << 1.17 * pow(10.0, -3.0), 0.0, 0.0, 0.0, 0.0, 1.02 * pow(10.0, -5.0), 0.0, -1.02 * pow(10.0, -5.0), 3.77 * pow(10.0, -3.0);
	AllInertias.block<3,3>(0,18) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2 * pow(10.0, -4.0);
	
	MatrixD Gravity(3, 1);
	Gravity << 0.0, 0.0, -9.80665;
	
	Dyn = new DynamicModelMDH( n, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity, MDHangles, MDHlengthsX, MDHlengthsZ );
	
	initTrajectory(trajectoryData);
	Traj = new Trajectory(trajectoryData);
}

PQPIKS_MDH::PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity )
: PQPIKS(inputData)
{
	int n = inputData.n;
	
	// TODO: add some error handling: if the sizes of input matrices are wrong, then fail gracefully
	
	Kin = new KinematicModelMDH( n, MDHangles, MDHlengthsX, MDHlengthsZ, inputData.obstacleAvoidanceDOF, inputData.obstacleAvoidanceNominalVel, inputData.obstacleSphereOfInfluence, inputData.obstacleUnityGain, inputData.obstacleAvoidancePoints );
	
	Dyn = new DynamicModelMDH( n, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity, MDHangles, MDHlengthsX, MDHlengthsZ );
	
	initTrajectory(trajectoryData);
	Traj = new Trajectory(trajectoryData);
}

PQPIKS_MDH::PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE )
: PQPIKS(inputData)
{
	int n = inputData.n;
	
	// TODO: add some error handling: if the sizes of input matrices are wrong, then fail gracefully
	
	Kin = new KinematicModelMDH( n, MDHangles, MDHlengthsX, MDHlengthsZ, angles_EE, r_EE, inputData.obstacleAvoidanceDOF, inputData.obstacleAvoidanceNominalVel, inputData.obstacleSphereOfInfluence, inputData.obstacleUnityGain, inputData.obstacleAvoidancePoints );
	
	Dyn = new DynamicModelMDH( n, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity, MDHangles, MDHlengthsX, MDHlengthsZ );
	
	initTrajectory(trajectoryData);
	Traj = new Trajectory(trajectoryData);
}

PQPIKS_MDH::~PQPIKS_MDH()
{
}

void PQPIKS_MDH::initTrajectory(TrajectoryData &trajectoryData)
{
	MatrixD T(4, 4);
	forwardKin(T, trajectoryData.q0);
	MatrixD x0 = T.block(0, 3, 3, 1);
	MatrixD points(3, trajectoryData.Points.cols() + 1);
	points.col(0) = x0;
	points.block(0, 1, 3, trajectoryData.Points.cols()) = trajectoryData.Points;
	trajectoryData.Points = points;
}

void PQPIKS_MDH::jacobian( MatrixD &J, const MatrixD &q )
{
	Kin->GetJacobian(J, q);
}

void PQPIKS_MDH::forwardKin( MatrixD &T, const MatrixD &q )
{
	Kin->ForwardKinematics(T, q);
}

void PQPIKS_MDH::errorEE( MatrixD &e, const MatrixD &x, const MatrixD &q )
{
	//	End effector position error for the CLIK
	MatrixD T(4, 4);
	Kin->ForwardKinematics(T, q);
	e = x - T.block<3,1>(0,3);
}

void PQPIKS_MDH::trajectoryGenerator( MatrixD &x, MatrixD &dxdt, double t )
{
	Traj->generator(x, dxdt, t);
}

void PQPIKS_MDH::inertiaMatrix( MatrixD &M, const MatrixD &q )
{
	Dyn->InertiaMatrix(M, q);
}

void PQPIKS_MDH::obstacleJacobianAndVelocity( MatrixD &Jo, double &dxdto, double &woSmoothing, const MatrixD &q, const MatrixD &obstacles )
{
	Kin->GetObstacleAvoidance(Jo, dxdto, woSmoothing, q, obstacles);
}

double PQPIKS_MDH::obstacleDistance( const MatrixD &q, const MatrixD &obstacles )
{
	MatrixD obstacleUnitVector(3, 1);
	int criticalPointNumber = 0;
	return Kin->GetObstacleDistance( obstacleUnitVector, criticalPointNumber, obstacles.col(criticalPointNumber), q );
}

double PQPIKS_MDH::getNominalMotionTime()
{
	return Traj->getNominalMotionTime();
}

void PQPIKS_MDH::setNewTool( const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE )
{
	Kin->SetNewTool( angles_EE, r_EE );
}

} /* namespace PQPIKS */
