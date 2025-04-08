// Based on the class PQPIKS provides an implementiation using modified Denavit-Hartenberg (MDH) parameters

#ifndef PQPIKS_MDH_H_
#define PQPIKS_MDH_H_

#include "PQPIKS.h"
#include "KinematicModelMDH.h"
#include "DynamicModelMDH.h"
//~ #include "SimpleTrajectory.h"
#include "Trajectory.h"

const double pi = M_PI;

namespace PQPIKS
{

class PQPIKS_MDH : public PQPIKS
{
	public:
		PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData ); // constructor with default kinematic and dynamic data (LWR 4+)
		PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE );	// parameterized constructor with only kinematic data (default dynamic data LWR 4+)
		PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity ); // constructor without an EE
		PQPIKS_MDH( PQPIKSdata inputData, TrajectoryData trajectoryData, const MatrixD &MDHangles, const MatrixD &MDHlengthsX, const MatrixD &MDHlengthsZ, const MatrixD &AllRotAxes, const MatrixD &Masses, const MatrixD &AllFirstMomentsOfMass, const MatrixD &AllInertias, const MatrixD &Gravity, const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE ); // full parameterized constructor
		
		virtual ~PQPIKS_MDH();
		
		void initTrajectory(TrajectoryData &trajectoryData);
		
		virtual void jacobian( MatrixD &J, const MatrixD &q );
		virtual void forwardKin( MatrixD &T, const MatrixD &q );
		virtual void errorEE( MatrixD &e, const MatrixD &x, const MatrixD &q );
		virtual void trajectoryGenerator( MatrixD &x, MatrixD &dxdt, double t );
		virtual void inertiaMatrix( MatrixD &M, const MatrixD &q );
		virtual void obstacleJacobianAndVelocity( MatrixD &Jo, double &dxdto, double &woSmoothing, const MatrixD &q, const MatrixD &obstacles );
		virtual double obstacleDistance( const MatrixD &q, const MatrixD &obstacles );
		
		double getNominalMotionTime();
		
		void setNewTool( const Eigen::Vector3d &angles_EE, const Eigen::Vector3d &r_EE );
		
	protected:
		KinematicModelMDH *Kin;
		DynamicModelMDH *Dyn;
		Trajectory *Traj;
};

} /* namespace PQPIKS */

#endif /* PQPIKS_MDH_H_ */
