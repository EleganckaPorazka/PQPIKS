#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <iostream>
#include <Eigen/Dense>
#include "typedefs.h"

namespace PQPIKS
{

struct TrajectoryData
{
	double totalTime{0.0};	// nominal motion time before scaling
	MatrixD Points;			// points defining the trajectory
	MatrixD SegmentTimes;
	MatrixD q0;				// joint space starting point
	double AccelerationPeriod{0.1}, DecelerationPeriod{0.1}, BlendingPeriod{0.0};
};

class Trajectory
{
	public:
		Trajectory( TrajectoryData trajectoryData );
		
		//~ Get the end effector position and velocity
		void generator(MatrixD &Pos, MatrixD &Vel, double t);
		
		void GetStandardTrajectory(MatrixD &Pos, MatrixD &Vel, double t);
		
		void GetBlendedTrajectory(MatrixD &Pos, MatrixD &Vel, double t);
		
		//~ Get the displacement along the path and the path speed when the velocity profile is parabolic
		void ParabolicVelProfile( double &PathDisp, double &PathVel, double SegmentLength, double LocalTime, double SegmentEndTime );
		
		//~ Get the displacement along the path and the path speed when the velocity profile is "sinusoidal"
		void SinusoidalVelProfile( double &PathDisp, double &PathVel, double SegmentLength, double LocalTime, double SegmentEndTime );
		
		double getNominalMotionTime();
		
		MatrixD getPathTimes();
		
		int getSegment(double t);
		
	protected:
		MatrixD Points;				// 3*NumberOfPoints matrix; (x,y,z) coordinates of the points making up the path (linear segments)
		MatrixD SegmentTimes;		// NumberOfPoints row vector with times needed to complete each segment; the index corresponds to the point being the end of the segment (index - 1 corresponds to the point starting the segment)
		int SegmentNumber;			// number of the current trajectory segment
		MatrixD PathTimes;			// time needed to complete the path from the start of the first segment to the end of the current segment
		double totalTime{0.0};		// total time of the trajectory
		double AccelerationPeriod;	// part of the segment time allocated for acceleration
		double DecelerationPeriod;	// part of the segment time allocated for deceleration
		double BlendingPeriod;		// the fraction of the segment's time to blend two linear segments of the trajectory
};

} //namespace PQPIKS

#endif /* TRAJECTORY_H_ */
