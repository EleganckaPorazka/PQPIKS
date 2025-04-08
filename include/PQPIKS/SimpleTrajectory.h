#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "typedefs.h"

namespace PQPIKS
{

struct TrajectoryData
{
	double totalTime;	// nominal motion time before scaling
	MatrixD Points;		// points defining the trajectory
	MatrixD q0;			// joint space starting point
};

class Trajectory
{
/* Defines a polynomial trajectory generator */
	public:
		Trajectory(TrajectoryData inputData);
		void generator(MatrixD &x, MatrixD &dxdt, double t);
		void timingLaw(double &u, double &dudt, double t);
		double getNominalMotionTime();
		
	private:
		double totalTime;	// total time
		MatrixD x0;			// start position
		MatrixD xEnd;		// end position
		MatrixD v;			// vector from x0 to xEnd

};

Trajectory::Trajectory(TrajectoryData inputData)
: totalTime(inputData.totalTime)
{
	if ( (inputData.Points.rows() != 3) or (inputData.Points.cols() != 2) )
	{
		std::cout << "Wrong size of the Points matrix. It should be 3 x 2.\n";
	}
	// TODO: error handling
	
	x0 = inputData.Points.col(0);
	xEnd = inputData.Points.col(1);
	v = xEnd - x0;
}

void Trajectory::generator(MatrixD &x, MatrixD &dxdt, double t)
{
	/* returns position x and velocity dxdt for a given time t */
	double u, dudt;
	timingLaw(u, dudt, t);
	x = x0 + v * u;
	dxdt = v * dudt;
}

void Trajectory::timingLaw(double &u, double &dudt, double t)
{
	if (t < 0.0)
	{
		u = 0.0;
		dudt = 0.0;
	}
	else if (t <= totalTime)
	{		
		u = 6.0 / pow(totalTime, 5.0) * pow(t, 5.0) - 15.0 / pow(totalTime, 4) * pow(t, 4.0) + 10.0 / pow(totalTime, 3.0) * pow(t, 3.0);
		dudt = 30.0 / pow(totalTime, 5.0) * pow(t, 4.0) - 60.0 / pow(totalTime, 4.0) * pow(t, 3.0) + 30.0 / pow(totalTime, 3.0) * pow(t, 2.0);
	}
	else
	{
		u = 1.0;
		dudt = 0.0;
	}
}

double Trajectory::getNominalMotionTime()
{
	return totalTime;
}

} /* namespace PQPIKS */
