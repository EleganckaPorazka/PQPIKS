#include "PQPIKS/Trajectory.h"

namespace PQPIKS
{
Trajectory::Trajectory( TrajectoryData trajectoryData )
: Points(trajectoryData.Points), SegmentTimes(trajectoryData.SegmentTimes), SegmentNumber(1), AccelerationPeriod(trajectoryData.AccelerationPeriod), DecelerationPeriod(trajectoryData.DecelerationPeriod), BlendingPeriod(trajectoryData.BlendingPeriod)
{
	//~ TODO: if the number of columns of Points is not equal to the number of elements in SegmentTimes, then FAIL
	if ( (SegmentTimes.rows() != 1 ) and (SegmentTimes.cols() != 1) )
	{
		std::cout << "Wrong size of SegmentTimes. It should be a row or a column vector.\n";
	}
	if ( (Points.rows() != 3) or (Points.cols() != SegmentTimes.size() ) )
	{
		std::cout << "Wrong size of Points. It should be 3 x " << SegmentTimes.size() << ".\n";
	}
	
	PathTimes = MatrixD::Zero(1, SegmentTimes.size());
	if (BlendingPeriod > 0.0)
	{
		PathTimes(1) = SegmentTimes(1);
		for (int i = 2; i < PathTimes.cols(); ++i)
		{
			PathTimes(i) = PathTimes(i - 1) + (1.0 - BlendingPeriod) * SegmentTimes(i);
		}
	}
	else
	{
		for (int i = 1; i < PathTimes.cols(); ++i)
		{
			PathTimes(i) = PathTimes(i - 1) + SegmentTimes(i);
		}
	}
	
	totalTime = PathTimes(PathTimes.cols() - 1);
}

void Trajectory::generator(MatrixD &Pos, MatrixD &Vel, double t)
{
	//~ Get the end effector position and velocity
	
	if (BlendingPeriod > 0.0)
	{
		GetBlendedTrajectory(Pos, Vel, t);
	}
	else
	{
		GetStandardTrajectory(Pos, Vel, t);
	}
}

void Trajectory::GetStandardTrajectory(MatrixD &Pos, MatrixD &Vel, double t)
{
	//~ The robot will stop at the end of the each linear segment.
	
	//~ Check whether it is time to increase the segment number (but only if it is not already the last segment)
	SegmentNumber = getSegment(t);
	double PathTime = PathTimes(SegmentNumber);
	
	//~ Get the start and end positions for the current segment; get the direction unit vector and the segment length; get the times
	Eigen::Vector3d StartPos = Points.col(SegmentNumber - 1);
	Eigen::Vector3d EndPos = Points.col(SegmentNumber);
	double SegmentLength = (EndPos - StartPos).norm();
	Eigen::Vector3d Direction = (EndPos - StartPos) / SegmentLength;
	double SegmentStartTime = PathTime - SegmentTimes(SegmentNumber);
	double SegmentEndTime = SegmentTimes(SegmentNumber);
	double LocalTime = t - SegmentStartTime;
	
	//~ Compute the end effector position and velocity for a given time
	double PathDisp, PathVel;
	//~ ParabolicVelProfile(PathDisp, PathVel, SegmentLength, LocalTime, SegmentEndTime);
	SinusoidalVelProfile(PathDisp, PathVel, SegmentLength, LocalTime, SegmentEndTime);

	Pos = StartPos + Direction * PathDisp;
	Vel = Direction * PathVel;
}

void Trajectory::GetBlendedTrajectory(MatrixD &Pos, MatrixD &Vel, double t)
{
	//~ The robot will not stop at the end of the each linear segment -- the segments will overlap each other.
	
	//~ Check whether it is time to increase the segment number (but only if it is not already the last segment)
	SegmentNumber = getSegment(t);
	double PathTime = PathTimes(SegmentNumber);
	
	//~ Get the start and end positions for the current segment; get the direction unit vector and the segment length; get the times
	Eigen::Vector3d StartPos = Points.col(SegmentNumber - 1);
	Eigen::Vector3d EndPos = Points.col(SegmentNumber);
	double SegmentLength = (EndPos - StartPos).norm();
	Eigen::Vector3d Direction = (EndPos - StartPos) / SegmentLength;
	double SegmentStartTime = PathTime - SegmentTimes(SegmentNumber);
	double SegmentEndTime = SegmentTimes(SegmentNumber);
	double LocalTime = t - SegmentStartTime;
	
	//~ Compute the end effector position and velocity for a given time
	double PathDisp, PathVel;
	SinusoidalVelProfile(PathDisp, PathVel, SegmentLength, LocalTime, SegmentEndTime);
	Pos = StartPos + Direction * PathDisp;
	Vel = Direction * PathVel;
	
	if (SegmentNumber < SegmentTimes.size() - 1) // there is no segment after the last to blend with
	{
		if ( LocalTime > SegmentEndTime - BlendingPeriod * SegmentTimes(SegmentNumber + 1) )
		{
			//~ Time to start the next segment
			int NextSegmentNumber = SegmentNumber + 1;
			
			Eigen::Vector3d NextStartPos = Points.col(NextSegmentNumber - 1);
			Eigen::Vector3d NextEndPos = Points.col(NextSegmentNumber);
			double NextSegmentLength = (NextEndPos - NextStartPos).norm();
			Eigen::Vector3d NextDirection = (NextEndPos - NextStartPos) / NextSegmentLength;
			//~ double NextSegmentStartTime = PathTime;
			double NextSegmentEndTime = SegmentTimes(NextSegmentNumber);
			double NextLocalTime = LocalTime - ( SegmentEndTime - BlendingPeriod * NextSegmentEndTime );
			
			//~ Compute the end effector position and velocity for a given time
			double NextPathDisp, NextPathVel;
			SinusoidalVelProfile(NextPathDisp, NextPathVel, NextSegmentLength, NextLocalTime, NextSegmentEndTime);
			Pos += NextDirection * NextPathDisp;
			Vel += NextDirection * NextPathVel;
		}
	
	}
}

void Trajectory::ParabolicVelProfile( double &PathDisp, double &PathVel, double SegmentLength, double LocalTime, double SegmentEndTime )
{
	//~ Get the displacement along the path and the path speed when the velocity profile is parabolic
	
	PathDisp = SegmentLength * ( 3.0 * (LocalTime / SegmentEndTime) * (LocalTime / SegmentEndTime) - 2.0 * (LocalTime / SegmentEndTime) * (LocalTime / SegmentEndTime) * (LocalTime / SegmentEndTime) );
	
	PathVel = 6.0 * SegmentLength * ( LocalTime / (SegmentEndTime * SegmentEndTime) - (LocalTime * LocalTime) / (SegmentEndTime * SegmentEndTime * SegmentEndTime) );
}

void Trajectory::SinusoidalVelProfile( double &PathDisp, double &PathVel, double SegmentLength, double LocalTime, double SegmentEndTime )
{
	//~ Get the displacement along the path and the path speed when the velocity profile is "sinusoidal"

	double MaxVel = 2.0 * SegmentLength / ( (2.0 - AccelerationPeriod - DecelerationPeriod) * SegmentEndTime );
	
	if (LocalTime < 0.0)
	{
		PathDisp = 0.0;
		PathVel = 0.0;
	}
	else if (LocalTime > SegmentEndTime)
	{
		PathDisp = SegmentLength;
		PathVel = 0.0;
	}
	else if (LocalTime < AccelerationPeriod * SegmentEndTime)
	{
		double x = M_PI / (AccelerationPeriod * SegmentEndTime) * LocalTime;
		PathDisp = MaxVel / ( 4.0 * M_PI * M_PI / (AccelerationPeriod * SegmentEndTime) ) * ( 2.0 * x * x + cos(2.0 * x) ) - (MaxVel * AccelerationPeriod * SegmentEndTime) / (4.0 * M_PI * M_PI);
		PathVel = (x - sin(x) * cos(x)) / M_PI * MaxVel;
	}
	else if (LocalTime > (1.0 - DecelerationPeriod) * SegmentEndTime)
	{
		double PathSoFar = MaxVel / ( 4.0 * M_PI * M_PI / (AccelerationPeriod * SegmentEndTime) ) * ( 2.0 * M_PI * M_PI + 1.0 ) - (MaxVel * AccelerationPeriod * SegmentEndTime) / (4.0 * M_PI * M_PI) + MaxVel * SegmentEndTime * ( 1.0 - (AccelerationPeriod + DecelerationPeriod) );
		double x = -M_PI / (DecelerationPeriod * SegmentEndTime) * ( LocalTime - SegmentEndTime * (1.0 - DecelerationPeriod) );
		double a1 = -M_PI / (DecelerationPeriod * SegmentEndTime);
		double a2 = SegmentEndTime * (1.0 - DecelerationPeriod);
		PathDisp = MaxVel * ( 2.0 * a1 * LocalTime * (a1 * (LocalTime - 2.0 * a2) / M_PI + 2.0) + cos( 2.0 * a1 * (a2 - LocalTime) ) / M_PI ) / (4.0 * a1) - ( MaxVel * ( 2.0 * a1 * a2 * (2.0 - a1 * a2 / M_PI) + (1 / M_PI) ) ) / (4.0 * a1) + PathSoFar;
		PathVel = ( (x - sin(x) * cos(x)) / M_PI + 1.0 ) * MaxVel;
	}
	else
	{
		double PathSoFar = MaxVel / ( 4.0 * M_PI * M_PI / (AccelerationPeriod * SegmentEndTime) ) * ( 2.0 * M_PI * M_PI + 1.0 ) - (MaxVel * AccelerationPeriod * SegmentEndTime) / (4.0 * M_PI * M_PI);
		PathDisp = MaxVel * (LocalTime - AccelerationPeriod * SegmentEndTime) + PathSoFar;
		PathVel = MaxVel;
	}
}

double Trajectory::getNominalMotionTime()
{
	return totalTime;
}

MatrixD Trajectory::getPathTimes()
{
	return PathTimes;
}

int Trajectory::getSegment(double t)
{
	//~ binary search
	int first = 0, last = PathTimes.size() - 1;
	while (first != last)
	{
		int mid = (first + last) / 2;
		if (PathTimes(mid) <= t)
		{
			first = mid + 1;
		}
		else
		{
			last = mid;
		}
	}
	
	if (last == 0)
	{
		last = 1;
	}
	
	return last;
}
	
} //namespace PQPIKS
