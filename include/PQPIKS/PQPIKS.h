#ifndef PQPIKS_H_
#define PQPIKS_H_

#include <iostream>
#include <Eigen/Dense>
#include <qpOASES.hpp>
#include "typedefs.h"

USING_NAMESPACE_QPOASES

namespace PQPIKS
{

struct PQPIKSdata
{
	// number of:
	int n;	// robot DOF
	int m;	// main task DOF
	int p;	// prediction steps
	
	double dt;			// cycle time
	double tHorizon;	// length of the prediction horizon in time
	
	MatrixD K;			// gains for CLIK (for each predicted step p)
	
	MatrixD qLim, dqdtLim, d2qdt2Lim;	// joint limits (qLim.col(0) = qMin, qLim.col(1) = qMax, etc.)
	// IMPORTANT: right now, joint limits have to be symmetric
	
	double dsdtMax{1.0};				// maximum path velocity; dsdtMax = 1 is the nominal velocity, dsdtMax < 1 forces a slowdown
	
	bool minimizeKineticEnergy{false};	// whether to include the minimization of the kinetic energy in the goal function or not
	
	// obstacle avoidance parameters
	int *obstacleAvoidanceDOF{0};							// how many DOF (starting from the robot base) to use for obstacle avoidance task
	double obstacleAvoidanceNominalVel{100.0}; //0.3?
	double obstacleSphereOfInfluence{0.4};
	double obstacleUnityGain{0.35};
	MatrixD obstacleAvoidancePoints{MatrixD::Zero(3, 1)};	// each point should correspond to DOF in obstacleAvoidanceDOF
	
	// qpOASES parameters
	real_t cputime0{1.0};	// for CPU time measurement
	int_t nWSRmax{100};		// maximum number of working set recalculations
};

class PQPIKS
{
	public:
		PQPIKS( PQPIKSdata inputData );
		
		virtual ~PQPIKS();
		
		void stateMatrices();
		void stateMatricesSizeReduction(MatrixD &Ssize, MatrixD &Sblocking, int N);
		void initConstraintsAndGoal();
		
		void constraints(const MatrixD &z0, const MatrixD &QEst, const MatrixD &dsdt, double s);
		
		void goalFunction(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt);
		void goalFunction(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt, double wKinE);
		
		void goalFunctionObstacleAvoidance(const MatrixD &z0, const MatrixD &QEst, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo);
		void goalFunctionObstacleAvoidance(const MatrixD &z0, const MatrixD &QEst, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo, double wKinE);
		
		// keep the results as members of the class and have a function getResults to return them?
		void solve(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, double wU, double wdsdt, double wdqdt, double wKinE = 0.0);
		
		void solveObstacleAvoidance(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo, double wKinE = 0.0);
		
		void goalFunctionWeightedVelocities(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt);
		void solveWeightedVelocities(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, double wU, double wdsdt, double wdqdt);
		
		void solve2();
		
		double setWeight(double t, double nominalMotionTime, double wMin, double wMax);	// TODO: modify it to be more general
		double setWeight(double t, double tStart, double tEnd, double wStart, double wEnd);	// overloaded version of setWeight
		
		void setOptionsQPOASES( Options myOptions );
		
		virtual void jacobian( MatrixD &J, const MatrixD &q ) = 0;
		virtual void forwardKin( MatrixD &T, const MatrixD &q ) = 0;
		virtual void errorEE( MatrixD &e, const MatrixD &x, const MatrixD &q ) = 0;
		virtual void trajectoryGenerator( MatrixD &x, MatrixD &dxdt, double t ) = 0;
		virtual void inertiaMatrix( MatrixD &M, const MatrixD &q ) = 0;
		virtual void obstacleJacobianAndVelocity( MatrixD &Jo, double &dxdto, double &woSmoothing, const MatrixD &q, const MatrixD &obstacles ) = 0;
		virtual double obstacleDistance( const MatrixD &q, const MatrixD &obstacles ) = 0;
		
		void testJacobian( const MatrixD &q );
		void testForwardKin( const MatrixD &q );
		void testErrorEE( const MatrixD &x, const MatrixD &q );
		void testStateMatrices();
		void testConstraintsAndGoalFunction();
		void testTrajectory(MatrixD &x, MatrixD &dxdt, double t);
		
	protected:
		// number of:
		int n;	// robot DOF
		int m;	// main task DOF
		int p;	// prediction steps
		int nV;	// control variables
		int nC;	// constraints
		
		double dt;			// cycle time
		double tHorizon;	// length of the prediction horizon in time
		
		MatrixD K;			// gains for CLIK (for each predicted step p)
		
		MatrixD A, B, A0, B0, A1, B1, B1t;	// state matrices
		
		MatrixD qLim, dqdtLim, d2qdt2Lim;	// joint limits (qLim.col(0) = qMin, qLim.col(1) = qMax, etc.)
		// IMPORTANT: right now, joint limits have to be symmetric
		
		MatrixD l;	// vector defining how each of the p steps is long in term of dt (i.e., l(2) = 3 means that dt2 = 3 * dt)
		
		double dsdtMax;		// maximum path velocity; dsdtMax = 1 is the nominal velocity, dsdtMax < 1 forces a slowdown
		
		bool minimizeKineticEnergy{false};	// whether to include the minimization of the kinetic energy in the goal function or not
		
		MatrixD QMin, QMax, dQdtMin, dQdtMax, d2Qdt2Min, d2Qdt2Max;	// joint limits for p predictions steps (i.e., QMin = [qLim.col(0); qLim.col(0); ... qLim.col(0)], QMin.size() = n*p)
		
		MatrixD lb, ub;			// lower and upper bounds of control variables: lb <= xi <= ub
		MatrixD lbC, ubC, C;	// lower and upper constraints, and constraint matrix: lbC <= xi <= ubC
		MatrixD Alpha, Beta;	// matrices needed for viability constraints computation
		
		MatrixD H, h;			// matrix and vector defining the goal function: 0.5 * xi' * H * 0.5 + xi' * h
		
		
		real_t *_lb, *_ub, *_C, *_lbC, *_ubC, *_H, *_h;	// pointers to arrays storing data of the above matrices
														// (qpOASES cannot use Eigen matrices directly, so _X = X.data())
		
		SQProblem *ikQP;		// QP solver from qpOASES library
		
		// qpOASES parameters
		real_t cputime0;		// for CPU time measurement
		int_t nWSRmax;			// maximum number of working set recalculations
		
		bool qpInitialized;		// signals whether SQProblem had already been solved once (true -- use hotstart() method) or not (false -- use init() method)
		
		// other stuff
		MatrixD Wvel;			// weight matrix for the weighted velocities task
};

} /* namespace PQPIKS */

#endif /* PQPIKS_H_ */
