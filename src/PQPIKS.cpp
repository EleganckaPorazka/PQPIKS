#include "PQPIKS/PQPIKS.h"

namespace PQPIKS
{
PQPIKS::PQPIKS( PQPIKSdata inputData )
: n(inputData.n), m(inputData.m), p(inputData.p), dt(inputData.dt), tHorizon(inputData.tHorizon), K(inputData.K), qLim(inputData.qLim), dqdtLim(inputData.dqdtLim), d2qdt2Lim(inputData.d2qdt2Lim), dsdtMax(inputData.dsdtMax), minimizeKineticEnergy(inputData.minimizeKineticEnergy), cputime0(inputData.cputime0), nWSRmax(inputData.nWSRmax)
{	
	nV = (n + 1) * p;
	nC = (m + 3 * n) * p;
	
	if ( (qLim.rows() != n) or (dqdtLim.rows() != n) or (d2qdt2Lim.rows() != n) or (qLim.cols() != 2) or (dqdtLim.cols() != 2) or (d2qdt2Lim.cols() != 2) )
	{
		// TODO: implement some error handling
		std::cout << "Error: wrong size of qLim or dqdtLim, or d2qdt2Lim. The right size of each matrix is n x 2.\n";
	}
	if ( (K.rows() != p) or (K.cols() != 1) )
	{
		// TODO: implement some error handling
		std::cout << "Error: wrong size of K. The right size is p x 1.\n";
	}
	
	stateMatrices();
	
	initConstraintsAndGoal();
	
	ikQP = new SQProblem(nV, nC);
	
	qpInitialized = false;
}

PQPIKS::~PQPIKS()
{
}

void PQPIKS::stateMatrices()
{
	/* Computes state matrices: MatrixD A, B, A0, B0, A1, B1 */
	
	A = MatrixD::Zero(2*n, 2*n);
	B = MatrixD::Zero(2*n, n);
	
	int N = std::ceil(tHorizon / dt); // number of prediction steps with uniform distribution (each step's length is equal to dt)
	
	// state matrices A and B
	A.block(0, 0, n, n) = MatrixD::Identity(n, n);
	A.block(0, n, n, n) = MatrixD::Identity(n, n) * dt;
	A.block(n, 0, n, n) = MatrixD::Zero(n, n);
	A.block(n, n, n, n) = MatrixD::Identity(n, n);
	B.block(0, 0, n, n) = MatrixD::Identity(n, n) * 0.5 * dt * dt;
	B.block(n, 0, n, n) = MatrixD::Identity(n, n) * dt;
	
	// state matrices for all N states AN and BN
	MatrixD AN = MatrixD::Zero(2*n*N, 2*n);
	MatrixD BN = MatrixD::Zero(2*n*N, n*N);
	
	for (int i = 0; i < N; ++i)
	{
		if (i == 0)
		{
			AN.block(2*n*i, 0, 2*n, 2*n) = A;
			BN.block(2*n*i, 0, 2*n, n) = B;
		}
		else
		{
			AN.block(2*n*i, 0, 2*n, 2*n) = AN.block(2*n*(i - 1), 0, 2*n, 2*n) * A;
			BN.block(2*n*i, 0, 2*n, n*N) = A * BN.block(2*n*(i - 1), 0, 2*n, n*N);
			BN.block(2*n*i, n*i, 2*n, n) = B;
		}
	}
	
	// selection matrices S0 and S1
	MatrixD S0 = MatrixD::Zero(n*N, 2*n*N);
	MatrixD S1 = MatrixD::Zero(n*N, 2*n*N);
	
	for (int i = 0; i < N; ++i)
	{
		S0.block(n*i, 2*n*i, n, n) = MatrixD::Identity(n, n);
		S1.block(n*i, 2*n*i + n, n, n) = MatrixD::Identity(n, n);
	}
	
	// reducing the size of the optimal control problem (Ssize) and implementing input blocking (Sblocking)
	MatrixD Ssize, Sblocking;
	if (p < N)
	{
		stateMatricesSizeReduction(Ssize, Sblocking, N);
	}
	else // p == N
	{
		l = MatrixD::Ones(p, 1);
		Ssize = MatrixD::Identity(n*N, n*N);
		Sblocking = MatrixD::Identity(n*N, n*N);
	}
	
	A0 = Ssize * S0 * AN;
	B0 = Ssize * S0 * BN * Sblocking;
	A1 = Ssize * S1 * AN;
	B1 = Ssize * S1 * BN * Sblocking;
	B1t = B1.transpose();
}

void PQPIKS::stateMatricesSizeReduction(MatrixD &Ssize, MatrixD &Sblocking, int N)
{
	// reducing the size of the optimal control problem (Ssize) and implementing input blocking (Sblocking)
	
	MatrixD tk(p, 1);	// vector of not uniformly distributed time instants; tk = [t1, t2, ..., tp]
	MatrixD Nk(p, 1);	// vector of uniform time steps leading to each time instant; Nk(i) = tk(i) / dt
	l = MatrixD::Ones(p, 1);	// vector defining how each of the p steps is long in term of dt (i.e., l(2) = 3 means that dt2 = 3 * dt)
	
	double a = (tHorizon - dt) / ( (p - 1) * (p - 1) );	// constant which will be useful below
	for (int k = 0; k < p; ++k)
	{
		tk(k) = a * (k + 1) * (k + 1) - 2 * a * (k + 1) + a + dt;	// not uniform time instances distribution
		Nk(k) = std::round(tk(k) / dt);			// it may cause problems if used not as intended, i.e. with p = N, like here: dt = 0.02, tHorizon = 0.1, p = 5
		if (k > 0)
		{
			l(k) = Nk(k) - Nk(k - 1);
		}
	}
		
	// compute Sblocking and Ssize
	Sblocking = MatrixD::Zero(n*N, n*p);
	int row = 0;
	int col = 0;
	for (int k = 0; k < p; ++k)
	{
		col = n*k;
		for (int i = 0; i < l(k); ++i)
		{
			Sblocking.block(row, col, n, n) = MatrixD::Identity(n, n);
			row += n;
		}
	}
	
	Ssize = MatrixD::Zero(n*p, n*N);
	for (int k = 0; k < p; ++k)
	{
		row = n*k;
		col = n*(Nk(k) - 1);
		Ssize.block(row, col, n, n) = MatrixD::Identity(n, n);
	}
	
}

void PQPIKS::initConstraintsAndGoal()
{
	// matrices and vectors to implement the viability constraints
	// IMPORTANT: symmetric joint limits assumed
	MatrixD alpha = MatrixD::Zero(n, 1);
	MatrixD beta = MatrixD::Zero(n, 1);
	for (int j = 0; j < n; ++j)
	{
		alpha(j) = -d2qdt2Lim(j, 1) / dqdtLim(j, 1);
		beta(j) = d2qdt2Lim(j, 1) / dqdtLim(j, 1) * qLim(j, 1);
	}
	
	Alpha = MatrixD::Zero(n*p, n*p);
	Beta = MatrixD::Zero(n*p, 1);
	for (int k = 0; k < p; ++k)
	{
		Alpha.block(n*k, n*k, n, n) = alpha.asDiagonal();
		Beta.block(n*k, 0, n, 1) = beta;
	}
		
	QMin = MatrixD::Zero(n*p, 1);
	QMax = MatrixD::Zero(n*p, 1);
	dQdtMin = MatrixD::Zero(n*p, 1);
	dQdtMax = MatrixD::Zero(n*p, 1);
	d2Qdt2Min = MatrixD::Zero(n*p, 1);
	d2Qdt2Max = MatrixD::Zero(n*p, 1);
	for (int k = 0; k < p; ++k)
	{
		// joint position constraints
		QMin.block(n*k, 0, n, 1) = qLim.col(0);
		QMax.block(n*k, 0, n, 1) = qLim.col(1);
		// joint velocity constraints
		dQdtMin.block(n*k, 0, n, 1) = dqdtLim.col(0);
		dQdtMax.block(n*k, 0, n, 1) = dqdtLim.col(1);
		// joint acceleration constraints
		d2Qdt2Min.block(n*k, 0, n, 1) = d2qdt2Lim.col(0);
		d2Qdt2Max.block(n*k, 0, n, 1) = d2qdt2Lim.col(1);
	}
	
	// lower and upper bounds of control variables: lb <= xi <= ub
	lb = MatrixD::Zero(nV, 1);
	ub = MatrixD::Zero(nV, 1);
	lb.block(0, 0, n*p, 1) = d2Qdt2Min;
	lb.block(n*p, 0, p, 1) = MatrixD::Zero(p, 1);
	ub.block(0, 0, n*p, 1) = d2Qdt2Max;
	ub.block(n*p, 0, p, 1) = MatrixD::Ones(p, 1) * dsdtMax;
	
	C = MatrixD::Zero(nC, nV);		// constraint matrix
	// constant part of the constraint matrix
	C.block(m*p, 0, n*p, n*p) = B0;							// joint position constraints
	C.block((m + n)*p, 0, n*p, n*p) = B1;					// joint velocity constraints
	C.block((m + 2*n)*p, 0, n*p, n*p) = B1 - Alpha * B0;	// viability constraints
	
	// lower and upper constraints: lbC <= C * xi <= ubC
	lbC = MatrixD::Zero(nC, 1);
	ubC = MatrixD::Zero(nC, 1);
	
	// goal function components: 0.5 * xi' * H * 0.5 + xi' * h
	H = MatrixD::Zero(nV, nV);
	h = MatrixD::Zero(nV, 1);
	
	// prepare data for the qpOASES solver
	_lb = lb.data();
	_ub = ub.data();
	_C = C.data();
	_lbC = lbC.data();
	_ubC = ubC.data();
	_H = H.data();
	_h = h.data();
	
	// init the matrix Wvel for the weighted velocities task
	Wvel = MatrixD::Zero(n*p, n*p);
	double smallestMaxVel = dqdtLim.col(1).minCoeff();
	for (int i = 0; i < (n*p); ++i)
	{
		Wvel(i, i) = ( smallestMaxVel / dQdtMax(i) ) * ( smallestMaxVel / dQdtMax(i) );
	}
}

void PQPIKS::constraints(const MatrixD &z0, const MatrixD &QEst, const MatrixD &dsdt, double s)
{

	MatrixD JEst = MatrixD::Zero(m*p, n*p);
	MatrixD DXDSEst = MatrixD::Zero(m*p, p);
	MatrixD dXdsEst = MatrixD::Zero(m*p, 1);
	
	MatrixD J(m, n), x(m, 1), dxds(m, 1), e(m, 1);
	double tEst = s;	// time in the predicted step
	for (int k = 0; k < p; ++k)
	{
		jacobian(J, QEst.block(n*k, 0, n, 1));
		JEst.block(m*k, n*k, m ,n) = J;
		tEst += l(k) * dt * dsdt(k);
		trajectoryGenerator(x, dxds, tEst);
		errorEE(e, x, QEst.block(n*k, 0, n, 1));
		DXDSEst.block(m*k, k, m, 1) = dxds + K(k) * e;
	}
	dXdsEst = DXDSEst * MatrixD::Ones(p, 1);
	
	// equality constraints
	C.block(0, 0, m*p, n*p) = JEst * B1;
	C.block(0, n*p, m*p, p) = DXDSEst;
	//~ lbC.block(0, 0, m*p, 1) = dXdsEst - JEst * A1 * z0;
	lbC.block(0, 0, m*p, 1) = ubC.block(0, 0, m*p, 1) = dXdsEst - JEst * A1 * z0;
	
	// inequality constraints
	// joint position constraints
	lbC.block(m*p, 0, n*p, 1) = QMin - A0 * z0;
	ubC.block(m*p, 0, n*p, 1) = QMax - A0 * z0;
	// joint velocity constraints
	lbC.block((m + n)*p, 0, n*p, 1) = dQdtMin - A1 * z0;
	ubC.block((m + n)*p, 0, n*p, 1) = dQdtMax - A1 * z0;
	// viability constraints
	lbC.block((m + 2*n)*p, 0, n*p, 1) = (Alpha * A0 - A1) * z0 - Beta;
	ubC.block((m + 2*n)*p, 0, n*p, 1) = (Alpha * A0 - A1) * z0 + Beta;
}

//~ void PQPIKS::goalFunction(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt)
//~ {
	//~ // goal function components
	
	//~ if (minimizeKineticEnergy)
	//~ {
		//~ MatrixD MEst = MatrixD::Zero(n*p, n*p);
		//~ MatrixD M(n, n);
		//~ for (int k = 0; k < p; ++k)
		//~ {
			//~ inertiaMatrix(M, QEst.block(n*k, 0, n, 1));
			//~ MEst.block(n*k, n*k, n, n) = M;
		//~ }
		//~ H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + wdqdt * B1.transpose() * MEst * B1;
		//~ H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
		//~ h.block(0, 0, n*p, 1) = wdqdt * B1.transpose() * MEst * A1 * z0;
	//~ }
	//~ else
	//~ {
		//~ H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + wdqdt * B1.transpose() * B1;
		//~ H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
		//~ h.block(0, 0, n*p, 1) = wdqdt * B1.transpose() * A1 * z0;
	//~ }
//~ }

void PQPIKS::goalFunction(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt)
{
	// goal function components
	
	H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + wdqdt * B1t * B1;
	H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
	h.block(0, 0, n*p, 1) = wdqdt * B1t * A1 * z0;
}

void PQPIKS::goalFunction(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt, double wKinE)
{
	// goal function components
	// overloaded version for kinetic energy minimization
	
	MatrixD MEst = MatrixD::Zero(n*p, n*p);
	MatrixD M(n, n);
	for (int k = 0; k < p; ++k)
	{
		inertiaMatrix(M, QEst.block(n*k, 0, n, 1));
		MEst.block(n*k, n*k, n, n) = M;
	}
	H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + B1t * ( wKinE * MEst + wdqdt * MatrixD::Identity(n*p, n*p) ) * B1;
	H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
	h.block(0, 0, n*p, 1) = B1t * ( wKinE * MEst + wdqdt * MatrixD::Identity(n*p, n*p) ) * A1 * z0;
}

void PQPIKS::solve(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, double wU, double wdsdt, double wdqdt, double wKinE)
{
	// based on the current state z0, and previous solution of U and dsdt, compute the predicted q, J, dxds for the next p steps
	MatrixD QEst = A0 * z0 + B0 * U;
	
	constraints(z0, QEst, dsdt, s);
	if (minimizeKineticEnergy)
	{
		goalFunction(z0, QEst, wU, wdsdt, wdqdt, wKinE);
	}
	else
	{
		goalFunction(z0, QEst, wU, wdsdt, wdqdt);
	}
	
	real_t xiOpt[nV];	// array for the QP solution
	// reset the cputime and nWSR parameters to default values:
	real_t cputime = cputime0;
	int_t nWSR = nWSRmax;
	
	if (qpInitialized == false)
	{
		ikQP->init( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
		qpInitialized = true;
	}
	else
	{
		ikQP->hotstart( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
	}
	
	s += dsdt(0) * dt;
	
	// get the results
	for (int i = 0; i < n*p; ++i)
	{
		U(i) = xiOpt[i];
	}
	for (int i = 0; i < p; ++i)
	{
		dsdt(i) = 1.0 - xiOpt[i + n*p];
	}
	
	// update the initial conditions for the next step
	z0 = A * z0 + B * U.block(0, 0, n, 1);
	//~ s += dsdt(0) * dt;
}

void PQPIKS::goalFunctionObstacleAvoidance(const MatrixD &z0, const MatrixD &QEst, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo)
{
	// goal function components including obstacle avoidance task
	
	MatrixD JoEst = MatrixD::Zero(p, n*p);
	MatrixD dXdtoEst(p, 1);
	
	MatrixD Jo = MatrixD::Zero(1, n);
	double dxdto;
	
	MatrixD woSmoothing = MatrixD::Zero(p, p);
	double woSmoothingFactor;
	
	for (int k = 0; k < p; ++k)
	{
		obstacleJacobianAndVelocity(Jo, dxdto, woSmoothingFactor, QEst.block(n*k, 0, n, 1), obstacles);
		JoEst.block(k, n*k, 1 ,n) = Jo;
		dXdtoEst(k, 0) = dxdto;
		woSmoothing(k, k) = woSmoothingFactor;
	}
	
	H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + wdqdt * B1t * B1 + wo * B1t * JoEst.transpose() * woSmoothing * JoEst * B1;
	H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
	h.block(0, 0, n*p, 1) = wdqdt * B1t * A1 * z0 + wo * B1t * JoEst.transpose() * woSmoothing * (JoEst * A1 * z0 - dXdtoEst);
}

void PQPIKS::goalFunctionObstacleAvoidance(const MatrixD &z0, const MatrixD &QEst, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo, double wKinE)
{
	// goal function components including obstacle avoidance task
	// overloaded version for kinetic energy minimization
	
	MatrixD JoEst = MatrixD::Zero(p, n*p);
	MatrixD dXdtoEst(p, 1);
	
	MatrixD Jo = MatrixD::Zero(1, n);
	double dxdto;
	
	MatrixD woSmoothing = MatrixD::Zero(p, p);
	double woSmoothingFactor;
	
	MatrixD MEst = MatrixD::Zero(n*p, n*p);
	MatrixD M(n, n);
	for (int k = 0; k < p; ++k)
	{
		obstacleJacobianAndVelocity(Jo, dxdto, woSmoothingFactor, QEst.block(n*k, 0, n, 1), obstacles);
		JoEst.block(k, n*k, 1 ,n) = Jo;
		//~ dXdtoEst.block(k, 0, 1, 1) = dxdto;
		dXdtoEst(k, 0) = dxdto;
		woSmoothing(k, k) = woSmoothingFactor;
		
		inertiaMatrix(M, QEst.block(n*k, 0, n, 1));
		MEst.block(n*k, n*k, n, n) = M;
	}
	H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + B1t * ( wKinE * MEst + wdqdt * MatrixD::Identity(n*p, n*p) ) * B1 + wo * B1t * JoEst.transpose() * woSmoothing * JoEst * B1;
	H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
	h.block(0, 0, n*p, 1) = B1t * ( wKinE * MEst + wdqdt * MatrixD::Identity(n*p, n*p) ) * A1 * z0 + wo * B1t * JoEst.transpose() * woSmoothing * (JoEst * A1 * z0 - dXdtoEst);
}

void PQPIKS::solveObstacleAvoidance(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, const MatrixD &obstacles, double wU, double wdsdt, double wdqdt, double wo, double wKinE)
{
	// based on the current state z0, and previous solution of U and dsdt, compute the predicted q, J, dxds for the next p steps
	MatrixD QEst = A0 * z0 + B0 * U;
	
	constraints(z0, QEst, dsdt, s);
	if (minimizeKineticEnergy)
	{
		goalFunctionObstacleAvoidance(z0, QEst, obstacles, wU, wdsdt, wdqdt, wo, wKinE);
	}
	else
	{
		goalFunctionObstacleAvoidance(z0, QEst, obstacles, wU, wdsdt, wdqdt, wo);
	}
	
	real_t xiOpt[nV];	// array for the QP solution
	// reset the cputime and nWSR parameters to default values:
	real_t cputime = cputime0;
	int_t nWSR = nWSRmax;
	
	if (qpInitialized == false)
	{
		ikQP->init( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
		qpInitialized = true;
	}
	else
	{
		ikQP->hotstart( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
	}
	
	s += dsdt(0) * dt;
	
	// get the results
	for (int i = 0; i < n*p; ++i)
	{
		U(i) = xiOpt[i];
	}
	for (int i = 0; i < p; ++i)
	{
		dsdt(i) = 1.0 - xiOpt[i + n*p];
	}
	
	// update the initial conditions for the next step
	z0 = A * z0 + B * U.block(0, 0, n, 1);
	//~ s += dsdt(0) * dt;
}

void PQPIKS::goalFunctionWeightedVelocities(const MatrixD &z0, const MatrixD &QEst, double wU, double wdsdt, double wdqdt)
{
	// goal function components
	
	H.block(0, 0, n*p, n*p) = wU * MatrixD::Identity(n*p, n*p) + wdqdt * B1t * Wvel * B1;
	H.block(n*p, n*p, p, p) = wdsdt * MatrixD::Identity(p, p);
	h.block(0, 0, n*p, 1) = wdqdt * B1t * Wvel * A1 * z0;
}

void PQPIKS::solveWeightedVelocities(MatrixD &z0, MatrixD &U, MatrixD &dsdt, double &s, double wU, double wdsdt, double wdqdt)
{
	// based on the current state z0, and previous solution of U and dsdt, compute the predicted q, J, dxds for the next p steps
	MatrixD QEst = A0 * z0 + B0 * U;
	
	constraints(z0, QEst, dsdt, s);
	
	goalFunctionWeightedVelocities(z0, QEst, wU, wdsdt, wdqdt);
	
	real_t xiOpt[nV];	// array for the QP solution
	// reset the cputime and nWSR parameters to default values:
	real_t cputime = cputime0;
	int_t nWSR = nWSRmax;
	
	if (qpInitialized == false)
	{
		ikQP->init( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
		qpInitialized = true;
	}
	else
	{
		ikQP->hotstart( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
		ikQP->getPrimalSolution( xiOpt );
	}
	
	s += dsdt(0) * dt;
	
	// get the results
	for (int i = 0; i < n*p; ++i)
	{
		U(i) = xiOpt[i];
	}
	for (int i = 0; i < p; ++i)
	{
		dsdt(i) = 1.0 - xiOpt[i + n*p];
	}
	
	// update the initial conditions for the next step
	z0 = A * z0 + B * U.block(0, 0, n, 1);
	//~ s += dsdt(0) * dt;
}

void PQPIKS::solve2()
{
	//~ TODO: write it properly, decide what should be class member and what an input data
	
	// penalty factors
	// for control vector U
	double wU = 1000.0;
	// for scaling factor dsdt
	double wdsdt = 10000.0;
	// for joint velocity vector dqdt
	//~ double wdqdt = 1000.0;
	double wdqdt = 1.0;
	
	double pi = M_PI;
	
	MatrixD z0(2*n, 1), U(n*p, 1), dsdt(p, 1);
	z0 << -pi/2.0, 0.0, 0.0, pi/2.0, 0.0, -pi/2.0, 0.0, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001;
	U = MatrixD::Zero(n*p, 1);
	dsdt = MatrixD::Ones(p, 1);
	double s = 0.0;
	
	// based on the current state z0, and previous solution of U and dsdt, compute the predicted q, J, dxds for the next p steps
	MatrixD QEst = A0 * z0 + B0 * U;
	
	constraints(z0, QEst, dsdt, s);
	goalFunction(z0, QEst, wU, wdsdt, wdqdt);
	
	/* Solve first QP ... */
	real_t cputime = cputime0;
	int_t nWSR = nWSRmax;
	ikQP->init( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
	real_t xiOpt[nV];
	ikQP->getPrimalSolution( xiOpt );
	std::cout << "xi* =\n";
	for (int i = 0; i < nV; ++i)
	{
		std::cout << xiOpt[i] << "\n";
	}
	std::cout << "\n";
	std::cout << "CPU time: " << cputime*1.0e6 << " microseconds\n\n";
	
	// get the results
	for (int i = 0; i < n*p; ++i)
	{
		U(i) = xiOpt[i];
	}
	for (int i = 0; i < p; ++i)
	{
		dsdt(i) = 1.0 - xiOpt[i + n*p];
	}
	
	// update the initial conditions for the next step
	z0 = A * z0 + B * U.block(0, 0, n, 1);
	s += dsdt(0) * dt;
	
	// based on the current state z0, and previous solution of U and dsdt, compute the predicted q, J, dxds for the next p steps
	QEst = A0 * z0 + B0 * U;
	
	constraints(z0, QEst, dsdt, s);
	goalFunction(z0, QEst, wU, wdsdt, wdqdt);
	
	/* Solve second QP */
	cputime = cputime0;
	nWSR = nWSRmax;
	ikQP->hotstart( _H, _h, _C, _lb, _ub, _lbC, _ubC, nWSR, &cputime );
	ikQP->getPrimalSolution( xiOpt );
	std::cout << "xi* =\n";
	for (int i = 0; i < nV; ++i)
	{
		std::cout << xiOpt[i] << "\n";
	}
	std::cout << "\n";
	std::cout << "CPU time: " << cputime*1.0e6 << " microseconds\n\n";
}

double PQPIKS::setWeight(double t, double nominalMotionTime, double wMin, double wMax)
{
	// TODO: make it more general
	
	/* Set the weight (penalty factor) for the task */
	
	double tSwitch = 0.9 * nominalMotionTime;
	
	double w;
	if (t < tSwitch)
	{
		w = wMin;
	}
	else if (t > nominalMotionTime)
	{
		w = wMax;
	}
	else
	{
		MatrixD a(3, 1);
		a << wMin, 3.0 * (wMax - wMin) / pow( (nominalMotionTime - tSwitch), 2.0 ), -2.0 * (wMax - wMin) / pow( (nominalMotionTime - tSwitch), 3.0 );
		w = a(0) + a(1) * pow( (t - tSwitch), 2.0 ) + a(2) * pow( (t - tSwitch), 3.0 );
	}
	return w;
}

double PQPIKS::setWeight(double t, double tStart, double tEnd, double wStart, double wEnd)
{
	/* Set the weight (penalty factor) for the task */
	
	double w;
	if (t < tStart)
	{
		w = wStart;
	}
	else if (t > tEnd)
	{
		w = wEnd;
	}
	else
	{
		MatrixD a(3, 1);
		a << wStart, 3.0 * (wEnd - wStart) / pow( (tEnd - tStart), 2.0 ), -2.0 * (wEnd - wStart) / pow( (tEnd - tStart), 3.0 );
		w = a(0) + a(1) * pow( (t - tStart), 2.0 ) + a(2) * pow( (t - tStart), 3.0 );
	}
	return w;
}

void PQPIKS::setOptionsQPOASES( Options myOptions )
{
	ikQP->setOptions( myOptions );
}


void PQPIKS::testJacobian( const MatrixD &q )
{
	//~ Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> J(m, n);
	MatrixD J(m, n);
	jacobian(J, q);
	std::cout << "J = \n" << J << "\n";
	std::cout << "J = \n";
	double *Jdata = J.data();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			std::cout << Jdata[i*n + j] << "\t";
		}
		std::cout << "\n";
	}
}

void PQPIKS::testForwardKin( const MatrixD &q )
{
	MatrixD T(4, 4);
	forwardKin(T, q);
	std::cout << "T = \n" << T << "\n";
}

void PQPIKS::testErrorEE( const MatrixD &x, const MatrixD &q )
{
	MatrixD e;
	errorEE(e, x, q);
	std::cout << "e = \n" << e << "\n";
}

void PQPIKS::testStateMatrices()
{
	std::cout << "A =\n" << A << "\n";
	std::cout << "B =\n" << B << "\n";
	std::cout << "A0 =\n" << A0 << "\n";
	std::cout << "B0 =\n" << B0 << "\n";
	std::cout << "A1 =\n" << A1 << "\n";
	std::cout << "B1 =\n" << B1 << "\n";
}

void PQPIKS::testConstraintsAndGoalFunction()
{
	std::cout << "ub =\n" << ub << "\n";
	std::cout << "lb =\n" << lb << "\n";
	std::cout << "ubC =\n" << ubC << "\n";
	std::cout << "lbC =\n" << lbC << "\n";
	std::cout << "C =\n" << C << "\n";
	std::cout << "H =\n" << H << "\n";
	std::cout << "h =\n" << h << "\n";
}

void PQPIKS::testTrajectory(MatrixD &x, MatrixD &dxdt, double t)
{
	trajectoryGenerator(x, dxdt, t);
	std::cout << "x = [" << x.transpose() << "], dxdt = [" << dxdt.transpose() << "]\n";
}

} /* namespace PQPIKS */
