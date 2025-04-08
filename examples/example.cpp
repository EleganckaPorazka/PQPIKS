#include <iostream>
#include <fstream>
#include <chrono>
#include "PQPIKS/PQPIKS_MDH.h"
#include "PQPIKS/typedefs.h"

void saveResults(const std::string &OutputFileName, const MatrixD &X)
{
	// Set the precision.
	//~ std::cout << std::setprecision(7);
	
	std::fstream OutFile;
	OutFile.open(OutputFileName, std::ios::out);
	if (OutFile.is_open()) {
		OutFile << X;
		OutFile.close();
		std::cout << "File written." << std::endl;
	}
	else {
		std::cout << "Could not create file: " << OutputFileName << std::endl;
	}
};

int main()
{
	PQPIKS::PQPIKSdata inputData;
	// number of:
	inputData.n = 7;	// robot DOF
	inputData.m = 3;	// main task DOF
	inputData.p = 3;	// prediction steps
	
	int n = inputData.n;
	int p = inputData.p;
	
	inputData.dt = 0.01;		// cycle time
	inputData.tHorizon = 0.25;	// length of the prediction horizon
	
	inputData.K = MatrixD::Zero(p, 1);	// gains for CLIK
	inputData.K << 25.0, 1.0, 1.0;
	
	// robot limits
	inputData.qLim = MatrixD::Zero(n, 2);
	inputData.dqdtLim = MatrixD::Zero(n, 2);
	inputData.d2qdt2Lim = MatrixD::Zero(n, 2);
	
	MatrixD qMax(n, 1);
	qMax << 165.0, 115.0, 165.0, 115.0, 165.0, 115.0, 165.0; qMax *= M_PI/180.0;
	inputData.qLim.col(0) = -qMax;
	inputData.qLim.col(1) = qMax;
	MatrixD dqdtMax(n, 1);
	//~ dqdtMax << 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0; dqdtMax *= M_PI/180.0;
	dqdtMax << 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0; dqdtMax *= M_PI/180.0;
	inputData.dqdtLim.col(0) = -dqdtMax;
	inputData.dqdtLim.col(1) = dqdtMax;
	MatrixD d2qdt2Max(n, 1);
	//~ d2qdt2Max << 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0; d2qdt2Max *= M_PI/180.0;
	d2qdt2Max << 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0; d2qdt2Max *= M_PI/180.0;
	inputData.d2qdt2Lim.col(0) = -d2qdt2Max;
	inputData.d2qdt2Lim.col(1) = d2qdt2Max;
	
	inputData.nWSRmax = 100;
	
	PQPIKS::TrajectoryData trajectoryData;
	
	MatrixD q0(n, 1);
	q0 << 0.0, 0.0, 0.0, -M_PI/2.0, 0.0, M_PI/2.0, 0.0;
	trajectoryData.q0 = q0;
	
	MatrixD Points(3, 5);
	Eigen::Vector3d P1(0.0, 0.3, 1.0);
	Eigen::Vector3d P2(0.0, -0.3, 1.0);
	Eigen::Vector3d P3(-0.5, 0.3, 0.6);
	Eigen::Vector3d P4(-0.5, -0.3, 0.6);
	Points.col(0) = P1;
	Points.col(1) = P2;
	Points.col(2) = P3;
	Points.col(3) = P4;
	Points.col(4) = P1;
	trajectoryData.Points = Points;
	MatrixD SegmentTimes = MatrixD::Zero(1, 6);
	SegmentTimes(1) = 2.0;
	SegmentTimes(2) = 1.5;
	SegmentTimes(3) = 2.0;
	SegmentTimes(4) = 1.5;
	SegmentTimes(5) = 2.0;
	
	trajectoryData.SegmentTimes = SegmentTimes;
	trajectoryData.AccelerationPeriod = 0.4;
	trajectoryData.DecelerationPeriod = 0.4;
	trajectoryData.BlendingPeriod = 0.15;
	
	Eigen::Vector3d angles_EE; angles_EE << 0.0, M_PI/2.0, M_PI/2.0;
	Eigen::Vector3d r_EE; r_EE << 0.1, 0.0, 0.078;
	
	//	Modified Denavit-Hartenberg parameters:
	MatrixD MDHangles(n, 1); // angles
	MDHangles << 0.0, pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2;
	MatrixD MDHlengthsX(n, 1); // lengths on X axis
	MDHlengthsX << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	MatrixD MDHlengthsZ(n, 1); // lengths on Z axis
	MDHlengthsZ << 0.31, 0.0, 0.4, 0.0, 0.39, 0.0, 0.0;
	
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
	
	//~ add end effector as a point mass m_EE = 2 kg at the r_EE; parameters modified according to the above paper
	double m_EE = 2;
	double d3 = MDHlengthsZ(2);
	double d5 = MDHlengthsZ(4);
	AllFirstMomentsOfMass(1, 1) += d3 * m_EE;
	AllFirstMomentsOfMass(1, 3) -= d5 * m_EE;
	AllFirstMomentsOfMass(1, 5) += r_EE(2) * m_EE;
	AllFirstMomentsOfMass(0, 6) += r_EE(0) * m_EE;
	AllFirstMomentsOfMass(1, 6) += r_EE(1) * m_EE;
	AllInertias(0, 3) += d3 * d3 * m_EE;
	AllInertias(0, 9) += d5 * d5 * m_EE;
	AllInertias(2, 11) += d5 * d5 * m_EE;
	AllInertias(0, 18) += r_EE(0) * r_EE(0) * m_EE;
	AllInertias(2, 20) += r_EE(2) * r_EE(2) * m_EE;
	
	MatrixD Gravity(3, 1);
	Gravity << 0.0, 0.0, -9.80665;
	
	//~ inputData.minimizeKineticEnergy = true;
	bool minimizeKineticEnergy = inputData.minimizeKineticEnergy;
	
	PQPIKS::PQPIKS_MDH myIK(inputData, trajectoryData, MDHangles, MDHlengthsX, MDHlengthsZ, AllRotAxes, Masses, AllFirstMomentsOfMass, AllInertias, Gravity, angles_EE, r_EE);
	
	MatrixD T = MatrixD::Zero(4, 4);
	MatrixD dqdt0 = MatrixD::Zero(n, 1);
	MatrixD z0(2*n, 1);
	z0.block(0, 0, n, 1) = q0;
	z0.block(n, 0, n, 1) = dqdt0;
	MatrixD U = MatrixD::Zero(n*p, 1);
	MatrixD dsdt = MatrixD::Ones(p, 1);
	double s = 0.0;
	double t = 0.0;
	double sEnd = myIK.getNominalMotionTime();
	int steps = 1500;
	int step = 0;
	
	// matrices for results
	MatrixD Q(steps, n), dQdt(steps, n), d2Qdt2(steps, n), S(steps, 1), dSdt(steps, 1), timeVec(steps, 1), msrTime(steps, 1), X(steps, 3), Xdes(steps, 3), kinE(steps, 1);
	
	MatrixD x(3, 1), dxds(3, 1);
	
	Options myOptions;	// qpOASES options
    // myOptions.setToMPC();
	myOptions.printLevel = PL_NONE;
	myIK.setOptionsQPOASES(myOptions);
	
	// penalty factors
	// for control vector U
	double wU = 1.0;
	// for scaling factor dsdt
	double wdsdt = 10000.0;
	// for joint velocity vector dqdt
	double wdqdtMin = 0.0, wdqdtMax = 0.0, wdqdt = 0.0, wKinE = 0.0;
	if (minimizeKineticEnergy)
	{
		wKinE = 800.0;
		wdqdtMin = 1.0;
		wdqdtMax = 1000.0;
		wdqdt = wdqdtMin;
	}
	else
	{
		wdqdt = 1000.0;
	}
	
	double nearZero = 0.00001;	// a small eps value to check whether joint velocities are sufficiently small to consider them zero
	MatrixD vel = z0.bottomRows(n).transpose() * z0.bottomRows(n);
	while ( (s <= sEnd) or (vel(0) > nearZero) )
	{
		//~ std::cout << "Step: " << step << "\n";
		// save z0, etc
		Q.block(step, 0, 1, n) = z0.block(0, 0, n, 1).transpose();
		dQdt.block(step, 0, 1, n) = z0.block(n, 0, n, 1).transpose();
		S(step) = s;
		timeVec(step) = t;
		myIK.forwardKin(T, Q.block(step, 0, 1, n).transpose());
		X.block(step, 0, 1, 3) = T.block(0, 3, 3, 1).transpose();
		myIK.trajectoryGenerator(x, dxds, s);
		Xdes.block(step, 0, 1, 3) = x.transpose();
		MatrixD M(n, n);
		myIK.inertiaMatrix(M, z0.block(0, 0, n, 1));
		kinE.block(step, 0, 1, 1) = 0.5 * z0.block(n, 0, n, 1).transpose() * M * z0.block(n, 0, n, 1);
		
		std::chrono::high_resolution_clock::time_point Tic = std::chrono::high_resolution_clock::now();
		if (minimizeKineticEnergy)
		{
			wdqdt = myIK.setWeight(s, 0.9*sEnd, sEnd, wdqdtMin, wdqdtMax);
			myIK.solve(z0, U, dsdt, s, wU, wdsdt, wdqdt, wKinE);
		}
		else
		{
			myIK.solve(z0, U, dsdt, s, wU, wdsdt, wdqdt);
		}
		std::chrono::high_resolution_clock::time_point Toc = std::chrono::high_resolution_clock::now();
		auto WallClockTime = std::chrono::duration_cast<std::chrono::microseconds>( Toc - Tic ).count();
		
		// save d2qdt2, dsdt
		d2Qdt2.block(step, 0, 1, n) = U.block(0, 0, n, 1).transpose();
		dSdt(step) = dsdt(0);
		msrTime(step) = WallClockTime;
		++step;
		t += inputData.dt;
		
		//~ std::cout << "Step: " << step << "\n";
		//~ std::cout << "z0 = " << z0.transpose() << "\n";
		//~ std::cout << "s = " << s << "\n";
		
		vel = z0.bottomRows(n).transpose() * z0.bottomRows(n);
		
		if (step >= steps)
		{
			s = 2 * sEnd;
			vel(0) = nearZero;
			break;
		}
	}
	
	std::cout << "Nominal motion time: " << sEnd << " seconds.\n";
	
	Q.conservativeResize(step, n);
	dQdt.conservativeResize(step, n);
	d2Qdt2.conservativeResize(step, n);
	S.conservativeResize(step, 1);
	dSdt.conservativeResize(step, 1);
	timeVec.conservativeResize(step, 1);
	msrTime.conservativeResize(step, 1);
	X.conservativeResize(step, 3);
	Xdes.conservativeResize(step, 3);
	kinE.conservativeResize(step, 1);
	
	saveResults("Q.txt", Q);
	saveResults("dQdt.txt", dQdt);
	saveResults("d2Qdt2.txt", d2Qdt2);
	saveResults("S.txt", S);
	saveResults("dSdt.txt", dSdt);
	saveResults("t.txt", timeVec);
	saveResults("wclock.txt", msrTime);
	saveResults("X.txt", X);
	saveResults("Xdes.txt", Xdes);
	saveResults("kinE.txt", kinE);
	
	return 0;
}
