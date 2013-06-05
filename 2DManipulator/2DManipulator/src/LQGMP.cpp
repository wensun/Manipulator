#include "LQGMP.h"


void LQGMP::createABVLK()
{
	int ns = (int) m_rrtpath.size();	
	pathlqg.clear();
	for(int i = 0; i < ns; i++){
		PathNodeLQG tmpnode;
		tmpnode.T = m_rrtpath[i].T;
		tmpnode.u = m_rrtpath[i].u;
		pathlqg.push_back(tmpnode);
	}
	for(int i = 0; i < ns; i++){
		RRT::PathNode tmp = m_rrtpath[i];
		Matrix<3> u = tmp.u;
		Matrix<3> X = tmp.T;

		Matrix<3,3> A;
		A.reset();
		Matrix<3,3> B;
		B.reset();
		Matrix<3,3> V;
		V.reset();
		A = identity<3>();
		B = dt * identity<3>();

		Matrix<2,3> H; H.reset();
		H(0,0) = X[0];			H(0,1) = 0;				H(0,2) = 0; 
		H(1,0) = len*cos(X[0]); H(1,1) = len*cos(X[1]); H(1,2) = len*cos(X[2]);

		V = B;
		pathlqg[i].A = A;
		pathlqg[i].B = B;
		pathlqg[i].V = V;
		pathlqg[i].H = H;
	}

	Matrix<3,3> S = C;
	int length = (int)pathlqg.size() - 1;
	for(int k = length - 1; k != -1; k--){
		pathlqg[k].L = -!(~pathlqg[k].B*S*pathlqg[k].B + D)*~pathlqg[k].B*S*pathlqg[k].A;
		S = C + ~pathlqg[k].A*S*pathlqg[k].A + ~pathlqg[k].A*S*pathlqg[k].B*pathlqg[k].L;
	}

	Matrix<3,3> P = P0;
	for(int k = 1; k <= length; k++){
		P = pathlqg[k-1].A*P*~pathlqg[k-1].A + pathlqg[k-1].V * M *~pathlqg[k-1].V;
		
		pathlqg[k].K = P * ~pathlqg[k].H* !(pathlqg[k].H*P*~pathlqg[k].H + N);
		P = (identity<3>() - pathlqg[k].K * pathlqg[k].H) * P;
	}
	pathlqg[0].Sigma = zeros<6,6>();
	pathlqg[0].Sigma.insert(0,0, P0);

	for(int k = 1; k <= length; k++){
		Matrix<6,6> F;
		F.insert(0,0, pathlqg[k-1].A);
		F.insert(0,3, pathlqg[k-1].B * pathlqg[k-1].L);
		F.insert(3,0, pathlqg[k].K * pathlqg[k].H * pathlqg[k-1].A);
		F.insert(3,3, pathlqg[k-1].A + pathlqg[k-1].B * pathlqg[k-1].L - pathlqg[k].K * pathlqg[k-1].H * pathlqg[k-1].A);
		pathlqg[k-1].F = F;

		Matrix<6, 5> G;
		G.insert(0,0, pathlqg[k-1].V);
		G.insert(0,3, zeros<3,2>());
		G.insert(3,0, pathlqg[k].K * pathlqg[k-1].H * pathlqg[k-1].V);
		G.insert(3,3, pathlqg[k].K);
		pathlqg[k-1].G = G;
		
		Matrix<5,5> Q;
		Q.reset();
		Q.insert(0,0, M);
		Q.insert(3,3, N);
		pathlqg[k].Sigma = F * pathlqg[k-1].Sigma * ~F + G * Q * ~G;
	}
}



void LQGMP::draw_prior_distribution(const int& cal_ellipse){

	createABVLK();
	for(int i = 0; i < (int)pathlqg.size(); i++){
		drawEllipse3d(pathlqg[i].T, pathlqg[i].Sigma.subMatrix<3,3>(0,0), cal_ellipse, true);
	}

}

double LQGMP::computeWSConfidence(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	Matrix<2,2> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);
	for(int i = 0; i < 2; i++){
		EVal(i,i) += 0.0001;
	}

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

	Matrix<2,2> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<2> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

	int num_pairs;
	CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete [] results;

	CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
	CAL_SetGroupScaling(cal_environment, 1,1,1);

	return distance;
}

double LQGMP::computeConfidence(const Matrix<3>& q, const Matrix<3,3>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	double mindis = 10000;
	int mlink = -1; 
	int mpoint = -1;
	double seg = 5;
	double seglen = len / seg;
	//find the maximum likelihodd collision point. the one that has the minimum times of standard deviation (1 after transformation.)
	for(int link = 1; link <=3; link++){
		for(int point = 0; point < seg; point++){
			double length = (point+1)*seglen; //length on the particular link.
			Matrix<2,3> Joc = zeros<2,3>();
			forward_jocabi(q, link, length, Joc);
			Matrix<2> c = forward_k(q, link, length) - Joc * q;
			Matrix<2,3> A = Joc; //A*q + c
			Matrix<2> xmean = forward_k(q, link, length);
			Matrix<2,2> xcov = A * R * ~A;
			double tmpdis = computeWSConfidence(xmean, xcov, cal_obstacles, cal_environment, cal_point); 
			if(tmpdis < mindis){
				mindis = tmpdis;
				mlink = link;
				mpoint = point;
			}
		} 
	}
	return mindis;
}

double LQGMP::computeProbability(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	P0 = initialCov;
	createABVLK();
	double log_prob = 0;
	for(int i = 0; i < (int)pathlqg.size(); i++)
	{
		Matrix<3,1> q = pathlqg[i].T;
		Matrix<3,3> qR = pathlqg[i].Sigma.subMatrix<3,3>(0,0);
		double conf = computeConfidence(q, qR, cal_obstacles, cal_environment, cal_point);
		double p = log(incompletegamma(1, 0.5*conf*conf));
		log_prob += p;
	}

	return exp(log_prob);
}



