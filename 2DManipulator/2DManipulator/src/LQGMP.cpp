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

double LQGMP::computeConfidence(const Matrix<3>& cpos, const Matrix<3,3>& cR, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	Matrix<2> p1 = forward_k(cpos, 1, len);
	Matrix<2> p2 = forward_k(cpos, 2, len);
	Matrix<2> p3 = forward_k(cpos, 3, len);
	int np[1] = {2};
	float line[18] = {0.0, 0.0, 0.0, p1[0], p1[1], 0.0, p1[0], p1[1], 0.0, p2[0], p2[1], 0.0, p2[0], p2[1], 0.0, p3[0], p3[1], 0.0};
	CAL_CreatePolyline(cal_point, 3, np, line);

	int num_pairs;
	CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
}

void LQGMP::draw_prior_distribution(const int& cal_ellipse){

	createABVLK();
	for(int i = 0; i < (int)pathlqg.size(); i++){
		drawEllipse2d(pathlqg[i].T.subMatrix<2,1>(0,0), pathlqg[i].Sigma.subMatrix<2,2>(0,0), cal_ellipse, false);
	}

}


