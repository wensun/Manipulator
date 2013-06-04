#ifndef _LQGMP_
#define _LQGMP_

#define _CRT_RAND_S
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <windows.h>
#include <algorithm>
#include <ppl.h>
#include <omp.h>
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "RRT.h"

class LQGMP{

public:
	struct PathNodeLQG{
		Matrix<3> T;
		Matrix<3> u; 

		Matrix<3,3> A; 
		Matrix<3,3> B;
		Matrix<3,3> V;
		Matrix<3,3> L;
		Matrix<3,2> K;
		Matrix<2,3> H;

		Matrix<6,6> F;
		Matrix<6,5> G;
		Matrix<6,6> Sigma;

		Matrix<6> y;
		Matrix<6,6> R;

		PathNodeLQG(){
			H.reset();
			T.reset();
			u.reset();
			A.reset();
			B.reset();
			V.reset();
			L.reset();
			K.reset();
			F.reset();
			G.reset();
			Sigma.reset();
			y.reset();
			R.reset();
		}
	};

	double dt;
	Matrix<3,3> P0; //initial covariance.
	Matrix<2,2> N; //sense noise
	Matrix<3,3> M; //process noise;
	Matrix<3,3> C;
	Matrix<3,3> D;
	std::vector<PathNodeLQG> m_rrtpathLQG;
	std::vector<PathNodeLQG> pathlqg;
	std::vector<RRT::PathNode> m_rrtpath;
	std::vector<PathNodeLQG> m_solutionLQG;
	double len;

	LQGMP(){};
	LQGMP(const std::vector<RRT::PathNode>& rawpath, const double& timestep, const Matrix<3,3>& initialCov){
		//m_rrtpath.clear();
		m_rrtpath = rawpath;
		dt = timestep;
		P0 = initialCov;

		Dynamics dy(dt);
		len = dy.len;

		int size = (int)m_rrtpath.size();
		pathlqg.clear();

		M = identity<3>() * 0.02*0.02;
		N = identity<2>() * 0.05*0.05;
		
		C = identity<3>();
		D = identity<3>();
	}

	void createABVLK();
	double computeConfidence(const Matrix<3>& cpos, const Matrix<3,3>& cR, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	void draw_prior_distribution(const int& cal_ellipse);
};


#endif