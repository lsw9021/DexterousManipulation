#ifndef __MUSCLE_OPTIMIZATION_H__
#define __MUSCLE_OPTIMIZATION_H__
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include "fem/fem.h"
#include "dart/dart.hpp"
#include "dart/simulation/simulation.hpp"

class IntegratedWorld;
class MuscleOptimization : public Ipopt::TNLP
{
	public:
										MuscleOptimization(const std::shared_ptr<IntegratedWorld>& iw);
										~MuscleOptimization();

	void								Update(const Eigen::VectorXd& qdd_desired);
	const Eigen::VectorXd&				GetSolution();	
	public:
	 bool					get_nlp_info(	Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
														Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) override;

	 bool					get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
														Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;	

	 bool					get_starting_point(	Ipopt::Index n, bool init_x, Ipopt::Number* x,
															bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
															Ipopt::Index m, bool init_lambda,
															Ipopt::Number* lambda) override;

	 bool					eval_f(	Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

	 bool					eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

	 bool					eval_g(	Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

	 bool					eval_jac_g( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
													Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
													Ipopt::Number* values) override;

	 bool					eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
												Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
												bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
												Ipopt::Index* jCol, Ipopt::Number* values) override;

	 void 					finalize_solution(	Ipopt::SolverReturn status,
															Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
															Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
															Ipopt::Number obj_value,
															const Ipopt::IpoptData* ip_data,
															Ipopt::IpoptCalculatedQuantities* ip_cq) override;
	protected:
										MuscleOptimization();
										MuscleOptimization(const MuscleOptimization&);
										MuscleOptimization& operator=(const MuscleOptimization&);

	void								UpdateConstraints(const Eigen::VectorXd& act);

	std::shared_ptr<IntegratedWorld>		mIntegratedWorld;

	double								mWeightTracking,mWeightEffort;

	Eigen::VectorXd						mQddDesired;
	Eigen::VectorXd						mSolution;

	Eigen::MatrixXd						mJt;
	Eigen::SparseMatrix<double>			mA;
	Eigen::VectorXd						mP;

	Eigen::MatrixXd						mM_minus_JtA;
	Eigen::VectorXd						mJtp_minus_c;

	int 								mSparseUpdateCount;
};
#endif