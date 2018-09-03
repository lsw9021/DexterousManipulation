#ifndef __IK_OPTIMIZATION_H__
#define __IK_OPTIMIZATION_H__
#include "dart/dart.hpp"
#include "dart/simulation/simulation.hpp"
#include <Eigen/Geometry>
#include <vector>
#include <utility>
#include <memory>
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>


typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
class IKOptimization : public Ipopt::TNLP
{
public:
	void ChangeSkeleton(const dart::dynamics::SkeletonPtr& skeleton);
	void AddTargetPositions(AnchorPoint ap,const Eigen::Vector3d& target);
	const std::vector<std::pair<AnchorPoint,Eigen::Vector3d>>& GetTargets();
	IKOptimization(const dart::dynamics::SkeletonPtr& skeleton);
	void ClearTarget();
	const Eigen::VectorXd& GetSolution();
	void SetSolution(const Eigen::VectorXd& sol);
	~IKOptimization();

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
									IKOptimization();
									IKOptimization(const IKOptimization&);
									IKOptimization& operator=(const IKOptimization&);

	dart::dynamics::SkeletonPtr								mSkeleton;
	Eigen::VectorXd											mSavePositions;
	Eigen::VectorXd											mInitialPositions;
	double 													w_reg;
	std::vector<std::pair<AnchorPoint,Eigen::Vector3d>> 	mTargets;
	Eigen::VectorXd											mSolution;
};

#endif