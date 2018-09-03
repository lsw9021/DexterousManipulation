#ifndef __FEM_COROTATE_FEM_CST_H__
#define __FEM_COROTATE_FEM_CST_H__
#include "../Tensor3333.h"
#include "Cst.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace Eigen {

typedef Matrix<double, 12, 1> Vector12d;
typedef Matrix<double, 12, 12> Matrix12d;
};
namespace FEM
{


class CorotateFEMCst : public Cst
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void EvaluatePotentialEnergy(const Eigen::VectorXd& x) override;
	void EvaluateGradient(const Eigen::VectorXd& x) override;
	void EvaluateHessian(const Eigen::VectorXd& x) override;
	void GetPotentialEnergy(double& e) override;
	void GetGradient(Eigen::VectorXd& g) override;
	void GetHessian(std::vector<Eigen::Triplet<double>>& h_triplets) override;
	void EvaluateDVector(const Eigen::VectorXd& x) override;
	void GetDVector(int& index,Eigen::VectorXd& d) override;
	void EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
	void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;
	int GetNumHessianTriplets() override;
	void AddOffset(int offset) override;

	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}
	int GetI3() {return mi3;}
	
	CorotateFEMCst(const CorotateFEMCst& other) = delete;
	CorotateFEMCst& operator=(const CorotateFEMCst& other) = delete;
	std::shared_ptr<Cst> Clone() override;
	static std::shared_ptr<CorotateFEMCst> Create(const std::string& name,double k,double poisson_ratio,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm);
protected:
	CorotateFEMCst(const std::string& name,double k,double poisson_ratio,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm);

	int mi0,mi1,mi2,mi3;
	double mVolume;
	double mMu,mLambda;
	double mPoissonRatio;
	Eigen::Matrix3d mInvDm;

	//For parallization
	double 		mE;
	Eigen::Vector12d mg;
	Eigen::Matrix12d mH;
	Eigen::Matrix3d md;
	Eigen::Matrix3d md_volume;

	//For Cache
	Eigen::Vector12d mX;
	Eigen::Matrix3d mF;
	Eigen::Matrix3d mR,mU,mV,mD;
	void ComputeF(const Eigen::VectorXd& x);
	void ComputeP(Eigen::Matrix3d& P);
	void ComputedPdF(Tensor3333& dPdF);
};
};
#endif