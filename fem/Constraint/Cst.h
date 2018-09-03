#ifndef __FEM_CST_H__
#define __FEM_CST_H__
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <string>
#include <memory>
namespace FEM
{
class Cst
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	virtual void EvaluatePotentialEnergy(const Eigen::VectorXd& x) = 0;
	virtual void EvaluateGradient(const Eigen::VectorXd& x) = 0;
	virtual void EvaluateHessian(const Eigen::VectorXd& x) = 0;

	virtual void GetPotentialEnergy(double& e) = 0;
	virtual void GetGradient(Eigen::VectorXd& g) = 0;
	virtual void GetHessian(std::vector<Eigen::Triplet<double>>& h_triplets) = 0;

	virtual void EvaluateDVector(const Eigen::VectorXd& x) = 0;
	virtual void GetDVector(int& index,Eigen::VectorXd& d) = 0;
	virtual void EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets) = 0;
	virtual void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) = 0;

	virtual int GetNumHessianTriplets() = 0;
	virtual void AddOffset(int offset) = 0;

	Cst(const Cst& other) = delete;
	Cst& operator=(const Cst& other) = delete;
	virtual std::shared_ptr<Cst> Clone() = 0;

	bool Equal(const std::shared_ptr<Cst>& other){return !(other->mName.compare(this->mName));}
	std::string GetName(){return mName;};
protected:
	Cst(const std::string& name,double k):mName(name),mStiffness(k){};
	double mStiffness;
	std::string mName;

};

};
#endif