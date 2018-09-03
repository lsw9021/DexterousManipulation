#ifndef __FEM_ATTACHMENT_CST_H__
#define __FEM_ATTACHMENT_CST_H__
#include "Cst.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{
class AttachmentCst : public Cst
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
	const Eigen::Vector3d& GetP() {return mp;}
	void SetP(const Eigen::Vector3d& p) {mp = p;}

	AttachmentCst(const AttachmentCst& other) = delete;
	AttachmentCst& operator=(const AttachmentCst& other) = delete;
	std::shared_ptr<Cst> Clone() override;
	static std::shared_ptr<AttachmentCst> Create(const std::string& name,double k,int i0,const Eigen::Vector3d& p);
protected:
	AttachmentCst(const std::string& name,double k,int i0,const Eigen::Vector3d& p);
	int mi0;
	Eigen::Vector3d mp;

	//For parallization
	double 		mE;
	Eigen::Vector3d mg;
	Eigen::Matrix3d mH;
	Eigen::Vector3d md;
};
}
#endif