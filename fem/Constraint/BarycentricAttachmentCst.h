#ifndef __BARYCENTRIC_ATTACHMENT_CST_H__
#define __BARYCENTRIC_ATTACHMENT_CST_H__

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


class BarycentricAttachmentCst : public Cst
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
	const Eigen::Vector3d& GetP() {return mP;}
	void SetP(const Eigen::Vector3d& p) {mP =p;}
	BarycentricAttachmentCst(const BarycentricAttachmentCst& other) = delete;
	BarycentricAttachmentCst& operator=(const BarycentricAttachmentCst& other) = delete;
	std::shared_ptr<Cst> Clone() override;
	static std::shared_ptr<BarycentricAttachmentCst> Create(const std::string& name,double k,const Eigen::Vector3d& p,const Eigen::Vector4d& c,int i0,int i1,int i2,int i3);
protected:
	BarycentricAttachmentCst(const std::string& name,double k,const Eigen::Vector3d& p,const Eigen::Vector4d& c,int i0,int i1,int i2,int i3);

	int mi0,mi1,mi2,mi3;
	Eigen::Vector4d mC;
	Eigen::Vector3d mP;

	//For parallization
	double 		mE;
	Eigen::Vector12d mg;
	Eigen::Matrix12d mH;
};
};

#endif