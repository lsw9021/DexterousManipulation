#ifndef __FEM_TENSOR_3333_H__
#define __FEM_TENSOR_3333_H__
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{
class Tensor3333
{
public:
	Tensor3333();
	Tensor3333(const Tensor3333& other);
	Tensor3333& operator=(const Tensor3333& other);

	Tensor3333 operator+() const;
	Tensor3333 operator-() const;
	Tensor3333 operator+(const Tensor3333& B) const;
	Tensor3333 operator-(const Tensor3333& B) const;
	Tensor3333 operator*(const Eigen::Matrix3d& m) const;
	Tensor3333 operator*(double a) const;

	Eigen::Matrix3d& operator()(int i, int j);

	void SetIdentity();
	void SetZero();
	Tensor3333 Transpose();

public:
	Eigen::Matrix3d A[3][3];
};

Tensor3333 operator*(double a,const Tensor3333& B);
Tensor3333 operator*(const Eigen::Matrix3d& m,const Tensor3333& B);
std::ostream& operator<<(std::ostream& os,const Tensor3333& B);
}
#endif