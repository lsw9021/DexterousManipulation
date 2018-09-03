#ifndef __BEZIER_CURVE_H__
#define __BEZIER_CURVE_H__

class BezierCurve
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	BezierCurve(){};
	BezierCurve(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& p2,double T)
		:mp0(p0),mp1(p1),mp2(p2),mInvT(1.0/T)
	{

	}

	void Initialize(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& p2,double T)
	{
		mp0 = p0;
		mp1 = p1;
		mp2 = p2;
		mInvT = 1.0/T;
	}
	Eigen::Vector3d GetPosition(double t) const
	{
		double s = t*mInvT;
		return
			mp0*(1-s)*(1-s)+
			mp1*2*s*(1-s)+
			mp2*s*s;
	}
public:
	Eigen::Vector3d mp0,mp1,mp2;
	double mInvT;

};
class BezierCurved
{
public:
	BezierCurved(){};
	void Initialize(double p0,double p1,double p2,double T)
	{
		mp0 = p0;
		mp1 = p1;
		mp2 = p2;
		mInvT = 1.0/T;
	}
	double GetPosition(double t) const
	{
		double s = t*mInvT;
		return
			mp0*(1-s)*(1-s)+
			mp1*2*s*(1-s)+
			mp2*s*s;
	}
public:
	double mp0,mp1,mp2;
	double mInvT;
};

class BSplined
{
public:
	BSplined(){};
	void Initialize(double p0,double p1,double p2,double p3,double T)
	{
		mp0 = p0;
		mp1 = p1;
		mp2 = p2;
		mp3 = p3;
		mInvT = 1.0/T;
	}
	double GetPosition(double t) const
	{
		double s = t*mInvT;
		double b0,b1,b2,b3;
		b0 = 1.0/6.0*(1-s)*(1-s)*(1-s);
		b1 = 1.0/6.0*(3*s*s*s - 6*s*s +4);
		b2 = 1.0/6.0*(-3*s*s*s +3*s*s +3*s +1);
		b3 = 1.0/6.0*s*s*s;
		return mp0*b0 + mp1*b1 + mp2*b2 + mp3*b3;
	}
public:
	double mp0,mp1,mp2,mp3;
	double mInvT;	
};
#endif