#ifndef __JUGGLING_INFO_H__
#define __JUGGLING_INFO_H__
#include <vector>
#include <string>
#include <Eigen/Core>
class JugglingInfo
{	
public:
	JugglingInfo(const std::vector<int>& V_sequences,int ball_size);

	std::string 	From();
	std::string 	To();

	
	void			CountPlusPlus();
	void			CountMinusMinus();

	Eigen::Vector3d GetTargetVelocity(const Eigen::Vector3d& from,const Eigen::Vector3d& to);
	double			GetT_free(){return std::max(T*(V[count]-2*D),0.0);};
	double			GetT_hold(){return T*D;};
	int 			GetBallIndex(){return ball_index[count];};
	int 			GetV(){return V[count];};
	double 			GetT(){return T;};
	double 			GetD(){return D;};
	int 			GetCount(){return count;}
	void			SetCount(int c){ count =c;}
private:
	std::string left,right;
	double T;
	double D;
	std::vector<int> V;

	//Automatically generated
	std::vector<int> ball_index;

	int count;
};



#endif