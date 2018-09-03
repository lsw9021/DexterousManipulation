#include "JugglingInfo.h"
#include <iostream>
JugglingInfo::
JugglingInfo(const std::vector<int>& V_sequences,int ball_size)
	:V(V_sequences),left("HandL"),right("HandR"),count(0),T(0.3),D(0.7)
{
	//Generate ball index.
	ball_index.resize(V.size(),-1);

	
	for(int i = 0;i<ball_size;i++)
	{
		int count_i = i;
		bool is_full = true;
		for(int j =i;j<V.size();j++)
		{
			if(ball_index[j]==-1)
			{
				count_i = j;
				is_full = false;
				break;
			}
		}

		if(is_full)
			break;
		
		while(ball_index[count_i]==-1)
		{
			ball_index[count_i] = i;
			if(V[count_i]==0)
				break;
			count_i += V[count_i];
			if(count_i>=V.size())
				break;
		}
	}

	std::cout<<"ball : ";
	for(auto b : ball_index)
		std::cout<<b<<" ";
	std::cout<<std::endl;

	count = 0;
}
void
JugglingInfo::
CountPlusPlus()
{
	count++;
}
void
JugglingInfo::
CountMinusMinus()
{
	count--;
}

Eigen::Vector3d
JugglingInfo::
GetTargetVelocity(const Eigen::Vector3d& from,const Eigen::Vector3d& to)
{	
	Eigen::Vector3d diff = to-from;
	double t_free = GetT_free();
	
	double dh = from[1] - to[1];
	double v_0y = 0.5*9.81*t_free - dh/t_free;

	double v_0x = diff[0]/t_free;
	double v_0z = diff[2]/t_free;

	return Eigen::Vector3d(v_0x,v_0y,v_0z);
}

std::string 
JugglingInfo::
From()
{
	if(count%2 == 0 )
		return right;
	else
		return left;
}
std::string 
JugglingInfo::
To()
{
	int next_count = count+GetV();
	
	if(next_count%2==0)
		return right;
	else
		return left;
}
