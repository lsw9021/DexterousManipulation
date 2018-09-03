#ifndef __RECORD_H__
#define __RECORD_H__
#include <string>
#include <map>
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

class MusculoSkeletalSystem;
class Controller;

namespace FEM
{
	class World;
};

class Record
{
public:
	Record(){};
	void Set(const dart::simulation::WorldPtr& rigid_world,
			 const std::shared_ptr<FEM::World>& soft_world,
			 const std::shared_ptr<MusculoSkeletalSystem>& musculo_skeletal_system);
	void Get(const dart::simulation::WorldPtr& rigid_world,
			 const std::shared_ptr<FEM::World>& soft_world,
			 const std::shared_ptr<MusculoSkeletalSystem>& musculo_skeletal_system);
	void Write(const std::string& path);
	void Read(const std::string& path);
	double			t;
	std::map<std::string,Eigen::VectorXd> rigid_body_positions;
	std::map<std::string,Eigen::VectorXd> rigid_body_velocities;
	std::vector<std::pair<std::string,std::string>> constraints;
	Eigen::VectorXd soft_body_positions;
	Eigen::VectorXd activation_levels;
};

#endif