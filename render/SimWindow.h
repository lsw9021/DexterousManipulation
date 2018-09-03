#ifndef __VMCON_SIM_WINDOW_H__
#define __VMCON_SIM_WINDOW_H__
#include "Camera.h"
#include "GLUTWindow.h"
#include "GLfunctions.h"
#include "DART_interface.h"
#include "FEM_interface.h"
#include <string>
class IntegratedWorld;
class SimWindow : public GUI::GLUTWindow
{
public:
	SimWindow(const std::string& record_path);

	std::shared_ptr<IntegratedWorld> mWorld;
	Eigen::VectorXd					mDisplacement;
	int 						mFrame;

	bool 						mIsPlay;
	bool 						mIsRotate;
	bool 						mIsDrag;
	bool 						mRenderDetail;
	bool						mCapture;
	double						mDisplayRatio;
	
	void LoadFromFolder(const std::string& path);
	void WriteTransformation(const std::string& path);
protected:
	void Display() override;
	void Keyboard(unsigned char key,int x,int y) override;
	void Mouse(int button, int state, int x, int y) override;
	void Motion(int x, int y) override;
	void Reshape(int w, int h) override;
	void Timer(int value) override;
	void Screenshot();
};

#endif