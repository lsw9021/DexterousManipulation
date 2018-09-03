#include "SimWindow.h"
#include "../vmcon/IntegratedWorld.h"
#include "../vmcon/MusculoSkeletalSystem.h"
#include "../vmcon/DART_helper.h"
#include "../vmcon/Ball.h"
#include "../vmcon/Record.h"
#include "dart/external/lodepng/lodepng.h"

#include <algorithm>
#include <fstream>
#include <boost/filesystem.hpp>
#include <GL/glut.h>
using namespace GUI;
using namespace FEM;
using namespace dart::simulation;
using namespace dart::dynamics;
SimWindow::
SimWindow(const std::string& record_path)
	:GLUTWindow(),mIsRotate(true),mIsDrag(false),mIsPlay(false),mFrame(0),mDisplayRatio(1.0),mRenderDetail(false),mCapture(false)
{
	mWorld = std::make_shared<IntegratedWorld>();
	mWorld->Initialize();
	LoadFromFolder(record_path);
	
	mDisplacement = mWorld->GetRecords()[0]->rigid_body_positions["HUMAN"];
	mDisplacement.setZero();
	mDisplayTimeout = 33;
}

void
SimWindow::
Display() 
{
	glClearColor(1.0, 1.0, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	initLights();
	mCamera->Apply();
	glDisable(GL_LIGHTING);
	// glLineWidth(2.0);
	// DrawLine(Eigen::Vector3d(0,0,0),Eigen::Vector3d(100,0,0),Eigen::Vector3d(1,0,0));
	// DrawLine(Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,100,0),Eigen::Vector3d(0,1,0));
	// DrawLine(Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,0,100),Eigen::Vector3d(0,0,1));
	// glLineWidth(1.0);
	// glColor3f(0,0,0);
	// glLineWidth(1.0);
	// glBegin(GL_LINES);
	// // for(double z = -0.2;z<=0.2;z+=0.1)
	// {
	// 	double z = 0.0;
	// 	for(double x=-2.0;x<=2.0;x+=0.1)
	// 	{
	// 		glVertex3f(x,-1.0,z);
	// 		glVertex3f(x,3.0,z);
	// 	}
	// 	for(double y=-1.0;y<=3.0;y+=0.1)
	// 	{
	// 		glVertex3f(-2.0,y,z);
	// 		glVertex3f(2.0,y,z);
	// 	}
	// }
	// glEnd();

	// glBegin(GL_LINES);
	// {
	// 	glLineWidth(5.0);
	// 	double z = 0.0;
	// 	for(double x=-2.0;x<=2.0;x+=0.1)
	// 	{
	// 		glVertex3f(x,z,-2.0);
	// 		glVertex3f(x,z,2.0);
	// 	}
	// 	for(double y=-1.0;y<=3.0;y+=0.1)
	// 	{
	// 		glVertex3f(-2.0,z,y);
	// 		glVertex3f(2.0,z,y);
	// 	}
	// 	glLineWidth(1.0);
	// }
	// glEnd();

	Eigen::Vector3d clr[5] =
	{
		Eigen::Vector3d(0.8,0.2,0.2),
		Eigen::Vector3d(0.2,0.8,0.2), 
		Eigen::Vector3d(0.2,0.2,0.8),
		Eigen::Vector3d(0.8,0.8,0.2),
		Eigen::Vector3d(0.2,0.8,0.8)
	};
	

	mWorld->GetRecords()[mFrame]->Get(mWorld->GetRigidWorld(),mWorld->GetSoftWorld(),mWorld->GetMusculoSkeletalSystem());
	// GUI::DrawStringOnScreen(0.8,0.2,std::to_string(mWorld->GetRigidWorld()->getTime()),true,Eigen::Vector3d(0,0,0));
	int ball_index = 0;
	for(int i =0;i<mWorld->GetRigidWorld()->getNumSkeletons();i++)
	{
		if(i==0)
		{
			mWorld->GetRigidWorld()->getSkeleton(i)->setPositions(mDisplacement+mWorld->GetRigidWorld()->getSkeleton(i)->getPositions());
			mWorld->GetRigidWorld()->getSkeleton(i)->computeForwardKinematics(true,false,false);
			if(mRenderDetail)
				mWorld->GetRigidWorld()->getSkeleton(i)->getPositions();
			DrawSkeleton(mWorld->GetRigidWorld()->getSkeleton(i),Eigen::Vector3d(0.8,0.8,0.8),!mRenderDetail);
		}
		else
		{
			// if(i==5)
				// continue;
			glDisable(GL_LIGHTING);
			// if(i==1)
			// {
				// glLineWidth(10.0);
				glColor3f(clr[ball_index][0],clr[ball_index][1],clr[ball_index][2]);
				
				Eigen::Vector3d prev_com = mWorld->GetRigidWorld()->getSkeleton(i)->getBodyNode(0)->getCOM();
				// glBegin(GL_LINE_STRIP);
				for(int j=0;j<60;j+=5){
					if(mFrame-j<0)
						break;
					mWorld->GetRecords()[mFrame-j]->Get(mWorld->GetRigidWorld(),mWorld->GetSoftWorld(),mWorld->GetMusculoSkeletalSystem());

					Eigen::Vector3d com = mWorld->GetRigidWorld()->getSkeleton(i)->getBodyNode(0)->getCOM();
					DrawArrow3D(com,(prev_com-com),(prev_com-com).norm(),0.005,(Eigen::Vector3d(1,1,1)-clr[ball_index])*(j/(double)60)+clr[ball_index],false);
					// glVertex3f(com[0],com[1],com[2]);
					prev_com = com;
					// DrawSkeleton(mWorld->GetRigidWorld()->getSkeleton(i),clr[ball_index]);
					// clr[ball_index]=(Eigen::Vector3d(1,1,1)-clr[ball_index])*0.1+clr[ball_index];
				}
				// glEnd();
			// }
			glEnable(GL_LIGHTING);
			// else
			// {
			mWorld->GetRecords()[mFrame]->Get(mWorld->GetRigidWorld(),mWorld->GetSoftWorld(),mWorld->GetMusculoSkeletalSystem());
			DrawSkeleton(mWorld->GetRigidWorld()->getSkeleton(i),clr[ball_index]);
			// }
			ball_index++;
		}
    
    }
	if(mRenderDetail)
	{
		DrawWorld(mWorld->GetSoftWorld());
		
		for(auto& mus: mWorld->GetMusculoSkeletalSystem()->GetMuscles()){
			DrawMuscleWayPoints(mus->origin_way_points);
			DrawMuscleWayPoints(mus->insertion_way_points);
		}
	}

	if(mCapture)
		Screenshot();
	glutSwapBuffers();
}
void
SimWindow::
Keyboard(unsigned char key,int x,int y) 
{
	switch(key)
	{
		case '-' : mDisplayRatio*=0.9;break;
		case '+' : mDisplayRatio*=1.1;break;
		// case '1' : mDisplacement[0+8] +=0.1; break;
		// case '2' : mDisplacement[1+8] +=0.1; break;
		// case '3' : mDisplacement[2+8] +=0.1; break;
		// case '4' : mDisplacement[3+8] +=0.1; break;
		// case '5' : mDisplacement[4+8] +=0.1; break;
		// case '6' : mDisplacement[5+8] +=0.1; break;
		// case '7' : mDisplacement[6+8] +=0.1; break;
		// case '8' : mDisplacement[7+8] +=0.1; break;
		// case '9' : mDisplacement[8+8] +=0.1; break;
		// case '0' : mDisplacement[9+8] +=0.1; break;
		// case 'q' : mDisplacement[0+8] -=0.1; break;
		// case 'w' : mDisplacement[1+8] -=0.1; break;
		// case 'e' : mDisplacement[2+8] -=0.1; break;
		// case 'r' : mDisplacement[3+8] -=0.1; break;
		// case 't' : mDisplacement[4+8] -=0.1; break;
		// case 'y' : mDisplacement[5+8] -=0.1; break;
		// case 'u' : mDisplacement[6+8] -=0.1; break;
		// case 'i' : mDisplacement[7+8] -=0.1; break;
		// case 'o' : mDisplacement[8+8] -=0.1; break;
		// case 'p' : mDisplacement[9+8] -=0.1; break;
		case 'a' : mRenderDetail =!mRenderDetail;break;
		case 'C' : mCapture =!mCapture;break;
		case ' ' : mIsPlay =!mIsPlay;break;
		case '`' : mIsRotate = !mIsRotate; break;
		case 'r' : mFrame = 0;break;
		case '[' : mFrame--;break;
		case ']' : mFrame++;break;
		case 27: exit(0);break;
		default : break;
	}
	std::cout<<mDisplacement.transpose()<<std::endl;
	if(mFrame<0)
		mFrame = mWorld->GetRecords().size()-1;
	if(mFrame>mWorld->GetRecords().size()-1)
		mFrame= 0;

	glutPostRedisplay();
}
void
SimWindow::
Mouse(int button, int state, int x, int y) 
{
	if (state == GLUT_DOWN)
	{
		mIsDrag = true;
		mMouseType = button;
		mPrevX = x;
		mPrevY = y;
	}
	else
	{
		mIsDrag = false;
		mMouseType = 0;
	}

	glutPostRedisplay();
}
void
SimWindow::
Motion(int x, int y) 
{
	if (!mIsDrag)
		return;

	int mod = glutGetModifiers();
	if (mMouseType == GLUT_LEFT_BUTTON)
	{
		if(!mIsRotate)
		mCamera->Translate(x,y,mPrevX,mPrevY);
		else
		mCamera->Rotate(x,y,mPrevX,mPrevY);
	}
	else if (mMouseType == GLUT_RIGHT_BUTTON)
	{
		switch (mod)
		{
		case GLUT_ACTIVE_SHIFT:
			mCamera->Zoom(x,y,mPrevX,mPrevY); break;
		default:
			mCamera->Pan(x,y,mPrevX,mPrevY); break;		
		}

	}
	mPrevX = x;
	mPrevY = y;
	glutPostRedisplay();
}
void
SimWindow::
Reshape(int w, int h) 
{
	glViewport(0, 0, w, h);
	mCamera->Apply();
}
void
SimWindow::
Timer(int value) 
{
	if(mIsPlay){
			if(mCapture)
				mFrame++;
			else
				mFrame+=(int)33.0*(mDisplayRatio/1.0);

			if(mFrame>=mWorld->GetRecords().size())
				mFrame = 0;
		}
		
	glutPostRedisplay();
	glutTimerFunc(mDisplayTimeout, TimerEvent,1);
}


void
SimWindow::
LoadFromFolder(const std::string& path)
{
	int count = 0;
	while(true)
	{
		std::string real_path = path+std::to_string(count++);
		if(!boost::filesystem::exists(real_path))
			break;

		auto new_rec = std::make_shared<Record>();
		new_rec->Read(real_path);
		mWorld->AddRecord(new_rec);
		
	}
}
void
SimWindow::
WriteTransformation(const std::string& path)
{
	boost::filesystem::create_directories(path);

	int count = 0;
	auto recs = mWorld->GetRecords();
	for(int i =0;i<recs.size();i++)
	{
		std::ofstream ofs(path+std::to_string(i));
		recs[count]->Get(mWorld->GetRigidWorld(),mWorld->GetSoftWorld(),mWorld->GetMusculoSkeletalSystem());
		for(int j=0;j<mWorld->GetRigidWorld()->getNumSkeletons();j++)
		{
			for(int k=0;k<mWorld->GetRigidWorld()->getSkeleton(j)->getNumBodyNodes();k++)
			{
				Eigen::Isometry3d T = mWorld->GetRigidWorld()->getSkeleton(j)->getBodyNode(k)->getTransform();
				ofs<<mWorld->GetRigidWorld()->getSkeleton(j)->getBodyNode(k)->getName()<<" ";
				for(int x=0;x<16;x++)
				{
					ofs<<T.data()[x]<<" ";
				}
				ofs<<std::endl;
			}
		}
		ofs<<"a "<<mWorld->GetMusculoSkeletalSystem()->GetActivationLevels().transpose()<<std::endl;

		ofs.close();
		count++;
	}


}
void SimWindow::
Screenshot() {
  static int count = 0;
  const char directory[8] = "frames";
  const char fileBase[8] = "Capture";
  char fileName[32];

  boost::filesystem::create_directories(directory);
  std::snprintf(fileName, sizeof(fileName), "%s%s%s%.4d.png",
                directory, "/", fileBase, count++);
  int tw = glutGet(GLUT_WINDOW_WIDTH);
  int th = glutGet(GLUT_WINDOW_HEIGHT);

  glReadPixels(0, 0,  tw, th, GL_RGBA, GL_UNSIGNED_BYTE, &mScreenshotTemp[0]);

  // reverse temp2 temp1
  for (int row = 0; row < th; row++) {
    memcpy(&mScreenshotTemp2[row * tw * 4],
           &mScreenshotTemp[(th - row - 1) * tw * 4], tw * 4);
  }

  unsigned result = lodepng::encode(fileName, mScreenshotTemp2, tw, th);

  // if there's an error, display it
  if (result) {
    std::cout << "lodepng error " << result << ": "
              << lodepng_error_text(result) << std::endl;
    return ;
  } else {
    std::cout << "wrote screenshot " << fileName << "\n";
    return ;
  }
}
