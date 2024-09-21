#pragma once

#include "openctmpp.h"

class vtkVectorText;
class vtkFollower;
class CUToothModel;
class CUvtkActor;
class vtkPolyData;
class vtkIntersectionPolyDataFilter;
class vtkSplineFilter;



class CTMIO
{
public:
	CTMIO();
	~CTMIO();
	bool writeCTMFile(vtkPolyData *polydata, const std::string fileName);
	bool writeCTMFile2(vtkPolyData* polydata, const std::string fileName);
	bool readCTMFile(const std::string fileName, vtkPolyData* polydata);
};

