
#include <array>
#include <fstream>
#include "CTMIO.h"
#include "vtkPolyData.h"
//#include "base/Vector3d.h"
//#include "vtkPointData.h"
//#include "vtkDoubleArray.h"
//#include <vtkOBJWriter.h>
//#include <vtkSTLWriter.h>
//#include <vtkCellArray.h>
//#include <vtkPointLocator.h>
//#include <map>
//#include <vtkCleanPolyData.h>
//#include <vtkFloatArray.h>



CTMIO::CTMIO()
{
	
}

CTMIO::~CTMIO()
{
}


/*
*   Public interface: it's used to compress polydata to ctm file by other modules.
*/
bool CTMIO::writeCTMFile(vtkPolyData *polydata, const std::string fileName)
{
	if (nullptr == polydata ||0 == polydata->GetNumberOfPoints() ||0 == polydata->GetNumberOfCells())
	{
		return false;
	}
	vtkPoints *points = polydata->GetPoints();
	CTMuint aVertCount = static_cast<unsigned int>(points->GetNumberOfPoints());
	CTMfloat *aVertices = new CTMfloat[3 * aVertCount];

	for (CTMuint i = 0; i < aVertCount; ++i)
	{
		CVector3d pt(points->GetPoint(i));
		for (CTMuint j = 0; j < 3; ++j)
		{
			aVertices[3 * i + j] = static_cast<CTMfloat>(pt[j]);
		}
	}
	CTMuint aTriCount = static_cast<CTMuint>(polydata->GetNumberOfCells());

	CTMuint *aIndices = new CTMuint[3 * aTriCount];
	int index = 0;
	for (CTMuint i = 0; i < aTriCount; ++i)
	{
		vtkIdList *ids = polydata->GetCell(i)->GetPointIds();
		if (ids->GetNumberOfIds() != 3)
			continue;
		for (CTMuint j = 0; j < 3; ++j)
		{
			aIndices[3 * index + j] = static_cast<CTMuint>(ids->GetId(j));
		}
		index++;
	}
	CTMfloat* aUVCoords;;
	vtkDoubleArray* uv = (vtkDoubleArray*)polydata->GetPointData()->GetTCoords();
	if (uv != NULL)
	{
		aUVCoords = new CTMfloat[2 * aVertCount];
		for (int i = 0; i < polydata->GetNumberOfPoints(); i++)
		{
			double* uvt = uv->GetTuple2(i);
			aUVCoords[2 * i + 0] = uvt[0];
			aUVCoords[2 * i + 1] = uvt[1];
		}
	}

	
	CTMexporter ctmExporter;
	
	try
	{
		// Define our mesh representation to OpenCTM (store references to it in
		// the context)
		ctmExporter.DefineMesh(aVertices, aVertCount, aIndices, index, NULL);
		if (uv != NULL) {
			CTMenum mm = ctmExporter.AddUVMap(aUVCoords, "texture", NULL);
		}
		// Save the OpenCTM file
		ctmExporter.Save(fileName.c_str());
		
	}
	catch (std::exception &e)
	{
		return false;
	}


	delete[] aVertices;
	aVertices = nullptr;
	delete[] aIndices;
	aIndices = nullptr;
	return true;
}

//����ʹ��
//bool CTMIO::readCTMFile(const std::string fileName, vtkPolyData* polydata)
//{
//	//CTMimporter ctmImporter;
//	//try
//	//{
//	//	ctmImporter.Load(fileName.c_str());
//	//}
//	//catch (CTMenum& e)
//	//{
//	//	return false;
//	//}
//
//	//CTMuint aTriCount = ctmImporter.GetInteger(CTM_TRIANGLE_COUNT);
//	//CTMuint aVertCount = ctmImporter.GetInteger(CTM_VERTEX_COUNT);
//	//const CTMfloat* aVertices = ctmImporter.GetFloatArray(CTM_VERTICES);
//	//const CTMuint* aIndices = ctmImporter.GetIntegerArray(CTM_INDICES);
//	//CTMenum uvMenum = ctmImporter.GetNamedUVMap("texture");
//	//const CTMfloat* aUVCoords;
//	//if (uvMenum != CTM_NONE)
//	//	aUVCoords = ctmImporter.GetFloatArray(CTM_UV_MAP_1);
//	//vSPNew(points, vtkPoints);
//	//for (CTMuint i = 0; i < aVertCount; ++i)
//	//{
//	//	double pt[3];
//	//	for (CTMuint j = 0; j < 3; ++j)
//	//	{
//	//		pt[j] = aVertices[3 * i + j];
//	//	}
//	//	points->InsertNextPoint(pt);
//	//}
//	//vSPNew(cellArray, vtkCellArray);
//	//for (CTMuint i = 0; i < aTriCount; ++i)
//	//{
//	//	vtkIdType cell[3];
//	//	for (CTMuint j = 0; j < 3; ++j)
//	//	{
//	//		cell[j] = aIndices[3 * i + j];
//	//	}
//	//	cellArray->InsertNextCell(3, cell);
//	//}
//	//polydata->SetPoints(points);
//	//polydata->SetPolys(cellArray);
//	//polydata->Modified();
//
//
//
//
//	//if (uvMenum != CTM_NONE)
//	//{
//	//	vSPNew(dataArray, vtkDoubleArray);
//	//	dataArray->SetName("texture");
//	//	dataArray->SetNumberOfComponents(2);
//	//	dataArray->SetNumberOfTuples(polydata->GetNumberOfPoints());
//
//	//	for (CTMuint i = 0; i < aVertCount; ++i)
//	//	{
//	//		dataArray->SetTuple2(i, aUVCoords[2 * i + 0], aUVCoords[2 * i + 1]);
//	//	}
//	//	polydata->GetPointData()->SetTCoords(dataArray);
//	//}
//
//	//return true;
//
//	CTMimporter ctmImporter;
//	try
//	{
//		ctmImporter.Load(fileName.c_str());
//	}
//	catch (CTMenum& e)
//	{
//		return false;
//	}
//
//	CTMuint aTriCount = ctmImporter.GetInteger(CTM_TRIANGLE_COUNT);
//	CTMuint aVertCount = ctmImporter.GetInteger(CTM_VERTEX_COUNT);
//	const CTMfloat* aVertices = ctmImporter.GetFloatArray(CTM_VERTICES);
//	const CTMuint* aIndices = ctmImporter.GetIntegerArray(CTM_INDICES);
//	CTMenum uvMenum = ctmImporter.GetNamedUVMap("texture");
//	const CTMfloat* aUVCoords = nullptr;
//	if (uvMenum != CTM_NONE)
//		aUVCoords = ctmImporter.GetFloatArray(CTM_UV_MAP_1);
//
//	// ���� vtkPoints �� vtkCellArray
//	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
//
//	// ����㵽 points
//	for (CTMuint i = 0; i < aVertCount; ++i)
//	{
//		double pt[3];
//		for (CTMuint j = 0; j < 3; ++j)
//		{
//			pt[j] = aVertices[3 * i + j];
//		}
//		points->InsertNextPoint(pt);  // ����㵽 vtkPoints
//	}
//
//	polydata->SetPoints(points);  // �������õ� polydata ��
//	polydata->Modified();         // ���� polydata ��ȷ����ȷ����
//
//
//
//	// ���� vtkPointLocator �������ռ�����
//	vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
//	pointLocator->SetDataSet(polydata);  // �� polydata ����Ϊ���ݼ�
//	pointLocator->BuildLocator();        // ���� locator ����
//
//	std::map<CTMuint, vtkIdType> pointIdMap;  // ӳ��ԭʼ�������� vtkPoint �� ID
//
//	// ���ҺͲ���㣬�����ظ�
//	for (CTMuint i = 0; i < aVertCount; ++i)
//	{
//		double pt[3];
//		for (CTMuint j = 0; j < 3; ++j)
//		{
//			pt[j] = aVertices[3 * i + j];
//		}
//
//		// ʹ�� vtkPointLocator ��������ĵ�
//		vtkIdType existingPtId = pointLocator->FindClosestPoint(pt);
//		if (existingPtId >= 0)  // �������Ƿ�ɹ�
//		{
//			double existingPt[3];
//			points->GetPoint(existingPtId, existingPt);
//			
//			double tolerance = 1.0E-4;
//			// �����Ƿ��㹻�ӽ�������ӽ���ʹ�����еĵ� ID
//			if (vtkMath::Distance2BetweenPoints(pt, existingPt) < tolerance)
//			{
//				pointIdMap[i] = existingPtId;  // ʹ���Ѵ��ڵ�� ID
//				continue;  // ���������µ�
//			}
//		}
//
//		// ���δ�ҵ��ӽ��ĵ㣬�����µ�
//		vtkIdType newPtId = points->InsertNextPoint(pt);
//		pointIdMap[i] = newPtId;  // �����µ�� ID
//		pointLocator->InsertPoint(newPtId, pt);  // ���¶�λ��
//	}
//
//	// ���������ε�Ԫ
//	for (CTMuint i = 0; i < aTriCount; ++i)
//	{
//		vtkIdType cell[3];
//		for (CTMuint j = 0; j < 3; ++j)
//		{
//			cell[j] = pointIdMap[aIndices[3 * i + j]];  // ʹ��ӳ���ĵ� ID
//		}
//		cellArray->InsertNextCell(3, cell);
//	}
//
//	polydata->SetPolys(cellArray);  // ������ε�Ԫ���õ� polydata ��
//	polydata->Modified();           // ���� polydata
//
//	// ���� UV ����
//	if (uvMenum != CTM_NONE)
//	{
//		vtkSmartPointer<vtkDoubleArray> dataArray = vtkSmartPointer<vtkDoubleArray>::New();
//		dataArray->SetName("texture");
//		dataArray->SetNumberOfComponents(2);
//		dataArray->SetNumberOfTuples(polydata->GetNumberOfPoints());
//
//		// ���� UV ����
//		for (CTMuint i = 0; i < aVertCount; ++i)
//		{
//			vtkIdType ptId = pointIdMap[i];  // ʹ��ӳ���ĵ� ID
//			dataArray->SetTuple2(ptId, aUVCoords[2 * i + 0], aUVCoords[2 * i + 1]);
//		}
//
//		polydata->GetPointData()->SetTCoords(dataArray);  // ���� UV ����
//	}
//
//
//	return true;
//}


//ԭʼ�汾
//bool CTMIO::readCTMFile(const std::string fileName, vtkPolyData* polydata)
//{
//CTMimporter ctmImporter;
//	try
//	{
//		ctmImporter.Load(fileName.c_str());
//	}
//	catch (CTMenum& e)
//	{
//		return false;
//	}
//
//	CTMuint aTriCount = ctmImporter.GetInteger(CTM_TRIANGLE_COUNT);
//	CTMuint aVertCount = ctmImporter.GetInteger(CTM_VERTEX_COUNT);
//	const CTMfloat* aVertices = ctmImporter.GetFloatArray(CTM_VERTICES);
//	const CTMuint* aIndices = ctmImporter.GetIntegerArray(CTM_INDICES);
//	CTMenum uvMenum = ctmImporter.GetNamedUVMap("texture");
//	const CTMfloat* aUVCoords;
//	if (uvMenum != CTM_NONE)
//		aUVCoords = ctmImporter.GetFloatArray(CTM_UV_MAP_1);
//	vSPNew(points, vtkPoints);
//	for (CTMuint i = 0; i < aVertCount; ++i)
//	{
//		double pt[3];
//		for (CTMuint j = 0; j < 3; ++j)
//		{
//			pt[j] = aVertices[3 * i + j];
//		}
//		points->InsertNextPoint(pt);
//	}
//	vSPNew(cellArray, vtkCellArray);
//	for (CTMuint i = 0; i < aTriCount; ++i)
//	{
//		vtkIdType cell[3];
//		for (CTMuint j = 0; j < 3; ++j)
//		{
//			cell[j] = aIndices[3 * i + j];
//		}
//		cellArray->InsertNextCell(3, cell);
//	}
//	polydata->SetPoints(points);
//	polydata->SetPolys(cellArray);
//	polydata->Modified();
//
//
//
//
//	if (uvMenum != CTM_NONE)
//	{
//		vSPNew(dataArray, vtkDoubleArray);
//		dataArray->SetName("texture");
//		dataArray->SetNumberOfComponents(2);
//		dataArray->SetNumberOfTuples(polydata->GetNumberOfPoints());
//
//		for (CTMuint i = 0; i < aVertCount; ++i)
//		{
//			dataArray->SetTuple2(i, aUVCoords[2 * i + 0], aUVCoords[2 * i + 1]);
//		}
//		polydata->GetPointData()->SetTCoords(dataArray);
//	}
//
//	return true;
//}




bool CTMIO::readCTMFile(const std::string fileName, vtkPolyData* polydata)
{
	CTMimporter importer;

	// Load the CTM file
	importer.Load(fileName.c_str());

	// Get the mesh information
	CTMuint vertexCount = importer.GetInteger(CTM_VERTEX_COUNT);
	const CTMfloat* vertices = importer.GetFloatArray(CTM_VERTICES);
	CTMuint triCount = importer.GetInteger(CTM_TRIANGLE_COUNT);
	const CTMuint* indices = importer.GetIntegerArray(CTM_INDICES);

	CTMenum uvMenum = importer.GetNamedUVMap("texture");

	const CTMfloat* aUVCoords = nullptr;
	if (uvMenum != CTM_NONE)
		aUVCoords = importer.GetFloatArray(CTM_UV_MAP_1);

	// Create vtkPoints and add the vertices
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (CTMuint i = 0; i < vertexCount; ++i) {
		points->InsertNextPoint(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2]);
	}

	// Create vtkCellArray and add the triangles
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	for (CTMuint i = 0; i < triCount; ++i) {
		triangles->InsertNextCell(3);
		triangles->InsertCellPoint(indices[i * 3]);
		triangles->InsertCellPoint(indices[i * 3 + 1]);
		triangles->InsertCellPoint(indices[i * 3 + 2]);
	}

	// Create vtkPolyData and set the points and triangles
	vSP<vtkPolyData> tempPd = vSP<vtkPolyData>::New();
	tempPd->SetPoints(points);
	tempPd->SetPolys(triangles);

	//vSP<vtkCleanPolyData> cleanedPd = vSP<vtkCleanPolyData>::New();
	//cleanedPd->SetInputData(tempPd);
	//cleanedPd->Update();

	
	polydata->DeepCopy(tempPd);
	polydata->Modified();

	// ���� UV ����
	if (uvMenum != CTM_NONE)
	{
		vtkSmartPointer<vtkDoubleArray> dataArray = vtkSmartPointer<vtkDoubleArray>::New();
		dataArray->SetName("texture");
		dataArray->SetNumberOfComponents(2);
		dataArray->SetNumberOfTuples(polydata->GetNumberOfPoints());

		// ���� UV ����
		for (CTMuint i = 0; i < vertexCount; ++i)
		{
			dataArray->SetTuple2(i, aUVCoords[2 * i + 0], aUVCoords[2 * i + 1]);
		}
		polydata->GetPointData()->SetTCoords(dataArray);  // ���� UV ����
	}
	return true;
}


bool CTMIO::writeCTMFile2(vtkPolyData* polyData, const std::string fileName)
{
	if (!polyData) {
		std::cerr << "Invalid vtkPolyData!" << std::endl;
		return false;
	}

	// Create CTM context
	CTMcontext context = ctmNewContext(CTM_EXPORT);

	// Extract points (vertices) from vtkPolyData
	vtkPoints* points = polyData->GetPoints();
	vtkIdType numPoints = points->GetNumberOfPoints();
	std::vector<CTMfloat> vertices(3 * numPoints);

	for (vtkIdType i = 0; i < numPoints; i++) {
		double p[3];
		points->GetPoint(i, p);
		vertices[3 * i + 0] = static_cast<CTMfloat>(p[0]);
		vertices[3 * i + 1] = static_cast<CTMfloat>(p[1]);
		vertices[3 * i + 2] = static_cast<CTMfloat>(p[2]);
	}

	// Extract faces (triangles) from vtkPolyData
	vtkCellArray* cells = polyData->GetPolys();
	vtkIdType numTriangles = cells->GetNumberOfCells();
	std::vector<CTMuint> indices(3 * numTriangles);

	vtkIdType triIndex = 0;

	vtkIdType npts;
	vtkIdType const* pts;

	cells->InitTraversal();
	while (cells->GetNextCell(npts, pts)) 
	{
		if (npts == 3) {  // Only handle triangles
			indices[3 * triIndex + 0] = static_cast<CTMuint>(pts[0]);
			indices[3 * triIndex + 1] = static_cast<CTMuint>(pts[1]);
			indices[3 * triIndex + 2] = static_cast<CTMuint>(pts[2]);
			triIndex++;
		}
	}

	// Define the mesh with vertices and indices
	ctmDefineMesh(context, vertices.data(), numPoints, indices.data(), numTriangles, nullptr);

	// Check if the vtkPolyData contains UV texture coordinates
	vtkDoubleArray* uvArray = dynamic_cast<vtkDoubleArray*>(polyData->GetPointData()->GetTCoords());
	if (uvArray) 
	{
		// Extract UV coordinates
		std::vector<CTMfloat> uvCoords(2 * numPoints);  // Two floats per vertex (U and V)
		for (vtkIdType i = 0; i < numPoints; i++) {
			double uv[2];
			uvArray->GetTuple(i, uv);
			uvCoords[2 * i + 0] = static_cast<CTMfloat>(uv[0]);  // U coordinate
			uvCoords[2 * i + 1] = static_cast<CTMfloat>(uv[1]);  // V coordinate
		}

		// Define UV coordinates in the CTM context
		//ctmDefineUVCoords(context, uvCoords.data());
		ctmAddUVMap(context, uvCoords.data(), "texture", nullptr);  //CTMenum mm = ctmExporter.AddUVMap(aUVCoords, "texture", NULL);
	}
	else 
	{
		std::cerr << "No UV texture coordinates found!" << std::endl;
	}

	/*
	vtkDoubleArray* uv = (vtkDoubleArray*)polydata->GetPointData()->GetTCoords();
	if (uv != NULL)
	{
		aUVCoords = new CTMfloat[2 * aVertCount];
		for (int i = 0; i < polydata->GetNumberOfPoints(); i++)
		{
			double* uvt = uv->GetTuple2(i);
			aUVCoords[2 * i + 0] = uvt[0];
			aUVCoords[2 * i + 1] = uvt[1];
		}
	}
	*/


	// Save the CTM file
	ctmSave(context, fileName.c_str());

	// Free the CTM context
	ctmFreeContext(context);
	return true;
}