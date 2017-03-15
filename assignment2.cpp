// Tên:  Vu Ngoc Linh
// MSSV: 51201929
#include "stdafx.h"
#include <Windows.h>
#include <math.h>
#include <GL.H>
#include <glut.h>
#include <iostream>
using namespace std;

#define PI			3.1415926
#define	COLORNUM	14
#define DEG2RAD		3.1415926f/180.0f

float	ColorArr[COLORNUM][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, { 0.0,  0.0, 1.0}, // red, green, blue
								 {1.0, 1.0,  0.0}, { 1.0, 0.0, 1.0},{ 0.0, 1.0, 1.0}, 
								 {0.3, 0.3, 0.3}, {0.5, 0.5, 0.5}, { 0.9,  0.9, 0.9},
								 {1.0, 0.5,  0.5}, { 0.5, 1.0, 0.5},{ 0.5, 0.5, 1.0},
								 {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};

class Point3 {
public:
	float x, y, z;

public:
	void set(float dx, float dy, float dz) {
		x = dx;
		y = dy;
		z = dz;
	}
	void set(Point3& p) {
		x = p.x;
		y = p.y;
		z = p.z;
	}

	// Constructor and Destructor
	Point3() { x = y = z = 0; };
	Point3(float dx, float dy, float dz) {
		x = dx;
		y = dy;
		z = dz;
	}
};

class Color3 {
public:
	float r, g, b;

public:
	void set(float red, float green, float blue) {
		r = red;
		g = green;
		b = blue;
	}
	void set(Color3& c) {
		r = c.r;
		g = c.g;
		b = c.b;
	}

	// Constructor and Destructor
	Color3() { r = g = b = 0; }
	Color3(float red, float green, float blue) {
		r = red;
		g = green;
		b = blue;
	}
};

class Vector3 {
public:
	float x, y, z;

public:
	void set(float dx, float dy, float dz) {
		x = dx;
		y = dy;
		z = dz;
	}
	void set(Vector3& v) {
		x = v.x;
		y = v.y;
		z = v.z;
	}
	void flip() {
		x = -x;
		y = -y;
		z = -z;
	}
	void normalize() {
		float temp = sqrt(x*x + y*y + z*z);
		x = x/temp;
		y = y/temp;
		z = z/temp;
	}
	Vector3 cross(Vector3 b) {
		Vector3 c(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
		return c;
	}
	float dot(Vector3 b) {
		return x*b.x + y*b.y + z*b.z;
	}

	// Constructor and Destructor
	Vector3() { x = y = z = 0; }
	Vector3(float dx, float dy, float dz) {
		x = dx;
		y = dy;
		z = dz;
	}
	Vector3(Vector3& v) {
		x = v.x;
		y = v.y;
		z = v.z;
	}
};

class VertexID {
public :
	int		vertIndex;
	int		colorIndex;
};

class Face {
public:
	int			nVerts;
	VertexID*	vert;
	Vector3		facenorm;

public:
	Face() {
		nVerts = 0;
		vert = NULL;
	}
	~Face() {
		if (vert != NULL) {
			delete[] vert;
			vert = NULL;
		}
		nVerts = 0;
	}
};

class Mesh {
public:
	int		numVerts;
	Point3* pt;

	int		numFaces;
	Face*	face;

	float	slideX, slideY, slideZ;
	float	rotateX, rotateY, rotateZ;

public:
	// CONSTRUCTOR AND DESTRUCTOR
	Mesh() {
		numVerts = 0;
		pt = NULL;
		numFaces = 0;
		face = NULL;
		slideX = slideY = slideZ = 0;
		rotateX = rotateY = rotateZ = 0;
	}
	~Mesh() {
		if (pt != NULL) {
			delete[] pt;
		}
		if (face != NULL) {
			delete[] face;
		}
		numVerts = 0;
		numFaces = 0;
	}

	void DrawWireframe();
	void DrawColor();
	void Draw();
	void SetColor(int colorIdx);

	void CreateCuboid(float fSizeX, float fSizeY, float fSizeZ);
	void CreateCylinder(int nSegment, float fHeight, float fRadius);
	void CreateSphere(int nSlice, int nStack, float radius);
	void CreateUShape(float	fSizeX, float fSizeY, float	fSizeZ, float fThick);
	void CreateYPlane(int row, int col, int size);
	void CreateDodecahedron();
	void CreateIsocahedron();
	void CreateTruncatedcube(float size); 

	void CalculateFacesNorm();
};

void Mesh::DrawWireframe() {
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int f = 0; f < numFaces; f++) {
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++) {
			int iv = face[f].vert[v].vertIndex;
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::DrawColor() {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (int f = 0; f < numFaces; f++) {
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++) {
			int iv = face[f].vert[v].vertIndex;
			int ic = face[f].vert[v].colorIndex;

			//ic = f % COLORNUM;
			glColor3f(ColorArr[ic][0], ColorArr[ic][1], ColorArr[ic][2]); 
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
};

void Mesh::Draw() {
	for (int f = 0; f < numFaces; f++){
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++){
			int		iv = face[f].vert[v].vertIndex;
			glNormal3f(face[f].facenorm.x, face[f].facenorm.y, face[f].facenorm.z);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::SetColor(int colorIdx) {
	for (int f = 0; f < numFaces; f++)
	{
		for (int v = 0; v < face[f].nVerts; v++)
		{
			face[f].vert[v].colorIndex = colorIdx;
		}
	}
}

void Mesh::CreateCuboid(float fSizeX, float fSizeY, float fSizeZ)
{
	int i;

	numVerts = 8;
	pt = new Point3[numVerts];
	pt[0].set(-fSizeX,  fSizeY,  fSizeZ);
	pt[1].set( fSizeX,  fSizeY,  fSizeZ);
	pt[2].set( fSizeX,  fSizeY, -fSizeZ);
	pt[3].set(-fSizeX,  fSizeY, -fSizeZ);
	pt[4].set(-fSizeX, -fSizeY,  fSizeZ);
	pt[5].set( fSizeX, -fSizeY,  fSizeZ);
	pt[6].set( fSizeX, -fSizeY, -fSizeZ);
	pt[7].set(-fSizeX, -fSizeY, -fSizeZ);

	numFaces= 6;
	face = new Face[numFaces];

	//Left face
	face[0].nVerts = 4;
	face[0].vert = new VertexID[face[0].nVerts];
	face[0].vert[0].vertIndex = 1;
	face[0].vert[1].vertIndex = 5;
	face[0].vert[2].vertIndex = 6;
	face[0].vert[3].vertIndex = 2;
	for(i = 0; i<face[0].nVerts ; i++)
		face[0].vert[i].colorIndex = 0;
	
	//Right face
	face[1].nVerts = 4;
	face[1].vert = new VertexID[face[1].nVerts];
	face[1].vert[0].vertIndex = 0;
	face[1].vert[1].vertIndex = 3;
	face[1].vert[2].vertIndex = 7;
	face[1].vert[3].vertIndex = 4;
	for(i = 0; i<face[1].nVerts ; i++)
		face[1].vert[i].colorIndex = 1;

	//top face
	face[2].nVerts = 4;
	face[2].vert = new VertexID[face[2].nVerts];
	face[2].vert[0].vertIndex = 0;
	face[2].vert[1].vertIndex = 1;
	face[2].vert[2].vertIndex = 2;
	face[2].vert[3].vertIndex = 3;
	for(i = 0; i<face[2].nVerts ; i++)
		face[2].vert[i].colorIndex = 2;

	//bottom face
	face[3].nVerts = 4;
	face[3].vert = new VertexID[face[3].nVerts];
	face[3].vert[0].vertIndex = 7;
	face[3].vert[1].vertIndex = 6;
	face[3].vert[2].vertIndex = 5;
	face[3].vert[3].vertIndex = 4;
	for(i = 0; i<face[3].nVerts ; i++)
		face[3].vert[i].colorIndex = 3;

	//near face
	face[4].nVerts = 4;
	face[4].vert = new VertexID[face[4].nVerts];
	face[4].vert[0].vertIndex = 4;
	face[4].vert[1].vertIndex = 5;
	face[4].vert[2].vertIndex = 1;
	face[4].vert[3].vertIndex = 0;
	for(i = 0; i<face[4].nVerts ; i++)
		face[4].vert[i].colorIndex = 4;

	//Far face
	face[5].nVerts = 4;
	face[5].vert = new VertexID[face[5].nVerts];
	face[5].vert[0].vertIndex = 3;
	face[5].vert[1].vertIndex = 2;
	face[5].vert[2].vertIndex = 6;
	face[5].vert[3].vertIndex = 7;
	for(i = 0; i<face[5].nVerts ; i++)
		face[5].vert[i].colorIndex = 5;
}

void Mesh::CreateCylinder(int nSegment, float fHeight, float fRadius)
{
	numVerts = 2 * nSegment + 2;
	pt = new Point3[numVerts];
	float sg = 360.0/nSegment;

	//Mat duoi
	for(int i = 0; i < nSegment; i++){
		pt[i].set(fRadius*cos(DEG2RAD*i*sg), -fHeight/2, fRadius*sin(DEG2RAD*i*sg));
	}
	pt[nSegment].set(0, -fHeight/2, 0);

	//Mat tren
	for(int i = nSegment+1; i < 2*nSegment+1; i++){
		pt[i].set(fRadius*cos(DEG2RAD*(i-nSegment-1)*sg), fHeight/2, fRadius*sin(DEG2RAD*(i-nSegment-1)*sg));
	}
	pt[2*nSegment + 1].set(0, fHeight/2, 0); 


	numFaces = 3 * nSegment;
	face = new Face[numFaces];

	
	for(int i = 0; i < nSegment; i++){
		// Mat ngang
		face[i].nVerts = 4;
		face[i].vert = new VertexID[4];
		face[i].vert[0].vertIndex = i;
		face[i].vert[1].vertIndex = nSegment + i + 1;
		if (i != nSegment - 1){
			face[i].vert[2].vertIndex = nSegment + i + 2;
			face[i].vert[3].vertIndex =  i + 1;
		}
		else {
			face[i].vert[2].vertIndex = nSegment + 1;
			face[i].vert[3].vertIndex = 0;
		}

		// Mat duoi
		face[i+nSegment].nVerts = 3;
		face[i+nSegment].vert = new VertexID[3];
		face[i+nSegment].vert[0].vertIndex = nSegment;
		face[i+nSegment].vert[1].vertIndex = i;
		if (i != nSegment - 1)
			face[i+nSegment].vert[2].vertIndex = i + 1;
		else
			face[i+nSegment].vert[2].vertIndex = 0;

		// Mat tren
		face[i+2*nSegment].nVerts = 3;
		face[i+2*nSegment].vert = new VertexID[3];
		face[i+2*nSegment].vert[0].vertIndex = 2*nSegment + 1;
		face[i+2*nSegment].vert[2].vertIndex = i + nSegment + 1;
		if (i != nSegment -1)
			face[i+2*nSegment].vert[1].vertIndex = i + nSegment + 2 ;		
		else
			face[i+2*nSegment].vert[1].vertIndex = nSegment + 1;
	}
}

void Mesh::CreateSphere(int nSlice, int nStack, float radius)
{
	numVerts = (nSlice)*(nStack +1);
	pt = new Point3[numVerts];

	for(int v = 0; v<= nStack ; v++){
		float bankinh = radius*sin(v*PI/nStack);
		for(int u = 0; u < nSlice; u++){
			pt[u+v*nSlice].set( bankinh*cos(2*PI*u/nSlice), bankinh*sin(2*PI*u/nSlice), radius*cos(PI*v/nStack));
		}
	}

	numFaces= nSlice * nStack;
	face = new Face[numFaces];

	for(int u = 0; u < nSlice; u++){
		for(int v = 0; v< nStack ; v++){
			int nFaces = u+v*nSlice;

			face[nFaces].nVerts = 4;
			face[nFaces].vert = new VertexID[face[0].nVerts];

			if(u == nSlice -1){
				face[nFaces].vert[2].vertIndex = (v+1)*nSlice;
				face[nFaces].vert[3].vertIndex = v*nSlice;
			} else {
				face[nFaces].vert[2].vertIndex = (v+1)*nSlice + (u+1);
				face[nFaces].vert[3].vertIndex = v*nSlice + u + 1;
			}

			face[nFaces].vert[0].vertIndex = u+v*nSlice ;
			face[nFaces].vert[1].vertIndex = (v+1)*nSlice + u;
		}
	}
}

void Mesh::CreateUShape(float fSizeX, float fSizeY, float fSizeZ, float fThick)
{
	int i, j;

	numVerts = 20;
	pt = new Point3[numVerts];
	pt[0].set(-fSizeX/2,         fSizeY/2,  fSizeZ/2);
	pt[1].set(-fSizeX/2+fThick,  fSizeY/2,  fSizeZ/2);
	pt[2].set(-fSizeX/2+fThick,  fSizeY/2, -fSizeZ/2+fThick);
	pt[3].set( fSizeX/2-fThick,  fSizeY/2, -fSizeZ/2+fThick);
	pt[4].set( fSizeX/2-fThick,  fSizeY/2,  fSizeZ/2);
	pt[5].set( fSizeX/2,         fSizeY/2,  fSizeZ/2);
	pt[6].set( fSizeX/2,         fSizeY/2, -fSizeZ/2);
	pt[7].set( fSizeX/2-fThick,  fSizeY/2, -fSizeZ/2);
	pt[8].set(-fSizeX/2+fThick,  fSizeY/2, -fSizeZ/2);
	pt[9].set(-fSizeX/2,         fSizeY/2, -fSizeZ/2);

	pt[10].set(-fSizeX/2,        -fSizeY/2,  fSizeZ/2);
	pt[11].set(-fSizeX/2+fThick, -fSizeY/2,  fSizeZ/2);
	pt[12].set(-fSizeX/2+fThick, -fSizeY/2, -fSizeZ/2+fThick);
	pt[13].set( fSizeX/2-fThick, -fSizeY/2, -fSizeZ/2+fThick);
	pt[14].set( fSizeX/2-fThick, -fSizeY/2,  fSizeZ/2);
	pt[15].set( fSizeX/2,        -fSizeY/2,  fSizeZ/2);
	pt[16].set( fSizeX/2,        -fSizeY/2, -fSizeZ/2);
	pt[17].set( fSizeX/2-fThick, -fSizeY/2, -fSizeZ/2);
	pt[18].set(-fSizeX/2+fThick, -fSizeY/2, -fSizeZ/2);
	pt[19].set(-fSizeX/2,        -fSizeY/2, -fSizeZ/2);

	numFaces= 16;
	face = new Face[numFaces];

	for (i = 0; i < 3; i++) {
		// Top face
		face[i].nVerts = 4;
		face[i].vert = new VertexID[face[i].nVerts];
		face[i].vert[0].vertIndex = 2*i;
		face[i].vert[1].vertIndex = 2*i + 1;
		face[i].vert[2].vertIndex = 9 - i - 1;
		face[i].vert[3].vertIndex = 9 - i;
		for(j = 0; j<face[i].nVerts ; j++)
			face[i].vert[j].colorIndex = 0;

		// Bottom face
		face[i+3].nVerts = 4;
		face[i+3].vert = new VertexID[face[i+3].nVerts];
		face[i+3].vert[3].vertIndex = 2*i + 10;
		face[i+3].vert[2].vertIndex = 2*i + 1 + 10;
		face[i+3].vert[1].vertIndex = 9 - i - 1 + 10;
		face[i+3].vert[0].vertIndex = 9 - i + 10;
		for(j = 0; j<face[i+3].nVerts ; j++)
			face[i+3].vert[j].colorIndex = 0;
	}

	for (i = 0; i < 10; i++) {
		face[i+6].nVerts = 4;
		face[i+6].vert = new VertexID[face[i+6].nVerts];
		if (i == 0) {
			face[i+6].vert[0].vertIndex = 0;
			face[i+6].vert[1].vertIndex = 9;
			face[i+6].vert[2].vertIndex = 19;
			face[i+6].vert[3].vertIndex = 10;
		} else {
			face[i+6].vert[0].vertIndex = i;
			face[i+6].vert[1].vertIndex = i - 1;
			face[i+6].vert[2].vertIndex = i + 9;
			face[i+6].vert[3].vertIndex = i + 10;
		}
		for(j = 0; j<face[i+6].nVerts ; j++)
			face[i+6].vert[j].colorIndex = 0;
	}
}


void Mesh::CreateYPlane(int row, int col, int size)
{
	numVerts = (row + 1)*(col + 1);
	pt = new Point3[numVerts];
	for (int i = 0; i <= row; i++)
		for (int j = 0; j <= col; j++)	
			pt[i + j + col * i].set((-col/2 + j)*size, 0, (-row/2 + i)*size);

	numFaces = row*col;
	face = new Face[numFaces];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			int index = i + j + (col - 1)*i;
			int vert = i + j + col * i; // toa do dinh dau tien cua mat index
			face[index].nVerts = 4;
			face[index].vert = new VertexID[4];
			face[index].vert[0].vertIndex = vert;
			face[index].vert[1].vertIndex = i + j + 1 + col*(i + 1);
			face[index].vert[2].vertIndex = i + j + 2 + col*(i + 1);
			face[index].vert[3].vertIndex = i + j + col*i + 1;

			for (int i = 0; i < 4; i++)
				if (vert % 2 == 0)
					face[index].vert[i].colorIndex = 12;
				else
					face[index].vert[i].colorIndex = 13;
		}
	}
}

void Mesh::CreateIsocahedron() 
{
	numVerts = 12;
	float t = (sqrt(5.0)-1)/2;
	pt = new Point3[numVerts];
	pt[0].set(0, 1, t);
	pt[1].set(0, 1, -t);
	pt[2].set(1, t, 0);
	pt[3].set(1, -t, 0);
	pt[4].set(0, -1, -t);
	pt[5].set(0, -1, t);
	pt[6].set(t, 0, 1);
	pt[7].set(-t, 0, 1);
	pt[8].set(t, 0, -1);
	pt[9].set(-t, 0, -1);
	pt[10].set(-1, t, 0);
	pt[11].set(-1, -t, 0);

	numFaces = 20;
	face = new Face[numFaces];
	face[0].nVerts = 3;
	face[0].vert = new VertexID[3];
	face[0].vert[0].vertIndex = 0;
	face[0].vert[1].vertIndex = 2;
	face[0].vert[2].vertIndex = 6;

	face[1].nVerts = 3;
	face[1].vert = new VertexID[3];
	face[1].vert[0].vertIndex = 6;
	face[1].vert[1].vertIndex = 2;
	face[1].vert[2].vertIndex = 3;

	face[2].nVerts = 3;
	face[2].vert = new VertexID[3];
	face[2].vert[0].vertIndex = 6;
	face[2].vert[1].vertIndex = 3;
	face[2].vert[2].vertIndex = 5;

	face[3].nVerts = 3;
	face[3].vert = new VertexID[3];
	face[3].vert[0].vertIndex = 7;
	face[3].vert[1].vertIndex = 6;
	face[3].vert[2].vertIndex = 5;

	face[4].nVerts = 3;
	face[4].vert = new VertexID[3];
	face[4].vert[0].vertIndex = 7;
	face[4].vert[1].vertIndex = 0;
	face[4].vert[2].vertIndex = 6;

	face[5].nVerts = 3;
	face[5].vert = new VertexID[3];
	face[5].vert[0].vertIndex = 2;
	face[5].vert[1].vertIndex = 8;
	face[5].vert[2].vertIndex = 3;

	face[6].nVerts = 3;
	face[6].vert = new VertexID[3];
	face[6].vert[0].vertIndex = 1;
	face[6].vert[1].vertIndex = 8;
	face[6].vert[2].vertIndex = 2;

	face[7].nVerts = 3;
	face[7].vert = new VertexID[3];
	face[7].vert[0].vertIndex = 0;
	face[7].vert[1].vertIndex = 1;
	face[7].vert[2].vertIndex = 2;

	face[8].nVerts = 3;
	face[8].vert = new VertexID[3];
	face[8].vert[0].vertIndex = 10;
	face[8].vert[1].vertIndex = 1;
	face[8].vert[2].vertIndex = 0;

	face[9].nVerts = 3;
	face[9].vert = new VertexID[3];
	face[9].vert[0].vertIndex = 10;
	face[9].vert[1].vertIndex = 9;
	face[9].vert[2].vertIndex = 1;

	face[10].nVerts = 3;
	face[10].vert = new VertexID[3];
	face[10].vert[0].vertIndex = 1;
	face[10].vert[1].vertIndex = 9;
	face[10].vert[2].vertIndex = 8;

	face[11].nVerts = 3;
	face[11].vert = new VertexID[3];
	face[11].vert[0].vertIndex = 3;
	face[11].vert[1].vertIndex = 8;
	face[11].vert[2].vertIndex = 4;

	face[12].nVerts = 3;
	face[12].vert = new VertexID[3];
	face[12].vert[0].vertIndex = 5;
	face[12].vert[1].vertIndex = 3;
	face[12].vert[2].vertIndex = 4;

	face[13].nVerts = 3;
	face[13].vert = new VertexID[3];
	face[13].vert[0].vertIndex = 11;
	face[13].vert[1].vertIndex = 5;
	face[13].vert[2].vertIndex = 4;

	face[14].nVerts = 3;
	face[14].vert = new VertexID[3];
	face[14].vert[0].vertIndex = 11;
	face[14].vert[1].vertIndex = 10;
	face[14].vert[2].vertIndex = 7;

	face[15].nVerts = 3;
	face[15].vert = new VertexID[3];
	face[15].vert[0].vertIndex = 7;
	face[15].vert[1].vertIndex = 10;
	face[15].vert[2].vertIndex = 0;

	face[16].nVerts = 3;
	face[16].vert = new VertexID[3];
	face[16].vert[0].vertIndex = 9;
	face[16].vert[1].vertIndex = 11;
	face[16].vert[2].vertIndex = 4;

	face[17].nVerts = 3;
	face[17].vert = new VertexID[3];
	face[17].vert[0].vertIndex = 9;
	face[17].vert[1].vertIndex = 4;
	face[17].vert[2].vertIndex = 8;

	face[18].nVerts = 3;
	face[18].vert = new VertexID[3];
	face[18].vert[0].vertIndex = 11;
	face[18].vert[1].vertIndex = 7;
	face[18].vert[2].vertIndex = 5;

	face[19].nVerts = 3;
	face[19].vert = new VertexID[3];
	face[19].vert[0].vertIndex = 11;
	face[19].vert[1].vertIndex = 9;
	face[19].vert[2].vertIndex = 10;

	int color = 1;
	for (int i = 0; i < numFaces; i++) {
		for (int j = 0; j < 3; j++)
			face[i].vert[j].colorIndex = color;
		color = (color+1)%14;
	}
}

void Mesh::CreateDodecahedron()
{	
	int i;

	numVerts = 20;
	pt = new Point3[numVerts];
	Mesh iso;
	iso.CreateIsocahedron();
	for (i = 0; i < numVerts; i++) {
		float x = ( iso.pt[iso.face[i].vert[0].vertIndex].x + 
					iso.pt[iso.face[i].vert[1].vertIndex].x +
					iso.pt[iso.face[i].vert[2].vertIndex].x ) / 3.0;
		float y = ( iso.pt[iso.face[i].vert[0].vertIndex].y + 
					iso.pt[iso.face[i].vert[1].vertIndex].y +
					iso.pt[iso.face[i].vert[2].vertIndex].y ) / 3.0;
		float z = ( iso.pt[iso.face[i].vert[0].vertIndex].z + 
					iso.pt[iso.face[i].vert[1].vertIndex].z +
					iso.pt[iso.face[i].vert[2].vertIndex].z ) / 3.0;
		pt[i].set(x, y, z);
	}

	numFaces = 12;
	face = new Face[numFaces];
	face[0].nVerts = 5;
	face[0].vert = new VertexID[5];
	face[0].vert[0].vertIndex = 0;
	face[0].vert[1].vertIndex = 4;
	face[0].vert[2].vertIndex = 15;
	face[0].vert[3].vertIndex = 8;
	face[0].vert[4].vertIndex = 7;

	face[1].nVerts = 5;
	face[1].vert = new VertexID[5];
	face[1].vert[0].vertIndex = 6;
	face[1].vert[1].vertIndex = 7;
	face[1].vert[2].vertIndex = 8;
	face[1].vert[3].vertIndex = 9;
	face[1].vert[4].vertIndex = 10;

	face[2].nVerts = 5;
	face[2].vert = new VertexID[5];
	face[2].vert[0].vertIndex = 5;
	face[2].vert[1].vertIndex = 1;
	face[2].vert[2].vertIndex = 0;
	face[2].vert[3].vertIndex = 7;
	face[2].vert[4].vertIndex = 6;

	face[3].nVerts = 5;
	face[3].vert = new VertexID[5];
	face[3].vert[0].vertIndex = 12;
	face[3].vert[1].vertIndex = 2;
	face[3].vert[2].vertIndex = 1;
	face[3].vert[3].vertIndex = 5;
	face[3].vert[4].vertIndex = 11;

	face[4].nVerts = 5;
	face[4].vert = new VertexID[5];
	face[4].vert[0].vertIndex = 13;
	face[4].vert[1].vertIndex = 12;
	face[4].vert[2].vertIndex = 11;
	face[4].vert[3].vertIndex = 17;
	face[4].vert[4].vertIndex = 16;

	face[5].nVerts = 5;
	face[5].vert = new VertexID[5];
	face[5].vert[0].vertIndex = 3;
	face[5].vert[1].vertIndex = 2;
	face[5].vert[2].vertIndex = 12;
	face[5].vert[3].vertIndex = 13;
	face[5].vert[4].vertIndex = 18;

	face[6].nVerts = 5;
	face[6].vert = new VertexID[5];
	face[6].vert[0].vertIndex = 0;
	face[6].vert[1].vertIndex = 1;
	face[6].vert[2].vertIndex = 2;
	face[6].vert[3].vertIndex = 3;
	face[6].vert[4].vertIndex = 4;

	face[7].nVerts = 5;
	face[7].vert = new VertexID[5];
	face[7].vert[0].vertIndex = 4;
	face[7].vert[1].vertIndex = 3;
	face[7].vert[2].vertIndex = 18;
	face[7].vert[3].vertIndex = 14;
	face[7].vert[4].vertIndex = 15;

	face[8].nVerts = 5;
	face[8].vert = new VertexID[5];
	face[8].vert[0].vertIndex = 11;
	face[8].vert[1].vertIndex = 5;
	face[8].vert[2].vertIndex = 6;
	face[8].vert[3].vertIndex = 10;
	face[8].vert[4].vertIndex = 17;

	face[9].nVerts = 5;
	face[9].vert = new VertexID[5];
	face[9].vert[0].vertIndex = 16;
	face[9].vert[1].vertIndex = 17;
	face[9].vert[2].vertIndex = 10;
	face[9].vert[3].vertIndex = 9;
	face[9].vert[4].vertIndex = 19;

	face[10].nVerts = 5;
	face[10].vert = new VertexID[5];
	face[10].vert[0].vertIndex = 15;
	face[10].vert[1].vertIndex = 14;
	face[10].vert[2].vertIndex = 19;
	face[10].vert[3].vertIndex = 9;
	face[10].vert[4].vertIndex = 8;

	face[11].nVerts = 5;
	face[11].vert = new VertexID[5];
	face[11].vert[0].vertIndex = 18;
	face[11].vert[1].vertIndex = 13;
	face[11].vert[2].vertIndex = 16;
	face[11].vert[3].vertIndex = 19;
	face[11].vert[4].vertIndex = 14;

	int color = 1;
	for (i = 0; i < numFaces; i++) {
		for (int j = 0; j < 5; j++)
			face[i].vert[j].colorIndex = color;
		color = (color+1)%14;
	}
}

void Mesh::CreateTruncatedcube(float size) {
	numVerts = 24;
	pt = new Point3[numVerts];
	pt[0].set(-size,   size,  size/2);
	pt[1].set(-size/2, size,  size);
	pt[2].set( size/2, size,  size);
	pt[3].set( size,   size,  size/2);
	pt[4].set( size,   size, -size/2);
	pt[5].set( size/2, size, -size);
	pt[6].set(-size/2, size, -size);
	pt[7].set(-size,   size, -size/2);

	pt[8].set( -size, size/2,  size);
	pt[9].set(  size, size/2,  size);
	pt[10].set( size, size/2, -size);
	pt[11].set(-size, size/2, -size);

	pt[12].set(-size, -size/2,  size);
	pt[13].set( size, -size/2,  size);
	pt[14].set( size, -size/2, -size);
	pt[15].set(-size, -size/2, -size);

	pt[16].set(-size,   -size,  size/2);
	pt[17].set(-size/2, -size,  size);
	pt[18].set( size/2, -size,  size);
	pt[19].set( size,   -size,  size/2);
	pt[20].set( size,   -size, -size/2);
	pt[21].set( size/2, -size, -size);
	pt[22].set(-size/2, -size, -size);
	pt[23].set(-size,   -size, -size/2);

	numFaces = 14;
	face = new Face[numFaces];
	face[0].nVerts = 8;
	face[0].vert = new VertexID[8];
	face[0].vert[0].vertIndex = 0;
	face[0].vert[1].vertIndex = 1;
	face[0].vert[2].vertIndex = 2;
	face[0].vert[3].vertIndex = 3;
	face[0].vert[4].vertIndex = 4;
	face[0].vert[5].vertIndex = 5;
	face[0].vert[6].vertIndex = 6;
	face[0].vert[7].vertIndex = 7;

	face[1].nVerts = 8;
	face[1].vert = new VertexID[8];
	face[1].vert[0].vertIndex = 1;
	face[1].vert[1].vertIndex = 8;
	face[1].vert[2].vertIndex = 12;
	face[1].vert[3].vertIndex = 17;
	face[1].vert[4].vertIndex = 18;
	face[1].vert[5].vertIndex = 13;
	face[1].vert[6].vertIndex = 9;
	face[1].vert[7].vertIndex = 2;

	face[2].nVerts = 8;
	face[2].vert = new VertexID[8];
	face[2].vert[0].vertIndex = 23;
	face[2].vert[1].vertIndex = 22;
	face[2].vert[2].vertIndex = 21;
	face[2].vert[3].vertIndex = 20;
	face[2].vert[4].vertIndex = 19;
	face[2].vert[5].vertIndex = 18;
	face[2].vert[6].vertIndex = 17;
	face[2].vert[7].vertIndex = 16;

	face[3].nVerts = 8;
	face[3].vert = new VertexID[8];
	face[3].vert[0].vertIndex = 22;
	face[3].vert[1].vertIndex = 15;
	face[3].vert[2].vertIndex = 11;
	face[3].vert[3].vertIndex = 6;
	face[3].vert[4].vertIndex = 5;
	face[3].vert[5].vertIndex = 10;
	face[3].vert[6].vertIndex = 14;
	face[3].vert[7].vertIndex = 21;

	face[4].nVerts = 8;
	face[4].vert = new VertexID[8];
	face[4].vert[0].vertIndex = 0;
	face[4].vert[1].vertIndex = 7;
	face[4].vert[2].vertIndex = 11;
	face[4].vert[3].vertIndex = 15;
	face[4].vert[4].vertIndex = 23;
	face[4].vert[5].vertIndex = 16;
	face[4].vert[6].vertIndex = 12;
	face[4].vert[7].vertIndex = 8;

	face[5].nVerts = 8;
	face[5].vert = new VertexID[8];
	face[5].vert[0].vertIndex = 4;
	face[5].vert[1].vertIndex = 3;
	face[5].vert[2].vertIndex = 9;
	face[5].vert[3].vertIndex = 13;
	face[5].vert[4].vertIndex = 19;
	face[5].vert[5].vertIndex = 20;
	face[5].vert[6].vertIndex = 14;
	face[5].vert[7].vertIndex = 10;

	face[6].nVerts = 3;
	face[6].vert = new VertexID[3];
	face[6].vert[0].vertIndex = 1;
	face[6].vert[1].vertIndex = 0;
	face[6].vert[2].vertIndex = 8;

	face[7].nVerts = 3;
	face[7].vert = new VertexID[3];
	face[7].vert[0].vertIndex = 3;
	face[7].vert[1].vertIndex = 2;
	face[7].vert[2].vertIndex = 9;

	face[8].nVerts = 3;
	face[8].vert = new VertexID[3];
	face[8].vert[0].vertIndex = 5;
	face[8].vert[1].vertIndex = 4;
	face[8].vert[2].vertIndex = 10;

	face[9].nVerts = 3;
	face[9].vert = new VertexID[3];
	face[9].vert[0].vertIndex = 7;
	face[9].vert[1].vertIndex = 6;
	face[9].vert[2].vertIndex = 11;

	face[10].nVerts = 3;
	face[10].vert = new VertexID[3];
	face[10].vert[0].vertIndex = 16;
	face[10].vert[1].vertIndex = 17;
	face[10].vert[2].vertIndex = 12;

	face[11].nVerts = 3;
	face[11].vert = new VertexID[3];
	face[11].vert[0].vertIndex = 18;
	face[11].vert[1].vertIndex = 19;
	face[11].vert[2].vertIndex = 13;

	face[12].nVerts = 3;
	face[12].vert = new VertexID[3];
	face[12].vert[0].vertIndex = 20;
	face[12].vert[1].vertIndex = 21;
	face[12].vert[2].vertIndex = 14;

	face[13].nVerts = 3;
	face[13].vert = new VertexID[3];
	face[13].vert[0].vertIndex = 22;
	face[13].vert[1].vertIndex = 23;
	face[13].vert[2].vertIndex = 15;

	int color = 1;
	for (int i = 0; i < numFaces; i++) {
		for (int j = 0; j < face[i].nVerts; j++)
			face[i].vert[j].colorIndex = color;
		color = (color+1)%14;
	}
}

void Mesh::CalculateFacesNorm() {
	for(int i = 0; i<numFaces; i++) {
		int n = face[i].nVerts;
		float x=0, y=0, z=0;
	
		for(int j = 0; j < n-1; j++){
			x+=(pt[face[i].vert[j].vertIndex].y-pt[face[i].vert[j+1].vertIndex].y)*(pt[face[i].vert[j].vertIndex].z+pt[face[i].vert[j+1].vertIndex].z);
			y+=(pt[face[i].vert[j].vertIndex].z-pt[face[i].vert[j+1].vertIndex].z)*(pt[face[i].vert[j].vertIndex].x+pt[face[i].vert[j+1].vertIndex].x);
			z+=(pt[face[i].vert[j].vertIndex].x-pt[face[i].vert[j+1].vertIndex].x)*(pt[face[i].vert[j].vertIndex].y+pt[face[i].vert[j+1].vertIndex].y);
		}
		x+=(pt[face[i].vert[n-1].vertIndex].y-pt[face[i].vert[0].vertIndex].y)*(pt[face[i].vert[n-1].vertIndex].z+pt[face[i].vert[0].vertIndex].z);
		y+=(pt[face[i].vert[n-1].vertIndex].z-pt[face[i].vert[0].vertIndex].z)*(pt[face[i].vert[n-1].vertIndex].x+pt[face[i].vert[0].vertIndex].x);
		z+=(pt[face[i].vert[n-1].vertIndex].x-pt[face[i].vert[0].vertIndex].x)*(pt[face[i].vert[n-1].vertIndex].y+pt[face[i].vert[0].vertIndex].y);

		face[i].facenorm.set(x, y, z);
		face[i].facenorm.normalize();
	}

}

int		screenWidth = 600;
int		screenHeight= 600;

bool	bWireFrame = false;
bool	b4View = false;
bool	isCartoon = false;
bool	light1 = false;

float	YPlanePos = 0;

float	cylBaseRadius = 1.2;
float	cylBaseHeight = 0.2; 
float	cylBaseRotateAngle = 0.0;

float	cuboidBaseSizeXZ = 0.3;
float	cuboidBaseSizeY = 3.0;
float	cuboidRotateAngle = 0;

float	cylAxisRadius = 0.1;
float	cylAxisHeight = 0.3; 
float	cylAxisOffset = 0.2;
float	cylAxisRotateAngle = 0;

Mesh	planeY;
Mesh	cylBase;
Mesh	cuboidBase;
Mesh    cylAxis;
Mesh	uShape;
Mesh	cubAxis;
Mesh	cylMain;
Mesh	sphere;

Mesh	iso;
Mesh	dode;
Mesh	sphe;
Mesh	trun;

float camera_angle;
float camera_height;
float camera_dis;
float camera_X, camera_Y, camera_Z;
float lookAt_X, lookAt_Y, lookAt_Z;

GLfloat	lightDiffuse[]={1.0f, 1.0f, 1.0f, 1.0f};
GLfloat	lightSpecular[]={1.0f, 1.0f, 1.0f, 1.0f};
GLfloat	lightAmbient[]={0.4f, 0.4f, 0.4f, 1.0f};
GLfloat light_position1[]={6.0f, 6.0f, 6.0f, 0.0f};
GLfloat light_position2[] ={-6.0f, 6.0f, -6.0f, 0.0f};

GLfloat ambient[] =			{0.0, 0.0, 0.0, 1.0};
GLfloat specular[] =		{1.0, 1.0, 1.0, 1.0};
GLfloat diffuseRed[] =		{1.0, 0.0, 0.0, 1.0};
GLfloat diffuseGreen[] =	{0.0, 1.0, 0.0, 1.0};
GLfloat diffuseBlue[] =		{0.0, 0.0, 1.0, 1.0};
GLfloat diffuseYellow[] =	{1.0, 1.0, 0.0, 1.0};
GLfloat diffuseMagenta[] =	{1.0, 0.0, 1.0, 1.0};
GLfloat diffuseCyan[] =		{0.0, 1.0, 1.0, 1.0};
GLfloat diffuseOrange[] =	{1.0, 0.65, 0.0, 1.0};
GLfloat diffuseWGray[] =	{0.7, 0.7, 0.7, 1.0};
GLfloat diffuseGray[] =		{0.4, 0.4, 0.4, 1.0};
GLfloat shiness = 100;

GLfloat floorPlane[4];
GLfloat floorShadow1[4][4];
GLfloat floorShadow2[4][4];
GLfloat floorVertices[4][3] = {
  { -5.0, 0.0, 5.0 },
  {  5.0, 0.0, 5.0 },
  {  5.0, 0.0, -5.0 },
  { -5.0, 0.0, -5.0 },
};

bool pickDode, pickIso, pickSphe, pickTrun; 
bool cred, cgreen, cblue, cpink;
int tmpDode =	0;
int tmpIso  =	0;
int tmpSphe =	0;
int tmpTrun =	0;


void mySetupCameraVolume(int nType)
{
	if(nType == 4)
	{
		glMatrixMode(GL_PROJECTION);			// set projection matrix current matrix
		glLoadIdentity();						// reset projection matrix

		// calculate aspect ratio of window
		gluPerspective(60.0f,(GLfloat)screenWidth/(GLfloat)screenHeight,1.0f,1000.0f);
	}
	else 
	{
		glMatrixMode(GL_PROJECTION);			// set projection matrix current matrix
		glLoadIdentity();						// reset projection matrix
		glOrtho(-5, 5, -5, 5, -1000, 1000);
	}
	
}

void changeCameraPos()
{
	camera_Y = camera_height;
	camera_X = camera_dis * cos(camera_angle);
	camera_Z = camera_dis * sin(camera_angle);
}


void mySpecialKeyboard(int theKey, int mouseX, int mouseY)
{
	switch(theKey) {

	case GLUT_KEY_UP:
		camera_height += 0.1;
		changeCameraPos();
		break;
	
	case GLUT_KEY_DOWN:
		camera_height -= 0.1;
		changeCameraPos();
		break;

	case GLUT_KEY_RIGHT:
		camera_angle += 0.1;
		changeCameraPos();
		break;

	case GLUT_KEY_LEFT:
		camera_angle -= 0.1;
		changeCameraPos();
		break;
	
	default:
		break;
	}

	glutPostRedisplay();
}

void myKeyboard(unsigned char key, int x, int y)
{
	float	fRInc;
	float	fAngle;
    switch(key)
    {
	case '1':
		cylBaseRotateAngle += 2;
		if(cylBaseRotateAngle > 360)
			cylBaseRotateAngle -= 360;
		break;
	case '2':
		cylBaseRotateAngle -= 2;
		if(cylBaseRotateAngle < 0)
			cylBaseRotateAngle += 360;
		break;
	case '3':
		cylAxisRotateAngle += 2;
		if(cylAxisRotateAngle > 360)
			cylAxisRotateAngle -= 360;
		break;
	case '4':
		cylAxisRotateAngle -= 2;
		if(cylAxisRotateAngle < 0)
			cylAxisRotateAngle += 360;
		break;
	case '5':
		cuboidRotateAngle += 2;
		if (cuboidRotateAngle > 360)
			cuboidRotateAngle -= 360;
		break;
	case '6':
		cuboidRotateAngle -= 2;
		if (cuboidRotateAngle < 0)
			cuboidRotateAngle += 360;
		break;
	case 'w':
	case 'W':
		bWireFrame = !bWireFrame;
		break;
	case 'a':
	case 'A':
		isCartoon = !isCartoon;
		break;
	case 'v':
	case 'V':
		b4View = !b4View;
		break;
	case '+' :
		camera_dis += 0.1;
		changeCameraPos();
		break;
	case '-' : 
		camera_dis -= 0.1;
		changeCameraPos();
		break;
	case 'l':
	case 'L':
		light1 = !light1;
		break;
	case 'r':
	case 'R':
		cred=true;
		cgreen=false;
		cblue=false;
		cpink= false;
		break;
	case 'g':
	case 'G':
		cred=false;
		cgreen=true;
		cblue=false;
		cpink = false;
		break;
	case 'b':
	case 'B':
		cred=false;
		cgreen=false;
		cblue=true;
		cpink=false;
		break;
	case 'p':
	case 'P':
		cred=false;
		cgreen=false;
		cblue=false;
		cpink=true;
		break;
	}
    glutPostRedisplay();
}

void setupMaterial(float ambient[], float diffuse[], float specular[], float shiness)
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shiness);
}

bool cylBaseRotate = true;
bool cylAxisRotate = true;
bool cuboidRotate = true;
void cartoonAction(int val) 
{
	if (isCartoon) {
		if (cylBaseRotate)
			cylBaseRotateAngle += val;
		else
			cylBaseRotateAngle -= val;
		if (cylAxisRotate)
			cylAxisRotateAngle += 2*val;
		else
			cylAxisRotateAngle -= 2*val;
		if (cuboidRotate)
			cuboidRotateAngle += 3*val;
		else
			cuboidRotateAngle -= 3*val;

		if (cylBaseRotateAngle > 360 || cylBaseRotateAngle < 0)
			cylBaseRotate = !cylBaseRotate;
		if (cylAxisRotateAngle > 360 || cylAxisRotateAngle < 0)
			cylAxisRotate = !cylAxisRotate;
		if (cuboidRotateAngle > 360 || cuboidRotateAngle < 0)
			cuboidRotate = !cuboidRotate;
	}
	glutPostRedisplay();
	glutTimerFunc(30, cartoonAction, 2);
}

Point3 nearPoint(Point3 A, Point3 B, Point3 C) 
{ 
	Vector3 AB;
	AB.set(B.x - A.x, B.y-A.y, B.z-A.z); 
	Vector3 AC;
	AC.set(C.x - A.x, C.y-A.y, C.z-A.z);

	float t = AC.dot(AB) / AB.dot(AB); 
	Point3 Q;
	Q.set(A.x + AB.x*t, A.y + AB.y*t, A.z + AB.z*t);

	return Q; 
}  

void leftClick(int button, int state, int x, int y){

	switch(button){
	case GLUT_LEFT_BUTTON:
		cred=false;
		cblue=false;
		cgreen=false;
		cpink=false;
		Point3 start, end;
		Point3 arrayPoint[4];
		arrayPoint[0].set( 3, 0.4,  3); //dode
		arrayPoint[1].set(-3, 0.4,  3); //iso
		arrayPoint[2].set(-3, 0.4, -3); //sphe
		arrayPoint[3].set( 3, 0.4, -3); //trun
		double k, l, m;
		double matModelView[16], matProjection[16]; 
		int viewport[4]; 
		glGetDoublev( GL_MODELVIEW_MATRIX, matModelView ); 
		glGetDoublev( GL_PROJECTION_MATRIX, matProjection ); 
		glGetIntegerv( GL_VIEWPORT, viewport ); 
		double winX = (double)x; 
		double winY = viewport[3] - (double)y; 
		gluUnProject(winX, winY, 0.0, matModelView, matProjection, viewport, &k, &l, &m); 
		start.set(k, l, m);
		gluUnProject(winX, winY, 1.0, matModelView, matProjection, viewport, &k, &l, &m);
		end.set(k, l, m);
		for( int i =0; i<4; i++){
			Point3 tmp;

			tmp=nearPoint(start, end, arrayPoint[i]);

			Vector3 v_distance;
			v_distance.set(tmp.x-arrayPoint[i].x, tmp.y-arrayPoint[i].y, tmp.z-arrayPoint[i].z);
			float distance= sqrt(v_distance.x*v_distance.x + v_distance.y*v_distance.y + v_distance.z*v_distance.z);
			if(distance<=0.4 &&  i==0){ //dode
				pickDode=true;
				pickIso=false;
				pickSphe=false;
				pickTrun=false;
				break;
			}
			else if(distance<=0.4 &&  i==1){ //iso
				pickDode=false;
				pickIso=true;
				pickSphe=false;
				pickTrun=false;
				break;

			}
			else if(distance<=0.4 &&  i==2){ //sphe
				pickDode=false;
				pickIso=false;
				pickSphe=true;
				pickTrun=false;
				break;
			}
			else if(distance<=0.4 &&  i==3){ //trun
				pickDode=false;
				pickIso=false;
				pickSphe=false;
				pickTrun=true;
				break;
			}
		}
		break;
	}
	glutPostRedisplay();
}
 

typedef	struct									
{
	GLubyte	* imageData;									// Image Data (Up To 32 Bits)
	GLuint	bpp;											// Image Color Depth In Bits Per Pixel
	GLuint	width;											// Image Width
	GLuint	height;											// Image Height
	GLuint	texID;											// Texture ID Used To Select A Texture
	GLuint	type;											// Image Type (GL_RGB, GL_RGBA)
} Texture;	

typedef struct
{
	GLubyte Header[12];									// TGA File Header
} TGAHeader;


typedef struct
{
	GLubyte		header[6];								// First 6 Useful Bytes From The Header
	GLuint		bytesPerPixel;							// Holds Number Of Bytes Per Pixel Used In The TGA File
	GLuint		imageSize;								// Used To Store The Image Size When Setting Aside Ram
	GLuint		temp;									// Temporary Variable
	GLuint		type;	
	GLuint		Height;									//Height of Image
	GLuint		Width;									//Width ofImage
	GLuint		Bpp;									// Bits Per Pixel
} TGA;


TGAHeader tgaheader;									// TGA header
TGA tga;										// TGA image data



GLubyte uTGAcompare[12] = {0,0,2, 0,0,0,0,0,0,0,0,0};	// Uncompressed TGA Header


bool LoadTGA(Texture * texture, char * filename)				// Load a TGA file
{
	FILE * fTGA;												// File pointer to texture file
	fTGA = fopen(filename, "rb");								// Open file for reading

	if(fTGA == NULL)											// If it didn't open....
	{
		return false;														// Exit function
	}

	if(fread(&tgaheader, sizeof(TGAHeader), 1, fTGA) == 0)					// Attempt to read 12 byte header from file
	{
		if(fTGA != NULL)													// Check to seeiffile is still open
		{
			fclose(fTGA);													// If it is, close it
		}
		return false;														// Exit function
	}

																	// an Uncompressed TGA image
	if(fread(tga.header, sizeof(tga.header), 1, fTGA) == 0)					// Read TGA header
	{										
		if(fTGA != NULL)													// if file is still open
		{
			fclose(fTGA);													// Close it
		}
		return false;														// Return failular
	}	

	texture->width  = tga.header[1] * 256 + tga.header[0];					// Determine The TGA Width	(highbyte*256+lowbyte)
	texture->height = tga.header[3] * 256 + tga.header[2];					// Determine The TGA Height	(highbyte*256+lowbyte)
	texture->bpp	= tga.header[4];										// Determine the bits per pixel
	tga.Width		= texture->width;										// Copy width into local structure						
	tga.Height		= texture->height;										// Copy height into local structure
	tga.Bpp			= texture->bpp;											// Copy BPP into local structure

	if((texture->width <= 0) || (texture->height <= 0) || ((texture->bpp != 24) && (texture->bpp !=32)))	// Make sure all information is valid
	{
		if(fTGA != NULL)													// Check if file is still open
		{
			fclose(fTGA);													// If so, close it
		}
		return false;														// Return failed
	}

	if(texture->bpp == 24)													// If the BPP of the image is 24...
		texture->type	= GL_RGB;											// Set Image type to GL_RGB
	else																	// Else if its 32 BPP
		texture->type	= GL_RGBA;											// Set image type to GL_RGBA

	tga.bytesPerPixel	= (tga.Bpp / 8);									// Compute the number of BYTES per pixel
	tga.imageSize		= (tga.bytesPerPixel * tga.Width * tga.Height);		// Compute the total amout ofmemory needed to store data
	texture->imageData	= (GLubyte *)malloc(tga.imageSize);					// Allocate that much memory

	if(texture->imageData == NULL)											// If no space was allocated
	{
		fclose(fTGA);														// Close the file
		return false;														// Return failed
	}

	if(fread(texture->imageData, 1, tga.imageSize, fTGA) != tga.imageSize)	// Attempt to read image data
	{
		if(texture->imageData != NULL)										// If imagedata has data in it
		{
			free(texture->imageData);										// Delete data from memory
		}
		fclose(fTGA);														// Close file
		return false;														// Return failed
	}

	// switch R and B
	for (int i = 0; i < tga.imageSize; i += tga.bytesPerPixel)
	{
		GLubyte temp = texture->imageData[i];
		texture->imageData[i] = texture->imageData[i+2];
		texture->imageData[i+2] = temp;
	}
	
		
	fclose(fTGA);															// Close file
	return true;															// All went well, continue on
}
Texture floorTex;

void loadTextures(void)	{
	bool status = LoadTGA(&floorTex, "marble.tga" );
	if(status) {
		glGenTextures( 1, &floorTex.texID );
		glBindTexture( GL_TEXTURE_2D, floorTex.texID );

		glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_LINEAR );
		glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_LINEAR );

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, floorTex.width, 
			floorTex.height , 0, 
			GL_RGB, GL_UNSIGNED_BYTE, floorTex.imageData);

		if(floorTex.imageData)
			free(floorTex.imageData);
	}
}

void drawFloor()
{
	glPushMatrix();
//	setupMaterial(ambient, diffuseGray, specular, shiness);
	glTranslated(planeY.slideX, planeY.slideY, planeY.slideZ);
	glDisable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D,floorTex.texID);
	glColor4f(1, 1, 1, 1.0);
	loadTextures();
	if (bWireFrame)
	{
		planeY.DrawWireframe();
	}
	else
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4f(1.0, 1.0, 1.0, 0.7);
		glBegin(GL_QUADS);

		glTexCoord2f(0,1);
		glVertex3fv(floorVertices[0]);
		glTexCoord2f(1,1);
		glVertex3fv(floorVertices[1]);
		glTexCoord2f(1,0);
		glVertex3fv(floorVertices[2]);
		glTexCoord2f(0,0);
		glVertex3fv(floorVertices[3]);
		glEnd();
		glDisable(GL_BLEND);
	}

	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);

		/*if(bWireFrame)
			planeY.DrawWireframe();
		else
			planeY.Draw();*/
	glPopMatrix();
}
void drawCylinderBase()
{
	glPushMatrix();
	setupMaterial(ambient, diffuseRed, specular, shiness);
		glRotated(cylBaseRotateAngle, 0, 1, 0);
		//////////////////////////////////
		glTranslated(0, cylBase.slideY, 0);
		if(bWireFrame)
			cylBase.DrawWireframe();
		else
			cylBase.Draw();
	glPopMatrix();
}
void drawCuboidBase()
{
	glPushMatrix();	
	setupMaterial(ambient, diffuseBlue, specular, shiness);
		glRotated(cylBaseRotateAngle, 0, 1, 0);
		//////////////////////////////////
		glTranslated(0, cuboidBase.slideY, 0);
		if(bWireFrame)
			cuboidBase.DrawWireframe();
		else
			cuboidBase.Draw();
	glPopMatrix();
}
void drawCylAxis()
{
	glPushMatrix();
		setupMaterial(ambient, diffuseGreen, specular, shiness);
		glRotated(cylBaseRotateAngle, 0, 1, 0);
		//////////////////////////////////
		glTranslated(0, cylAxis.slideY , cylAxis.slideZ);
		glRotatef(90, 1, 0, 0);
		glRotated(cylAxisRotateAngle, 0, 1, 0);
		if(bWireFrame)
			cylAxis.DrawWireframe();
		else
			cylAxis.Draw();
	glPopMatrix();
}
void drawUShape()
{
	glPushMatrix();
	setupMaterial(ambient, diffuseYellow, specular, shiness);
	glRotated(cylBaseRotateAngle, 0, 1, 0);
	//////////////////////////////////
	glTranslated(0, uShape.slideY, uShape.slideZ);
	glRotated(cylAxisRotateAngle, 0, 0, 1);
	if (bWireFrame)
		uShape.DrawWireframe();
	else
		uShape.Draw();
	glPopMatrix();

}
void drawCubAxis()
{
	glPushMatrix();
	setupMaterial(ambient, diffuseOrange, specular, shiness);
	glRotated(cylBaseRotateAngle, 0, 1, 0);
	//////////////////////////////////
	glTranslated(0, cubAxis.slideY, cubAxis.slideZ);
	glRotated(cylAxisRotateAngle, 0, 0, 1);
	glRotated(cuboidRotateAngle, 1, 0, 0);
	if (bWireFrame)
		cubAxis.DrawWireframe();
	else
		cubAxis.Draw();
	glPopMatrix();
}
void drawCylMain()
{
	glPushMatrix();
	setupMaterial(ambient, diffuseWGray, specular, shiness);
	glRotated(cylBaseRotateAngle, 0, 1, 0);
	//////////////////////////////////
	glTranslated(0, cylMain.slideY, cylMain.slideZ);
	glRotatef(90, 0, 0, 1);
	
	glRotated(cylAxisRotateAngle, 0, 0, 1);
	glRotated(cuboidRotateAngle, 0, 1, 0);
	if (bWireFrame)
		cylMain.DrawWireframe();
	else
		cylMain.Draw();
	glPopMatrix();
}
void drawSphere(int number)
{
	glPushMatrix();
	if (number == 0) setupMaterial(ambient, diffuseRed, specular, shiness);
	else if (number == 1) setupMaterial(ambient, diffuseGreen, specular, shiness);
	else if (number == 2) setupMaterial(ambient, diffuseBlue, specular, shiness);
	glRotated(cylBaseRotateAngle, 0, 1, 0);
	//////////////////////////////////
	glTranslated(0, sphere.slideY, sphere.slideZ);
	glRotated(cylAxisRotateAngle, 0, 0, 1);
	glRotated(-cuboidRotateAngle, 1, 0, 0);
		
	glRotated(90, 0, 0, 1);
	glRotated(number * 120, 0, 1, 0);
	glTranslatef(0, 0, 0.5);
	if (bWireFrame)
		sphere.DrawWireframe();
	else
		sphere.Draw();
	glPopMatrix();
}

void drawFrame(float fSize) {
	glLineWidth(2.3);
	setupMaterial(diffuseBlue, diffuseBlue, diffuseBlue, 0);

	glBegin(GL_LINE_LOOP);
	glVertex3f(-fSize / 2, fSize / 2, fSize / 2);
	glVertex3f(-fSize / 2, fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, fSize / 2, fSize / 2);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(-fSize / 2, fSize / 2, fSize / 2);
	glVertex3f(-fSize / 2, -fSize / 2, fSize / 2);
	glVertex3f(-fSize / 2, fSize / 2, -fSize / 2);
	glVertex3f(-fSize / 2, -fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, -fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, fSize / 2, fSize / 2);
	glVertex3f(fSize / 2, -fSize / 2, fSize / 2);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(fSize / 2, -fSize / 2, fSize / 2);
	glVertex3f(-fSize / 2, -fSize / 2, fSize / 2);
	glVertex3f(-fSize / 2, -fSize / 2, -fSize / 2);
	glVertex3f(fSize / 2, -fSize / 2, -fSize / 2);
	glEnd();
}

void drawDodecahedron() 
{
	glPushMatrix();
	//////////////////////////////////
	glTranslated(dode.slideX, dode.slideY, dode.slideZ);
	glScaled(0.45, 0.45, 0.45);
	if (bWireFrame)
		dode.DrawWireframe();
	else {
		if (pickDode == true)
		{
			if (cred   == true) tmpDode = 1;
			if (cgreen == true) tmpDode = 2;
			if (cblue  == true) tmpDode = 3;
			if (cpink  == true) tmpDode = 4;
			drawFrame(1.8);
		}
		if (tmpDode == 0)     setupMaterial(ambient, diffuseBlue,    specular, shiness);
		else if(tmpDode == 1) setupMaterial(ambient, diffuseRed,     specular, shiness);
		else if(tmpDode == 2) setupMaterial(ambient, diffuseGreen,   specular, shiness);
	    else if(tmpDode == 3) setupMaterial(ambient, diffuseBlue,    specular, shiness);
		else if(tmpDode == 4) setupMaterial(ambient, diffuseMagenta, specular, shiness);
		dode.Draw();
	}
	glPopMatrix();
}
void drawIsocahedron() 
{
	glPushMatrix();
	//////////////////////////////////
	glTranslated(iso.slideX, iso.slideY, iso.slideZ);
	glScaled(0.4, 0.4, 0.4);
	if (bWireFrame)
		iso.DrawWireframe();
	else {
		if (pickIso == true)
		{
			if (cred   == true) tmpIso = 1;
			if (cgreen == true) tmpIso = 2;
			if (cblue  == true) tmpIso = 3;
			if (cpink  == true) tmpIso = 4;
			drawFrame(1.9);
		}
		if (tmpIso == 0)     setupMaterial(ambient, diffuseGreen,   specular, shiness);
		else if(tmpIso == 1) setupMaterial(ambient, diffuseRed,     specular, shiness);
		else if(tmpIso == 2) setupMaterial(ambient, diffuseGreen,   specular, shiness);
	    else if(tmpIso == 3) setupMaterial(ambient, diffuseBlue,    specular, shiness);
		else if(tmpIso == 4) setupMaterial(ambient, diffuseMagenta, specular, shiness);
		iso.Draw();
	}
	glPopMatrix();
}
void drawSphere()
{
	glPushMatrix();
	//////////////////////////////////
	glTranslated(sphe.slideX, sphe.slideY, sphe.slideZ);
	if (bWireFrame)
		sphe.DrawWireframe();
	else {
		if (pickSphe == true)
		{
			if (cred   == true) tmpSphe = 1;
			if (cgreen == true) tmpSphe = 2;
			if (cblue  == true) tmpSphe = 3;
			if (cpink  == true) tmpSphe = 4;
			drawFrame(0.7);
		}
		if (tmpSphe == 0)     setupMaterial(ambient, diffuseMagenta, specular, shiness);
		else if(tmpSphe == 1) setupMaterial(ambient, diffuseRed,     specular, shiness);
		else if(tmpSphe == 2) setupMaterial(ambient, diffuseGreen,   specular, shiness);
	    else if(tmpSphe == 3) setupMaterial(ambient, diffuseBlue,    specular, shiness);
		else if(tmpSphe == 4) setupMaterial(ambient, diffuseMagenta, specular, shiness);
		sphe.Draw();
	}
	glPopMatrix();
}
void drawTruncatedcube() 
{
	glPushMatrix();
	//////////////////////////////////
	glTranslated(trun.slideX, trun.slideY, trun.slideZ);
	if (bWireFrame)
		trun.DrawWireframe();
	else {
		if (pickTrun == true)
		{
			if (cred   == true) tmpTrun = 1;
			if (cgreen == true) tmpTrun = 2;
			if (cblue  == true) tmpTrun = 3;
			if (cpink  == true) tmpTrun = 4;
			drawFrame(0.8);
		}
		if (tmpTrun == 0)     setupMaterial(ambient, diffuseGreen,   specular, shiness);
		else if(tmpTrun == 1) setupMaterial(ambient, diffuseRed,     specular, shiness);
		else if(tmpTrun == 2) setupMaterial(ambient, diffuseGreen,   specular, shiness);
	    else if(tmpTrun == 3) setupMaterial(ambient, diffuseBlue,    specular, shiness);
		else if(tmpTrun == 4) setupMaterial(ambient, diffuseMagenta, specular, shiness);
		trun.Draw();
	}
	glPopMatrix();
}

void drawMecha() {
	drawCylinderBase();
	drawCuboidBase();
	drawCylAxis();
	drawUShape();
	drawCubAxis();
	drawCylMain();
	sphere.SetColor(0);
	drawSphere(0);
	sphere.SetColor(1);
	drawSphere(1);
	sphere.SetColor(7);
	drawSphere(2);

	drawDodecahedron();
	drawIsocahedron();
	drawSphere();
	drawTruncatedcube();
}

void shadMatrix(GLfloat shadowMat[4][4], GLfloat groundplane[4], GLfloat lightpos[4])
{
  GLfloat dot;

  /* Find dot product between light position vector and ground plane normal. */
  dot = groundplane[0] * lightpos[0] +
    groundplane[1] * lightpos[1] +
    groundplane[2] * lightpos[2] +
    groundplane[3] * lightpos[3];

  shadowMat[0][0] = dot - lightpos[0] * groundplane[0];
  shadowMat[1][0] = 0.f - lightpos[0] * groundplane[1];
  shadowMat[2][0] = 0.f - lightpos[0] * groundplane[2];
  shadowMat[3][0] = 0.f - lightpos[0] * groundplane[3];

  shadowMat[0][1] = 0.f - lightpos[1] * groundplane[0];
  shadowMat[1][1] = dot - lightpos[1] * groundplane[1];
  shadowMat[2][1] = 0.f - lightpos[1] * groundplane[2];
  shadowMat[3][1] = 0.f - lightpos[1] * groundplane[3];

  shadowMat[0][2] = 0.f - lightpos[2] * groundplane[0];
  shadowMat[1][2] = 0.f - lightpos[2] * groundplane[1];
  shadowMat[2][2] = dot - lightpos[2] * groundplane[2];
  shadowMat[3][2] = 0.f - lightpos[2] * groundplane[3];

  shadowMat[0][3] = 0.f - lightpos[3] * groundplane[0];
  shadowMat[1][3] = 0.f - lightpos[3] * groundplane[1];
  shadowMat[2][3] = 0.f - lightpos[3] * groundplane[2];
  shadowMat[3][3] = dot - lightpos[3] * groundplane[3];

}

void fiPlane(GLfloat plane[4], GLfloat v0[3], GLfloat v1[3], GLfloat v2[3])
{
  GLfloat vec0[3], vec1[3];

  /* Need 2 vectors to find cross product. */
  vec0[0] = v1[0] - v0[0];
  vec0[1] = v1[1] - v0[1];
  vec0[2] = v1[2] - v0[2];

  vec1[0] = v2[0] - v0[0];
  vec1[1] = v2[1] - v0[1];
  vec1[2] = v2[2] - v0[2];

  /* find cross product to get A, B, and C of plane equation */
  plane[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
  plane[1] = -(vec0[0] * vec1[2] - vec0[2] * vec1[0]);
  plane[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];

  plane[3] = -(plane[0] * v0[0] + plane[1] * v0[1] + plane[2] * v0[2]);
}

void DisplayOneView(int nType, int left, int right, int top, int bottom)
{
	mySetupCameraVolume(nType);
	glViewport(left, top, right - left, bottom - top);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	//code here
	changeCameraPos();
	if (nType == 1) 
		gluLookAt(0, camera_dis, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	else if (nType == 2) 
		gluLookAt(0, 0.0, camera_dis, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	else if (nType == 3) 
		gluLookAt(camera_dis, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	else 
		gluLookAt(camera_X, camera_Y, camera_Z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	
	/////////////////////////////////////////////////
	fiPlane(floorPlane, floorVertices[1], floorVertices[2], floorVertices[3]);
	shadMatrix(floorShadow1, floorPlane, light_position1);
	shadMatrix(floorShadow2, floorPlane, light_position2);
	
	glPushMatrix();
    glEnable(GL_NORMALIZE);
	glCullFace(GL_FRONT);
	glScalef(1.0, -1.0, 1.0);
	drawMecha();
	glDisable(GL_NORMALIZE);
	glCullFace(GL_BACK);
	glPopMatrix();
	
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);  /* Force the 50% black. */
	glDisable(GL_DEPTH_TEST);
    glColor4f(0.0, 0.0, 0.0, 0.5);

	glPushMatrix();
	glMultMatrixf((GLfloat *) floorShadow1);
	drawMecha();
    glPopMatrix();

	if (light1) {
		glPushMatrix();
		glMultMatrixf((GLfloat *) floorShadow2);
		drawMecha();
		glPopMatrix();
	}

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	
	drawMecha();
	drawFloor();
}

void myDisplay()
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(light1 == false)
	{
		glDisable(GL_LIGHT1);
	}
	else glEnable(GL_LIGHT1);
	
	if (b4View == false)
	{
		DisplayOneView(4, 0, screenWidth, 0, screenHeight);
	}
	else
	{
		DisplayOneView(1, 0, screenWidth / 2, 0, screenHeight / 2);
		DisplayOneView(2, 0, screenWidth / 2, screenHeight / 2, screenHeight);
		DisplayOneView(3, screenWidth / 2, screenWidth, screenHeight / 2, screenHeight);
		DisplayOneView(4, screenWidth / 2, screenWidth, 0, screenHeight / 2);
	}

	glFlush();
    glutSwapBuffers();

}

void myInit()
{
	glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

	glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	float	fHalfSize = 3;

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glColor3f(0.0f, 0.0f, 0.0f);

	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
//	glEnable(GL_CULL_FACE);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-fHalfSize, fHalfSize, -fHalfSize, fHalfSize, -1000, 1000);

	camera_angle = 0.0; 
	camera_height = 2.0; 
	camera_dis = 11; 
	pickDode = pickIso = pickSphe = pickTrun = false;

	glMatrixMode(GL_MODELVIEW);
}
int main(int argc, char* argv[])
{
	cout << "1, 2: Rotate the base" << endl;
	cout << "3, 4: Rotate axes 1 (cylinder axes)" << endl;
	cout << "5, 6: Rotate axes 2 (cuboid axes)" << endl;
	cout << "W, w: Switch between wireframe and solid mode" << endl;
	cout << "V, v: to switch between 1 and 4 views." << endl;
	cout << "+   : to increase camera distance." << endl;
	cout << "-   : to decrease camera distance." << endl;
	cout << "up arrow  : to increase camera height." << endl;
	cout << "down arrow: to decrease camera height." << endl;
	cout << "<-        : to rotate camera clockwise." << endl;
	cout << "->        : to rotate camera counterclockwise." << endl;

	glutInit(&argc, (char**)argv); //initialize the tool kit
	glutInitDisplayMode(GLUT_DOUBLE |GLUT_RGB);//set the display mode
	glutInitWindowSize(screenWidth, screenHeight); //set window size
	glutInitWindowPosition(100, 100); // set window position on screen
	glutCreateWindow("Assignment1 - Vu Ngoc Linh - 51201929"); // open the screen window

	// San nha
	planeY.CreateYPlane(10, 10, 1);
	planeY.slideY = YPlanePos;
	planeY.CalculateFacesNorm();

	// Cai de xanh la
	cylBase.CreateCylinder(20, cylBaseHeight, cylBaseRadius);
	cylBase.slideY = YPlanePos + cylBaseHeight/2.0;
	cylBase.CalculateFacesNorm();
	cylBase.SetColor(1);

	// Cai than do
	cuboidBase.CreateCuboid(cuboidBaseSizeXZ/2.0, cuboidBaseSizeY/2.0, cuboidBaseSizeXZ/2.0);
	cuboidBase.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY/2.0;
	cuboidBase.CalculateFacesNorm();
	cuboidBase.SetColor(0);

	// Truc quay 1 xanh duong
	cylAxis.CreateCylinder(20, cylAxisHeight, cylAxisRadius);
	cylAxis.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY - cylAxisOffset;
	cylAxis.slideZ = cylAxisHeight/2+cuboidBaseSizeXZ/2;
	cylAxis.CalculateFacesNorm();
	cylAxis.SetColor(2);

	// Hinh chu U
	uShape.CreateUShape(1.6, 0.3, 1.4, 0.1);
	uShape.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY - cylAxisOffset;
	uShape.slideZ = 1.4/2 + cylAxisHeight + cuboidBaseSizeXZ/2;
	uShape.CalculateFacesNorm();
	uShape.SetColor(3);

	// Truc hong
	cubAxis.CreateCuboid(0.7, 0.04, 0.04);
	cubAxis.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY - cylAxisOffset;
	cubAxis.slideZ = 1.0 + cylAxisHeight + cuboidBaseSizeXZ/2;
	cubAxis.CalculateFacesNorm();
	cubAxis.SetColor(4);

	// Dia tron o giua
	cylMain.CreateCylinder(20, 0.1, 0.7);
	cylMain.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY - cylAxisOffset;
	cylMain.slideZ = 1.0 + cylAxisHeight + cuboidBaseSizeXZ/2;
	cylMain.CalculateFacesNorm();
	cylMain.SetColor(5);

	// May vien bi
	sphere.CreateSphere(10, 10, 0.15);
	sphere.slideY = YPlanePos + cylBaseHeight + cuboidBaseSizeY - cylAxisOffset;
	sphere.slideZ = 1.0 + cylAxisHeight + cuboidBaseSizeXZ/2;
	sphere.CalculateFacesNorm();

	// Cuc 12 mat
	dode.CreateDodecahedron();
	dode.slideX = 3;
	dode.slideY = YPlanePos + 0.83*0.45;
	dode.slideZ = 3;
	dode.CalculateFacesNorm();

	// Cuc 20 mat
	iso.CreateIsocahedron();
	iso.slideX = -3;
	iso.slideY = YPlanePos + 1*0.4;
	iso.slideZ = 3;
	iso.CalculateFacesNorm();

	// Cuc hinh cau
	sphe.CreateSphere(10, 10, 0.4);
	sphe.slideX = -3;
	sphe.slideY = YPlanePos + 0.4;
	sphe.slideZ = -3;
	sphe.CalculateFacesNorm();

	// Cuc vuong mai goc
	trun.CreateTruncatedcube(0.4);
	trun.slideX = 3;
	trun.slideY = YPlanePos + 0.4;
	trun.slideZ = -3;
	trun.CalculateFacesNorm();

	myInit();

	glutKeyboardFunc(myKeyboard);
    glutDisplayFunc(myDisplay);
	glutSpecialFunc(mySpecialKeyboard);
	glutTimerFunc(30, cartoonAction, 2);
	glutMouseFunc(leftClick);
	  
	glutMainLoop();
	return 0;
}