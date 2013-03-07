/*
Spencer Egart
1/30/2013

Reads in a PDB file and displays the protein in a rotating point cloud.
*/


#include <GL\freeglut.h>
#include "Geometry.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace csc480;
using namespace std;

//Rotation speed
const double SECONDS_PER_ROT= 16 ;

//Display size
const double BOX_SIZE = 1;

//The indices of the important parts of the PDB entries for an atom
const int X_INDEX = 30, Y_INDEX = 38, Z_INDEX = 46, ELEMENT_INDEX = 77, COORD_LENGTH = 8;

//Protein metadata
string header;
vector<string> titles;

//Display lists
int G_dlistID_axes, G_dlistID_protein;

//Holds the relevant data for an atom
struct Atom
{
	Point3 point; char element;

	explicit Atom(Point3 point, char element) : point(point), element(element) {}

	//Scale the coordinates
	void scale(const double scale)
	{
		point.x *= scale;
		point.y *= scale;
		point.z *= scale;
	}

	//Translate the coordinates
	void translate(const double x, const double y, const double z)
	{
		point.x += x;
		point.y += y;
		point.z += z;
	}

	double distance(const double x, const double y, const double z)
	{
		double distX = point.x - x;
		double distY = point.y - y;
		double distZ = point.z - z;
		return sqrt(distX * distX + distY * distY + distZ * distZ);
	}
};

//Using a vector of Atom structs to represent a Protein
typedef vector<Atom> Protein;

//Operator overload for printing an atom
inline std::ostream& operator<<( std::ostream& out, const Atom& atom ) // out << atom
{ out << atom.element << "( " << atom.point.x << " " << atom.point.y << " " << atom.point.z << " )" ; return out ; }

//Translate the Protein so that its center of mass is on the origin
void centerAtoms(Protein &protein)
{
	double centerX = 0, centerY = 0, centerZ = 0;
	for(Protein::size_type i = 0; i < protein.size(); i++)
	{
		centerX += protein[i].point.x;
		centerY += protein[i].point.y;
		centerZ += protein[i].point.z;
	}
	centerX /= protein.size();
	centerY /= protein.size();
	centerZ /= protein.size();
	for(Protein::size_type i = 0; i < protein.size(); i++)
	{
		protein[i].translate(-centerX, -centerY, -centerZ);
	}
}

//Determine the max distance from the origin in the Protein
double maxDistance(Protein &protein)
{
	
	double rMax = 0;
	for(Protein::size_type i = 0; i < protein.size(); i++)
	{
		Atom atom = protein[i];

		//Distance from (0,0,0)
		double dist = atom.distance(0, 0, 0);

		if(dist > rMax)
		{
			rMax = dist;
		}
	}

	return rMax;
}

//Scale a vector of protein to our viewport
void scaleAtoms(Protein &protein)
{
	double rMax = maxDistance(protein);
	//Scale all the protein accordingly
	double scale = 1/(rMax/(BOX_SIZE/2));
	for(Protein::size_type i = 0; i < protein.size(); i++)
	{
		protein[i].scale(scale);
	}
}

//Read a pdb file, then scale the results to our viewport.
Protein readPDB(string filename)
{
	Protein protein;

	//Open the file
	ifstream pdbFile;
	pdbFile.open(filename, ios::in);

	string line;

	//Go through the file and read in the protein
	if(pdbFile.is_open())
	{
		while(pdbFile.good())
		{
			//Read a line
			getline(pdbFile, line);

			if(line.size() > 0)
			{
				//We have an atom
				if(line.substr(0,4).compare("ATOM") == 0)
				{
					//Extract the relevant info based on the PDB format
					double x = atof(line.substr(X_INDEX, COORD_LENGTH).c_str());
					double y = atof(line.substr(Y_INDEX, COORD_LENGTH).c_str());
					double z = atof(line.substr(Z_INDEX, COORD_LENGTH).c_str());
					char element = line[ELEMENT_INDEX];

					protein.push_back( Atom( Point3(x, y, z), element) );
				}
				//We're reading in the header
				else if(line.substr(0, 6).compare("HEADER") == 0)
				{
					header += line.substr(10, 40);
				}
				//Reading in a title
				else if(line.substr(0, 5).compare("TITLE") == 0)
				{
					titles.push_back(line.substr(10));
				}
			}
		}
		//Close the file
		pdbFile.close();
	}

	return protein;
}

void normalize(Protein &protein)
{
	//Center and scale the protein to our viewport
	centerAtoms(protein);
	scaleAtoms(protein);

	//Make sure everything is in the box
	for(Atom atom: protein)
	{
		assert(atom.point.x <= BOX_SIZE/2);
		assert(atom.point.x >= -BOX_SIZE/2);
		assert(atom.point.y <= BOX_SIZE/2);
		assert(atom.point.y >= -BOX_SIZE/2);
		assert(atom.point.z <= BOX_SIZE/2);
		assert(atom.point.z >= -BOX_SIZE/2);
	}
}

//Write a text label
void drawText(string text, double x, double y)
{
	glRasterPos2d(x, y);
	glColor3d(0, 0, 0);
	for(string::size_type i = 0; i < text.size(); i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, text[i]);
}

//Draw the axis planes (not messerschmitts)
void drawAxes()
{
	glColor3d(0, 0, 0);
	//x-y plane
	glBegin(GL_LINE_LOOP);
	glVertex3d(-1* (BOX_SIZE/2), (BOX_SIZE/2), 0);
	glVertex3d((BOX_SIZE/2), (BOX_SIZE/2), 0);
	glVertex3d((BOX_SIZE/2), -1* (BOX_SIZE/2), 0);
	glVertex3d(-1* (BOX_SIZE/2), -1* (BOX_SIZE/2), 0);
	glEnd();

	//x-z plane
	glBegin(GL_LINE_LOOP);
	glVertex3d(-1* (BOX_SIZE/2), 0, (BOX_SIZE/2));
	glVertex3d((BOX_SIZE/2), 0, (BOX_SIZE/2));
	glVertex3d((BOX_SIZE/2), 0, -1* (BOX_SIZE/2));
	glVertex3d(-1* (BOX_SIZE/2), 0, -1* (BOX_SIZE/2));
	glEnd();

	//y-z plane
	glBegin(GL_LINE_LOOP);
	glVertex3d(0, -1* (BOX_SIZE/2), (BOX_SIZE/2));
	glVertex3d(0, (BOX_SIZE/2), (BOX_SIZE/2));
	glVertex3d(0, (BOX_SIZE/2), -1* (BOX_SIZE/2));
	glVertex3d(0, -1* (BOX_SIZE/2), -1* (BOX_SIZE/2));
	glEnd();
}

//Draw a point for each atom
void drawAtoms(const Protein &protein)
{

	glBegin(GL_POINTS);
	for(Atom atom: protein)
	{
		

		//Get the color for the atom
		switch(atom.element)
		{
		case 'O': glColor3d(1, 0, 0); break;
		case 'C': glColor3d(0, 1, 0); break;
		case 'N': glColor3d(0, 0, 1); break;
		case 'S': glColor3d(1, 1, 0); break;
		default: glColor3d(.5, .5, .5); break; //Some of the files had other elements
		}
		//Draw the atom
		glVertex3d(atom.point.x, atom.point.y, atom.point.z);
	}
	glEnd();
}

//Make a display list for the axes
int makeDisplayList_Axes()
{
	int dlistID = glGenLists(1);
	glNewList(dlistID, GL_COMPILE);
	drawAxes();
	glEndList();
	return dlistID;
}

//Make a display list for the protein
int makeDisplayList_Protein(const Protein &protein)
{
	int dlistID = glGenLists(1);
	glNewList(dlistID, GL_COMPILE);
	drawAtoms(protein);
	glEndList();
	return dlistID;
}

//Draws the PDB classification header and title line(s)
void writeMetadata()
{
	//Determine the line height based on the bitmap font.
	double height = (glutBitmapHeight(GLUT_BITMAP_HELVETICA_10))*1.5/glutGet(GLUT_WINDOW_HEIGHT);
	
	//Find a starting Y that will let every line fit.
	double textYPos = (-1 + height)  + (height * (titles.size() +1) );

	//Header
	drawText(header, -.9, textYPos);
	
	//Title(s)
	for(vector<string>::size_type i = 0; i < titles.size() ; ++i)
	{
		drawText(titles[i], -.9, textYPos + (i+1)*(-height));
	}
}

//Main display method, where the magic happens
void display()
{
	//Reset and clear
	glLoadIdentity() ;
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Write out the PDB metadata, for fun
	writeMetadata();


	//Tumble
	double degrees = glutGet( GLUT_ELAPSED_TIME ) * 0.001* 360 / SECONDS_PER_ROT ;
	glRotated( degrees, 1,1, 1 ) ;
	glRotated( degrees, 1,0, 0 ) ;

	//Draw everything
	glCallList( G_dlistID_axes ) ;
	glCallList( G_dlistID_protein);

	//update
	glutSwapBuffers();
}

//Just redraw every frame
void animate()
{
	glutPostRedisplay() ;
}

//Set up the window
void glInit()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

	glutInitWindowSize(600, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("CSC 580 Program 1");
}

//This is where we start, at the bottom because of C++
int main(int argc, char** argv)
{
	//Pass commandline args to glut, not that we use them
	glutInit(&argc, argv);


	//Get user input
	string filename;
	cout << "Enter a PDB filename: ";
	cin >> filename;

	//Load the file
	Protein protein = readPDB(filename);

	//Check if file reading was successful
	if(protein.size() > 0)
	{
		//Normalize the data
		normalize(protein);

		//GL init stuff
		glInit();

		//Make display lists
		G_dlistID_axes = makeDisplayList_Axes();
		G_dlistID_protein = makeDisplayList_Protein(protein);

		//Antialiasing stuff
		glEnable (GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		//Set background and point size
		glClearColor(1, 1, 1, 1);
		glPointSize(2);

		//Set display functions
		glutDisplayFunc(display);
		glutIdleFunc(animate);

		//Start it up
		glutMainLoop();

	}
	else //Something went wrong
	{
		cout << "Error reading PDB file " << filename << endl;
		system("PAUSE");
	}

	//Just because
	return 0;
}
