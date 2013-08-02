#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <cmath>

using namespace std;

int pbc_round(double input)
{
	int i  = input;

	if (abs(input - i) >= 0.5)
	{
		if (input > 0) {i += 1;}
		if (input < 0) {i -= 1;}
	}
return i;
}

int main()
{

	string infile;
	int nooa, noha;	
	double xlattice, ylattice, zlattice;

	cout << "XYZ file\n==> ";
	cin >> infile;
	cout << "Number of O atoms\n==> ";
	cin >> nooa;
	cout << "Number of H atoms\n==> ";
	cin >> noha;
	cout << "Lattice constants (x y z)\n==> ";
	cin >> xlattice >> ylattice >> zlattice;
	
	ifstream input;
	input.open(infile.c_str());

	string atom;
	double ox, oy, oz, hx, hy, hz;
	vector <double> oxs, oys, ozs, hxs, hys, hzs;

	while (! infile.eof())
	{
		infile >> atom;
		if (atom == "O")
		{
			infile >> ox >> oy >> oz;
			oxs.push_back(ox);
			oys.push_back(oy);
			ozs.push_back(oz);			
		}
		if (atom == "H")
		{
			infile >> hx >> hy >> hz;		
			hxs.push_back(hx);
			hys.push_back(hy);
			hzs.push_back(hz);
		}
	}

	int nframes = oxs.size() / nooa;
	
	for (int i = 0; i < nframes; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			for (int k = 0; k < noha; k ++)
			{
				double dx = oxs[j + i*nooa] - hxs[k + i*noha];
				double dy = oys[j + i*nooa] - hys[k + i*noha];
				double dz = ozs[j + i*nooa] - hzs[k + i*noha];
				
				dx -= xlattice*pbc_round(dx/xlattice);
        	                dy -= ylattice*pbc_round(dy/ylattice);
	                        dz -= zlattice*pbc_round(dz/zlattice);
			
				double distance = sqrt ( dx*dx + dy*dy + dz*dz );
				
				if (distance < 1.3








