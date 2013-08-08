#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <cmath>

using namespace std;
//pbc-round function to deal with periodic boundary conditions
//#############################################
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
//############################################

int main()
{
	// Main menu -- get info from user
	//######################################
	string infile;
	int nooa, noha;	
	double xlattice, ylattice, zlattice, timestep;

	cout << "XYZ file\n==> ";
	cin >> infile;
	cout << "Number of O atoms\n==> ";
	cin >> nooa;
	cout << "Number of H atoms\n==> ";
	cin >> noha;
	cout << "Lattice constants (x y z)\n==> ";
	cin >> xlattice >> ylattice >> zlattice;
	cout << "Timestep (fs)\n==> ";
	cin >> timestep;
	//#######################################
	
	// Read inputfile and place data into respective vectors
	//####################################################
	ifstream input;
	input.open(infile.c_str());

	string atom;
	double ox, oy, oz, hx, hy, hz;
	vector <double> oxs, oys, ozs, hxs, hys, hzs;

	while (! input.eof())
	{
		input >> atom;
		if (atom == "O")
		{
			input >> ox >> oy >> oz;
			oxs.push_back(ox);
			oys.push_back(oy);
			ozs.push_back(oz);			
		}
		if (atom == "H")
		{
			input >> hx >> hy >> hz;		
			hxs.push_back(hx);
			hys.push_back(hy);
			hzs.push_back(hz);
		}
	}	

	int nframes = oxs.size() / nooa;
	//######################################################

	// Create and initialize the O-H partners array
	//#############################################
	int partners[nooa*nframes][4];
	for (int i = 0; i < nooa*nframes; i ++)
	{
		for (int j = 0; j < 4; j ++)
		{
			partners[i][j] = -1;
		}
	}
	//#############################################

	// Go through the inputfiles data and perform designated calculations.
	// This is the main portion of the code.
	//####################################################################
	for (int i = 0; i < nframes; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			int hcount = 0;

			for (int k = 0; k < noha; k ++)
			{
				double dx = oxs[j + i*nooa] - hxs[k + i*noha];
				double dy = oys[j + i*nooa] - hys[k + i*noha];
				double dz = ozs[j + i*nooa] - hzs[k + i*noha];
				
				dx -= xlattice*pbc_round(dx/xlattice);
        	                dy -= ylattice*pbc_round(dy/ylattice);
	                        dz -= zlattice*pbc_round(dz/zlattice);
			
				double distance = sqrt ( dx*dx + dy*dy + dz*dz );
				
				if (distance < 1.2)
				{
					partners[j + i*nooa][hcount] = k;
					hcount ++;
				}
			}
		}
	}

	vector <int> oindex, frame;
	ofstream output;
        output.open("atomcenter.xyz");

	for (int i = 0; i < nframes; i ++)
	{

		for (int k = 0; k < nooa; k ++)
		{
			int count = 0;
		
			for (int j = 0; j < 4; j ++)
			{
				if (partners[k + i*nooa][j] != -1)
				{
				count ++;
				}
			}
		
			if (count == 3)
			{
				frame.push_back(i);
				oindex.push_back(k);
			}
		}
	}

 	vector <int> newframe, newoindex;
	for (int i = 1; i < frame.size(); i ++)
	{
		int o1 = oindex[i-1];
		int o2 = oindex[i];
		int oframe = frame[i];

		double dx = oxs[o2 + oframe*nooa] - oxs[o1 + oframe*nooa];
                double dy = oys[o2 + oframe*nooa] - oys[o1 + oframe*nooa];
                double dz = ozs[o2 + oframe*nooa] - ozs[o1 + oframe*nooa];

                dx -= xlattice*pbc_round(dx/xlattice);
                dy -= ylattice*pbc_round(dy/ylattice);
                dz -= zlattice*pbc_round(dz/zlattice);

                double distance = sqrt ( dx*dx + dy*dy + dz*dz );

                if (distance < 3.6)
		{
			newframe.push_back(oframe);
			newoindex.push_back(o1);			
		}
	}

	
	
	for (int i = 0; i < nframes; i ++)
	{
		double minusx = 0, minusy = 0, minusz = 0;
		output << nooa + noha << endl;
		output << "frame: " << i << endl;

		for (int j = 0; j < newframe.size(); j++)
		{
			if (i == newframe[j])
			{
				minusx = oxs[newoindex[j] + i*nooa];
				minusy = oys[newoindex[j] + i*nooa];
				minusz = ozs[newoindex[j] + i*nooa];
			}
		}
		
		for (int k = 0; k < nooa; k ++)
		{
			oxs[k + i*nooa] -= minusx + 0.5*xlattice;
			oys[k + i*nooa] -= minusy + ylattice*0.5;
			ozs[k + i*nooa] -= minusz + zlattice*0.5;
			
			if (oxs[k + i*nooa] > xlattice)
			{
				oxs[k + i*nooa] -= xlattice;
			}
			if (oxs[k + i*nooa] < 0.0)
			{
				oxs[k + i*nooa] += xlattice;
			}
			if (oys[k + i*nooa] > ylattice)
                        {
                                oys[k + i*nooa] -= ylattice;
                        }
                        if (oys[k + i*nooa] < 0.0)
                        {
                                oys[k + i*nooa] += ylattice;
                        }
			if (ozs[k + i*nooa] > zlattice)
                        {
                                ozs[k + i*nooa] -= zlattice;
                        }
                        if (ozs[k + i*nooa] < 0.0)
                        {
                                ozs[k + i*nooa] += zlattice;
                        }
			output << "O  " << oxs[k + i*nooa] << "  " << oys[k + i*nooa] << "  " << ozs[k + i*nooa] << endl;
		}
		for (int n = 0; n < noha; n ++)
		{
			hxs[n + i*noha] -= minusx + 0.5*xlattice;
			hys[n + i*noha] -= minusy + 0.5*ylattice;
			hzs[n + i*noha] -= minusz + 0.5*zlattice;
			
			if (hxs[n + i*noha] > xlattice)
                        {
                                hxs[n + i*noha] -= xlattice;
                        }
                        if (hxs[n + i*noha] < 0.0)
                        {
                                hxs[n + i*noha] += xlattice;
                        }
                        if (hys[n + i*noha] > ylattice)
                        {
                                hys[n + i*noha] -= ylattice;
                        }
                        if (hys[n + i*noha] < 0.0)
                        {
                                hys[n + i*noha] += ylattice;
                        }
                        if (hzs[n + i*noha] > zlattice)
                        {
                                hzs[n + i*noha] -= zlattice;
                        }
                        if (hzs[n + i*noha] < 0.0)
                        {
                                hzs[n + i*noha] += zlattice;
                        }
			output << "H  " << hxs[n + i*noha] << "  " << hys[n + i*noha] << "  " << hzs[n + i*noha] << endl;
		}	
			
				
	}	

	// close files and terminate program
	//################################
	input.close();
	output.close();
	return 0;
	//################################
}

	//#########################################################################	





