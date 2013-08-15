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

	// Wrap the coordinates
	//######################################################
	for (int i = 0; i < oxs.size(); i ++)
	{
		if (oxs[i] > xlattice)
		{
			oxs[i] -= xlattice;
		}
		if (oxs[i] < 0.0)
		{
			oxs[i] += xlattice;
		}
		if (oys[i] > ylattice)
                {
                        oys[i] -= ylattice;
                }
                if (oys[i] < 0.0)
                {
                        oys[i] += ylattice;
                }
		if (ozs[i] > zlattice)
                {
                        ozs[i] -= zlattice;
                }
                if (ozs[i] < 0.0)
                {
                        ozs[i] += zlattice;
                }
	}
	for (int i = 0; i < hxs.size(); i ++)
	{
		if (hxs[i] > xlattice)
                {
                        hxs[i] -= xlattice;
                }
                if (hxs[i] < 0.0)
                {
                        hxs[i] += xlattice;
                }
                if (hys[i] > ylattice)
                {
                        hys[i] -= ylattice;
                }
                if (hys[i] < 0.0)
                {
                        hys[i] += ylattice;
                }
                if (hzs[i] > zlattice)
                {
                        hzs[i] -= zlattice;
                }
                if (hzs[i] < 0.0)
                {
                        hzs[i] += zlattice;
                }
        }
	//#######################################################



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
	ofstream output;
	output.open("atomcenter.xyz");
	int oindex;
	double xinc, yinc, zinc;

	for (int i = 0; i < nframes; i ++)
	{
		int oindex = -1;
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
					hcount ++;
				}
			}
		
			if (hcount == 3)
			{
				oindex = j;
				xinc = oxs[j + i*nooa];
				yinc = oys[j + i*nooa];
				zinc = ozs[j + i*nooa];
			}
			
		}
		output << nooa + noha << "\n\n";
		if (oindex != -1)
		{
		for ( int j = 0; j < nooa; j ++)
		{
			oxs[j + i*nooa] -= xinc;
			oxs[j + i*nooa] += 0.5*xlattice;
			oys[j + i*nooa] -= yinc;
			oys[j + i*nooa] += 0.5*ylattice;
			ozs[j + i*nooa] -= zinc; 
			ozs[j + i*nooa] += 0.5*zlattice;
		}
		for ( int j = 0; j < noha; j ++)
                {
                        hxs[j + i*noha] -= xinc;   
                        hxs[j + i*noha] += 0.5*xlattice;
                        hys[j + i*noha] -= yinc;
                        hys[j + i*noha] += 0.5*ylattice;
                        hzs[j + i*noha] -= zinc; 
                        hzs[j + i*noha] += 0.5*zlattice;
                }
		for (int j = 0; j < nooa; j ++)
        	{
                if (oxs[j + i*nooa] > xlattice)
                {
                        oxs[j + i*nooa] -= xlattice;
                }
                if (oxs[j + i*nooa] < 0.0)
                {
                        oxs[j + i*nooa] += xlattice;
                }
                if (oys[j + i*nooa] > ylattice)
                {
                        oys[j + i*nooa] -= ylattice;
                }
                if (oys[j + i*nooa] < 0.0)
                {
                        oys[j + i*nooa] += ylattice;
                }
                if (ozs[j + i*nooa] > zlattice)
                {
                        ozs[j + i*nooa] -= zlattice;
                }
                if (ozs[j + i*nooa] < 0.0)
                {
                        ozs[j + i*nooa] += zlattice;
                }
        }
        for (int j = 0; j < noha; j ++)
        {
                if (hxs[j + i*noha] > xlattice)
                {
                        hxs[j + i*noha] -= xlattice;
                }
                if (hxs[j + i*noha] < 0.0)
                {
                        hxs[j + i*noha] += xlattice;
                }
                if (hys[j + i*noha] > ylattice)
                {
                        hys[j + i*noha] -= ylattice;
                }
                if (hys[j + i*noha] < 0.0)
                {
                        hys[j + i*noha] += ylattice;
                }
                if (hzs[j + i*noha] > zlattice)
                {
                        hzs[j + i*noha] -= zlattice;
                }
                if (hzs[j + i*noha] < 0.0)
                {
                        hzs[j + i*noha] += zlattice;
                }
        }

		for ( int j = 0; j < nooa; j++)
		{
			output << "O  " << oxs[j + i*nooa] << "  " << oys[j +i*nooa] << "  " << ozs[j + i*nooa] << endl;
		}
		for ( int j = 0; j < noha; j++)
                {
                        output << "H  " << hxs[j + i*noha] << "  " << hys[j +i*noha] << "  " << hzs[j + i*noha] << endl;
                }

}
else
{
	for ( int j = 0; j < nooa; j++)
                {
                        output << "O  " << oxs[j + i*nooa] << "  " << oys[j +i*nooa] << "  " << ozs[j + i*nooa] << endl;
                }
                for ( int j = 0; j < noha; j++)
                {
                        output << "H  " << hxs[j + i*noha] << "  " << hys[j +i*noha] << "  " << hzs[j + i*noha] << endl;
                }
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





