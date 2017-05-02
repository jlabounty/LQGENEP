#include <fstream>
#include <string>

int Lepto
{

	int I, KS, KF, ORIG;
	double px, py, pz, E, m;

	std::ifstream file("TestOut.txt");
	std::string str; 
	while (std::getline(file, str))
	{
		if(str >> I >> KS >> KF >> ORIG >> px >> py >> pz >> E >> m;)
		{
			if((KS != 21) && (KF == 15))
			{
				cout << KF << "    " << E << endl;
			}
		}
	}
	


	return 0;
}
