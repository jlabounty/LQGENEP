int BuildTree()
{
	std::string inputFile = "./TestOut.txt";
	std::string outdir = "./";
	int events = 1000000;

	gSystem->Load("libeicsmear");
	BuildTree(inputFile.c_str(), outdir.c_str(), events);
	cout << "********************************************" << endl;
	cout << "File: " << inputFile << " processed in directory " << outdir << endl;

	return 0;
}
