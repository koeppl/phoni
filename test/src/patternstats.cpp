#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main(const int argc, const char *const argv[]) {
    if(argc < 2) {
        std::cerr << "usage: " << argv[0] << " file " << endl;
        return 1;
    }

    double avg = 0;
    size_t max = 0;
    size_t reducible = 0;

    size_t len = 0;

    {
        std::string lenfile = argv[1];
        lenfile += ".lengths";
        ifstream is(lenfile);
        while(is) {
            std::string line;
            getline(is, line);
            if(line[0] == '>') continue;
            std::istringstream iss(line);
			size_t number;
			for (std::istringstream numbers_iss(line); numbers_iss >> number; ) {
				if(number > max) { max = number; }
				avg += number;
				++len;
	//		std::cout << number << endl;
			}
        }
        is.close();
    }
    {
        std::string reffile = argv[1];
        reffile += ".pointers";
        ifstream is(reffile);
        size_t prev = -1;
        while(is) {
            std::string line;
            getline(is, line);
            if(line[0] == '>') continue;
            std::istringstream iss(line);
			size_t number;
			for (std::istringstream numbers_iss(line); numbers_iss >> number; ) {
				if(number == prev+1) {
					++reducible;
				}
				prev = number;
				//std::cout << number << endl;
			}
                        prev = -1;


        }
    }
    
     std::cout << " size= " << len << " avg= " << (avg/len) << " max=" << max << " reducible=" << reducible << std::endl;

    return 0;
}
