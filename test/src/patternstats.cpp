#include <iostream>
#include <fstream>
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
        lenfile += ".binrev.length";
        ifstream is(lenfile);
        while(is) {
            size_t i;
            is.read((char*)&i, sizeof(i));
            if(!is) { break; }

            if(i > max) { max = i; }
            avg += i;
            ++len;
        }
        is.close();
    }
    {
        std::string reffile = argv[1];
        reffile += ".binrev.pointers";
        ifstream is(reffile);
        size_t prev = 0;
        while(is) {
            size_t i;
            is.read((char*)&i, sizeof(i));
            if(!is) { break; }
            if(i == prev-1) {
                ++reducible;
            }
            prev = i;
            std::cout << i << endl;
        }
    }
    
     std::cout << " size= " << len << " avg= " << (avg/len) << " max=" << max << " reducible=" << reducible << std::endl;

    return 0;
}
