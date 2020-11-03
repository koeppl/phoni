#include <fstream>
#include <iostream>

using namespace std;

    inline static void write_int(ostream& os, const size_t i){
      os.write(reinterpret_cast<const char*>(&i), 5);
    }

int main(int argc, char *const argv[]) {
	if(argc < 2) {
		std::cerr << "Usage: " << argv[0] << " textfilename " << std::endl;
	}
	const std::string s = argv[1];
	ifstream bwtfile(s + ".bwt");

	ofstream headfile(s + ".bwt.heads", std::ios::binary);
	ofstream lenfile(s + ".bwt.len", std::ios::binary);

	char oldchar = bwtfile.get();
	size_t runlength = 1;
	
	auto print_run = [&] () {
		headfile.put(oldchar);
		write_int(lenfile, runlength);
	};

	while(bwtfile) {
		const char c = bwtfile.get();
		if(!bwtfile) { break; }
		if(c == oldchar) {
			++runlength;
		} else {
			print_run();
			oldchar = c; runlength = 1;
		}
	}
	print_run();

	return 0;
}
