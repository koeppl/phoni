/* matching_statistics - Computes the matching statistics from BWT and Thresholds
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file matching_statistics.cpp
   \brief matching_statistics.cpp Computes the matching statistics from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <phoni.hpp>

#include <malloc_count.h>


std::vector<std::string> read_pattern_desc(const std::string& patternpath) {
  std::vector<std::string> descs;
  ifstream is(patternpath + "desc.txt");
  std::string line;
  while (std::getline(is, line)) {
    descs.push_back(line);
  }
  return descs;
}

   inline static void read_int(istream& is, size_t& i){
      is.read(reinterpret_cast<char*>(&i), sizeof(size_t));
    }

int main(int argc, char *const argv[]) {
  Args args;
  parseArgs(argc, argv, args);

#ifdef NDEBUG
  verbose("RELEASE build");
#else
  verbose("DEBUG build");
#endif

  verbose("Memory peak: ", malloc_count_peak());

  verbose("Deserializing the PHONI index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();


  ms_pointers<> ms;
  {
    ifstream in(args.filename + ".phoni");
    ms.load(in, args.filename);
  }

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PHONI index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Reading patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();
  const std::string patterndir = args.patterns + ".dir/";

  std::vector<std::string> patterndescs = read_pattern_desc(patterndir);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::ofstream f_pointers(args.patterns + ".pointers");
  std::ofstream f_lengths(args.patterns + ".lengths");

  if (!f_pointers.is_open())
    error("open() file " + std::string(args.filename) + ".pointers failed");

  if (!f_lengths.is_open())
    error("open() file " + std::string(args.filename) + ".lengths failed");

  for(size_t patternid = 0; patternid < patterndescs.size(); ++patternid){
    const std::string& patterndesc = patterndescs[patternid];
    verbose("Processing pattern ", patterndesc);
    f_lengths << ">" << patterndesc << " " << endl;
    f_pointers << ">" << patterndesc << " " << endl;
    const std::string patternfilename = patterndir +  std::to_string(patternid);
    //auto [lengths,refs] = 
    
    t_insert_start = std::chrono::high_resolution_clock::now();
    const size_t patternlength = ms.query(patternfilename, std::string(args.filename) + ".binrev.length",  std::string(args.filename) + ".binrev.pointers");
    t_insert_end = std::chrono::high_resolution_clock::now();
    verbose("Finished processing pattern ", patterndesc);
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    {
      const size_t* len_file;
      const size_t* ref_file;
      size_t len_size;
      size_t ref_size;

      const std::string len_filename = std::string(args.filename) + ".binrev.length";
      const std::string ref_filename = std::string(args.filename) + ".binrev.pointers";

      map_file(len_filename.c_str(), len_file, len_size);
      map_file(ref_filename.c_str(), ref_file, ref_size);

      // ifstream len_file(std::string(args.filename) + ".binrev.length", std::ios::binary);
      // ifstream ref_file(std::string(args.filename) + ".binrev.pointers", std::ios::binary);
      for(size_t i = 0; i < patternlength; ++i) {
        const size_t len = len_file[patternlength-i-1];
        const size_t ref = ref_file[patternlength-i-1];
        // size_t len;
        // size_t ref;
        // len_file.seekg((patternlength-i-1)*sizeof(size_t), std::ios_base::beg);
        // ref_file.seekg((patternlength-i-1)*sizeof(size_t), std::ios_base::beg);
        // read_int(len_file, len);
        // read_int(ref_file, ref);
        // DCHECK_EQ(lengths[i], len);
        // DCHECK_EQ(refs[i], ref);
        f_lengths << len << " ";
        f_pointers << ref << " ";
      }
    f_lengths << endl;
    f_pointers << endl;
    }
  }

  f_pointers.close();
  f_lengths.close();

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    verbose("Thresholds size (bytes): ", space);
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}
