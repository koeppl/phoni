/* shapedslp_test - Test if the ShapedSLP generates the whole original text
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
   \file shapedslp_test.cpp
   \brief shapedslp_test.cpp Test if the ShapedSLP generates the whole original text.
   \author Massimiliano Rossi
   \date 18/09/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>


int main(int argc, char *const argv[])
{
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;

  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index
  verbose("Building random access");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // pfp_ra ra(args.filename, args.w);
  std::string filename_slp = args.filename + ".slp";
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;
  ifstream fs(filename_slp);
  ra.load(fs);

  size_t n = ra.getLen();

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Reading text");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::vector<uint8_t> text;
  if(args.is_fasta)
    read_fasta_file(args.filename.c_str(),text);
  else
    read_file(args.filename.c_str(),text);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Checking if equals");
  t_insert_start = std::chrono::high_resolution_clock::now();

  if(n != text.size())
    error("Text size is different", " ra: ", n, " text: ", text.size());

  for(size_t i = 0; i < n; ++i)
    if(ra.charAt(i) != text[i])
      error("Different character in position ", i, " ra: ", ra.charAt(i), " text: ", text[i]);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}