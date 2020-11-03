/* rlbwt - Computes the run-length encoding of the BWT
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
   \file rlbwt.cpp
   \brief rlbwt.cpp Computes the run-length encoding of the BWT.
   \author Massimiliano Rossi
   \date 02/11/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <malloc_count.h>

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index

  verbose("Compressing the BWT");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Open output files
  std::string bwt_filename = args.filename;
  std::string bwt_heads_filename = args.filename + ".heads";
  std::string bwt_len_filename = args.filename + ".len";

  FILE *bwt;
  FILE *bwt_heads;
  FILE *bwt_len;

  if ((bwt = fopen(bwt_filename.c_str(), "r")) == nullptr)
    error("open() file " + std::string(bwt_filename) + " failed");

  if ((bwt_heads = fopen(bwt_heads_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(bwt_heads_filename) + " failed");

  if ((bwt_len = fopen(bwt_len_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(bwt_len_filename) + " failed");

  // Get the BWT length
  struct stat filestat;
  int fn = fileno(bwt);

  if (fstat(fn, &filestat) < 0)
      error("stat() file " + std::string(bwt_filename) + " failed" );

  if(filestat.st_size % sizeof(uint8_t) != 0)
      error("invilid file " + std::string(bwt_filename));

  size_t n = filestat.st_size;

  // Start processing
  uint8_t prev, curr;
  size_t len = 1;

  if ((fread(&prev, sizeof(uint8_t), 1, bwt)) != 1)
    error("fread() file " + std::string(bwt_filename) + " failed");
  
  for(size_t i = 0; i < n-1; ++i)
  {
    if ((fread(&curr, sizeof(uint8_t), 1, bwt)) != 1)
      error("fread() file " + std::string(bwt_filename) + " failed");
    if( curr == prev )
    {
      len ++;
    }
    else
    {
      if ((fwrite(&prev, sizeof(uint8_t), 1, bwt_heads)) != 1)
      error("fwrite() file " + std::string(bwt_heads_filename) + " failed");

      if ((fwrite(&len, 5, 1, bwt_len)) != 1)
        error("fwrite() file " + std::string(bwt_len_filename) + " failed");
      
      prev = curr;
      len = 1;
    }
  }

  if ((fwrite(&prev, sizeof(uint8_t), 1, bwt_heads)) != 1)
        error("fwrite() file " + std::string(bwt_heads_filename) + " failed");

  if ((fwrite(&len, 5, 1, bwt_len)) != 1)
    error("fwrite() file " + std::string(bwt_len_filename) + " failed");

  fclose(bwt);
  fclose(bwt_heads);
  fclose(bwt_len);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

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