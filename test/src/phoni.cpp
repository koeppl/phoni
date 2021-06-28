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

#include <kseq.h>
#include <zlib.h>



template<class matchingstats> 
void run(const Args& args) {
  verbose("Deserializing the PHONI index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  matchingstats ms;
  {
    ifstream in(args.filename + ".phoni");
    ms.load(in, args.filename);
  }

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PHONI index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  FILE *out_fd;

  const std::string tmp_out_filename = std::string(args.patterns) + ".ms.tmp.out";

  if ((out_fd = fopen(tmp_out_filename.c_str(), "w")) == nullptr)
    error("open() file " + tmp_out_filename + " failed");

  gzFile fp = gzopen(args.patterns.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0)
  {
    auto res = ms.query(seq->seq.s, seq->seq.l);

    size_t q_length = res.first.size();
    fwrite(&q_length, sizeof(size_t), 1, out_fd);
    fwrite(res.first.data(), sizeof(size_t), q_length, out_fd);
    fwrite(res.second.data(), sizeof(size_t), q_length, out_fd);
  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(out_fd);

  t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Printing plain output");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::ofstream f_pointers(args.patterns + ".pointers");
  std::ofstream f_lengths(args.patterns + ".lengths");

  if (!f_pointers.is_open())
    error("open() file " + std::string(args.patterns) + ".pointers failed");

  if (!f_lengths.is_open())
    error("open() file " + std::string(args.patterns) + ".lengths failed");

  FILE *in_fd;

  if ((in_fd = fopen(tmp_out_filename.c_str(), "r")) == nullptr)
    error("open() file " + tmp_out_filename + " failed");

  size_t n_seq = 0;
  size_t length = 0;
  size_t m = 100; // Reserved size for pointers and lengths
  size_t *mem = (size_t *)malloc(m * sizeof(size_t));
  while (!feof(in_fd) and fread(&length, sizeof(size_t), 1, in_fd) > 0)
  {
    if (m < length)
    {
      // Resize lengths and pointers
      m = length;
      mem = (size_t *)realloc(mem, m * sizeof(size_t));
    }

    if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
      error("fread() file " + std::string(tmp_out_filename) + " failed");

    f_lengths << ">" + std::to_string(n_seq) << endl;
    for (size_t i = 0; i < length; ++i)
      f_lengths << mem[length - 1 - i] << " ";
    f_lengths << endl;

    if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
      error("fread() file " + std::string(tmp_out_filename) + " failed");

    f_pointers << ">" + std::to_string(n_seq) << endl;
    for (size_t i = 0; i < length; ++i)
      f_pointers << mem[length - 1 - i] << " ";
    f_pointers << endl;

    n_seq++;
  }
  fclose(in_fd);

  if (std::remove(tmp_out_filename.c_str()) != 0)
    error("remove() file " + tmp_out_filename + " failed");

  f_pointers.close();
  f_lengths.close();

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
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

  if(args.grammar == "naive") {
  verbose("using naive grammar");
  run<ms_pointers<PlainSlp<var_t, Fblc, Fblc>> >(args);
  } else {
  verbose("using default grammar");
  run<ms_pointers<> >(args);
    }


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
