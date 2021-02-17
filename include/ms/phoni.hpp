/* ms_pointers - Computes the matching statistics pointers from BWT and Thresholds 
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
   \file ms_pointers.hpp
   \brief ms_pointers.hpp Computes the matching statistics pointers from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _MS_POINTERS_HH
#define _MS_POINTERS_HH


/** FLAGS **/
#define MEASURE_TIME 1  //measure the time for LCE and backward search?
#define NAIVE_LCE_SCHEDULE //stupidly execute two LCEs without heurstics
#ifndef NAIVE_LCE_SCHEDULE //apply a heuristic
//#define SORT_BY_DISTANCE_HEURISTIC 1 // apply a heuristic to compute the LCE with the closer BWT position first
#endif

#ifndef DCHECK_HPP
#define DCHECK_HPP
#include <string>
#include <sstream>
#include <stdexcept>

#ifndef DCHECK
#ifdef NDEBUG
#define ON_DEBUG(x)
#define DCHECK_(x,y,z)
#define DCHECK(x) 
#define DCHECK_EQ(x, y) 
#define DCHECK_NE(x, y) 
#define DCHECK_LE(x, y) 
#define DCHECK_LT(x, y) 
#define DCHECK_GE(x, y) 
#define DCHECK_GT(x, y) 
#else//NDEBUG
#define ON_DEBUG(x) x
#define DCHECK_(x,y,z) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK(x) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x))
#define DCHECK_EQ(x, y) DCHECK_((x) == (y), x,y)
#define DCHECK_NE(x, y) DCHECK_((x) != (y), x,y)
#define DCHECK_LE(x, y) DCHECK_((x) <= (y), x,y)
#define DCHECK_LT(x, y) DCHECK_((x) < (y) ,x,y)
#define DCHECK_GE(x, y) DCHECK_((x) >= (y),x,y )
#define DCHECK_GT(x, y) DCHECK_((x) > (y) ,x,y)
#endif //NDEBUG
#endif //DCHECK
#endif /* DCHECK_HPP */


#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include<ms_rle_string.hpp>

#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "ShapedSlp.hpp"
#include "ShapedSlpV2.hpp"
#include "SelfShapedSlp.hpp"
#include "SelfShapedSlpV2.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "IncBitLenCode.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"

#ifdef MEASURE_TIME
struct Stopwatch {
	std::chrono::high_resolution_clock::time_point m_val = std::chrono::high_resolution_clock::now();

	void reset() {
		m_val = std::chrono::high_resolution_clock::now();
	}
	double seconds() const { 
		return std::chrono::duration<double, std::ratio<1>>(std::chrono::high_resolution_clock::now() - m_val).count();
	}

};
#endif//MEASURE_TIME

using var_t = uint32_t;
using Fblc = FixedBitLenCode<>;
using SelSd = SelectSdvec<>;
using SelMcl = SelectMcl<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;


template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class SlpT = SelfShapedSlp<var_t, DagcSd, DagcSd, SelSd>
          >
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:

    // std::vector<size_t> thresholds;
    SlpT slp;

    // std::vector<ulint> samples_start;
    int_vector<> samples_start;
    // int_vector<> samples_end;
    // std::vector<ulint> samples_last;

    // static const uchar TERMINATOR = 1;
    // bool sais = true;
    // /*
    //  * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
    //  */
    // //F column of the BWT (vector of 256 elements)
    // std::vector<ulint> F;
    // //L column of the BWT, run-length compressed
    // rle_string_t bwt;
    // ulint terminator_position = 0;
    // ulint r = 0; //number of BWT runs

    typedef size_t size_type;

    ms_pointers()
        : ri::r_index<sparse_bv_type, rle_string_t>()
        {}

    void build(const std::string& filename) 
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");

        if (true)
        {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        }
        else
        {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }



        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        // int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = " , this->r);
        verbose("Rate n/r = " , double(this->bwt.size()) / this->r);
        verbose("log2(r) = " , log2(double(this->r)));
        verbose("log2(n/r) = " , log2(double(this->bwt.size()) / this->r));



        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);


        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose(3);

//load_grammar(filename);

        verbose("text length: ", slp.getLen());
        verbose("bwt length: ", this->bwt.size());
    }
  

    void load_grammar(const std::string& filename) {
        {
            verbose("Load Grammar");
            std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

            ifstream fs(filename + ".slp");
            slp.load(fs);

            std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
            verbose("Memory peak: ", malloc_count_peak());
            verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        }
        DCHECK_EQ(slp.getLen()+1, this->bwt.size());
    }



    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples)
    {
        int log_n = bitsize(uint64_t(n));
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2*SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1,fd))
        {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    void write_int(ostream& os, const size_t& i){
      os.write(reinterpret_cast<const char*>(&i), sizeof(size_t));
    }

    // Computes the matching statistics pointers for the given pattern
    //std::pair<std::vector<size_t>, std::vector<size_t>> 
    size_t query(const std::string& patternfile, const std::string& len_filename, const std::string& ref_filename) {

      const char* p;
      size_t m;
      map_file(patternfile.c_str(), p, m);

        auto pattern_at = [&] (const size_t pos) {
          return p[pos];
        };

        // ifstream p(patternfile);
        // p.seekg(0,ios_base::end);
        // const size_t m = p.tellg();
        //
        // auto pattern_at = [&] (const size_t pos) {
        //   p.seekg(pos, ios_base::beg);
        //   return p.get();
        // };


        ofstream len_file(len_filename, std::ios::binary);
        ofstream ref_file(ref_filename, std::ios::binary);

        const size_t n = slp.getLen();
        verbose("pattern length: ", m);
        
        size_t last_len;
        size_t last_ref;

        auto write_len = [&] (const size_t i) { 
          last_len = i;
          write_int(len_file, last_len);
        };
        auto write_ref = [&] (const size_t i) { 
          last_ref = i;
          write_int(ref_file, last_ref);
        };

        //TODO: we could allocate the file here and store the numbers *backwards* !
        ON_DEBUG(std::vector<size_t> ms_references(m,0)); //! stores at (m-i)-th entry the text position of the match with pattern[m-i..] when computing the matching statistics of pattern[m-i-1..]
        ON_DEBUG(std::vector<size_t> ms_lengths(m,1));
        ON_DEBUG(ms_lengths[m-1] = 1);
        write_len(1);


        //!todo: we need here a while look in case that we want to support suffixes of the pattern that are not part of the text
        DCHECK_GT(this->bwt.number_of_letter(pattern_at(m-1)), 0);

        //! Start with the last character
        auto pos = this->bwt.select(1, pattern_at(m-1));
        {
            const ri::ulint run_of_j = this->bwt.run_of_position(pos);
            ON_DEBUG(ms_references[m-1] = samples_start[run_of_j]);
            write_ref(samples_start[run_of_j]);
            DCHECK_EQ(slp.charAt(ms_references[m-1]),  pattern_at(m-1));
        }
        pos = LF(pos, pattern_at(m-1));

		#ifdef MEASURE_TIME
		double time_lce = 0;
		double time_backwardstep = 0;
		#endif

        for (size_t i = 1; i < m; ++i) {
            const auto c = pattern_at(m - i - 1);
            DCHECK_EQ(ms_lengths[m-i], last_len);
            ON_DEBUG(DCHECK_EQ(ms_references[m-i], last_ref));

			const size_t number_of_runs_of_c = this->bwt.number_of_letter(c);
            if(number_of_runs_of_c == 0) {
                ON_DEBUG(ms_lengths[m-i-1] = 0);
                write_len(0);
                write_int(len_file, last_len);
                // std::cout << "2 letter " << c  << " not found for " << (m-i-1) << " : " << ms_lengths[m-i-1] << std::endl;
                ON_DEBUG(ms_references[m - i - 1] = 1);
                write_ref(1);
            } 
            else if (pos < this->bwt.size() && this->bwt[pos] == c) {
                DCHECK_NE(i, 0);
                ON_DEBUG(ms_lengths[m-i-1] = ms_lengths[m-i]+1);
                write_len(last_len+1);
                // std::cout << "0 Len for " << (m-i-1) << " : " << ms_lengths[m-i-1] << std::endl;

                DCHECK_GT(ms_references[m - i], 0);
                ON_DEBUG(ms_references[m - i - 1] = ms_references[m - i]-1);
                DCHECK_GT(last_ref, 0);
                write_ref(last_ref-1);
            }
            else {
                const ri::ulint rank = this->bwt.rank(pos, c);

				size_t sa0, sa1;

				if(rank > 0) {
					sa0 = this->bwt.select(rank-1, c);
					DCHECK_LT(sa0, pos);
				}
				if(rank < number_of_runs_of_c) {
					sa1 = this->bwt.select(rank, c);
					DCHECK_GT(sa1, pos);
				}
				
				struct Triplet {
					size_t sa, ref, len;
				};

				auto compute_succeeding_lce = [&] () -> Triplet {
                // if(rank < this->bwt.number_of_letter(c)) { // LCE with the succeeding position
					DCHECK_LT(rank, number_of_runs_of_c);

					const ri::ulint run1 = this->bwt.run_of_position(sa1);

                    const size_t textposStart = this->samples_start[run1];
					#ifdef MEASURE_TIME
					Stopwatch s;
					#endif
                    const size_t lenStart = textposStart+1 >= n ? 0 : lceToRBounded(slp, textposStart+1, last_ref, last_len);
					#ifdef MEASURE_TIME
					time_lce += s.seconds();
					#endif
                    // ON_DEBUG(
                    //         const size_t textposLast = this->samples_last[run1];
                    //         const size_t lenLast = lceToRBounded(slp, textposLast+1, ms_references[m-i]);
                    //         DCHECK_GT(lenStart, lenLast);
                    // )
                    const size_t ref1 = textposStart;
                    const size_t len1 = lenStart;
					return {sa1, ref1, len1};
                };

				auto compute_preceding_lce = [&] () -> Triplet {
                // if(rank > 0) { // LCE with the preceding position
					DCHECK_GT(rank, 0);

					const ri::ulint run0 = this->bwt.run_of_position(sa0);

                    const size_t textposLast = this->samples_last[run0];
					#ifdef MEASURE_TIME
					Stopwatch s;
					#endif
                    const size_t lenLast = textposLast+1 >= n ? 0 : lceToRBounded(slp, textposLast+1, last_ref, last_len);
					#ifdef MEASURE_TIME
					time_lce += s.seconds();
					#endif
                    // ON_DEBUG( //sanity check
                    //         const size_t textposStart = this->samples_start[run0];
                    //         const size_t lenStart = lceToRBounded(slp, textposStart+1, ms_references[m-i], ms_lengths[m-1]);
                    //         DCHECK_LE(lenStart, lenLast);
                    // )
                    const size_t ref0 = textposLast;
                    const size_t len0 = lenLast;
					return {sa0, ref0, len0};
                };

				const Triplet c = [&] () -> Triplet {
					if(rank == 0) {
						return compute_succeeding_lce();
					}
					else if(rank >= number_of_runs_of_c) {
						return compute_preceding_lce();
					}
#ifdef NAIVE_LCE_SCHEDULE 
					{
						const Triplet a = compute_preceding_lce();
						const Triplet b = compute_succeeding_lce();
						if(a.len < b.len) { return b; }
						return a;
					}
#else //NAIVE_LCE_SCHEDULE
#ifdef SORT_BY_DISTANCE_HEURISTIC
					if(pos - sa0 > sa1 - pos) {
#else
					if(true) {
#endif//SORT_BY_DISTANCE_HEURISTIC
						auto eval_first = &compute_preceding_lce;
						auto eval_second = &compute_succeeding_lce;
						const Triplet a = (*eval_first)();
						if(last_len <= a.len) {
							return a;
						} 
						const Triplet b = (*eval_second)();
						if(b.len > a.len) { return b; }
						return a;
					} 
					auto eval_first = &compute_succeeding_lce;
					auto eval_second = &compute_preceding_lce;
					const Triplet a = (*eval_first)();
					if(last_len <= a.len) {
						return a;
					} 
					const Triplet b = (*eval_second)();
					if(b.len > a.len) { return b; }
					return a;
#endif //NAIVE_LCE_SCHEDULE
				}();

                // if(len0 < len1) {
                //     len0 = len1;
                //     ref0 = ref1;
                //     sa0 = sa1;
                // }
                    
				DCHECK_GT(c.ref, 0);
                ON_DEBUG(ms_lengths[m-i-1] = 1 + std::min(ms_lengths[m-i], c.len));
                write_len(1 + std::min(last_len, c.len));
                ON_DEBUG(ms_references[m-i-1] = c.ref);
                write_ref(c.ref);
                pos = c.sa;

                DCHECK_GT(ms_lengths[m-i-1], 0);
                // std::cout << "1 Len for " << (m-i-1) << " : " << ms_lengths[m-i-1] << " with LCE " << std::max(len1,len0) << std::endl;
            }
			#ifdef MEASURE_TIME
			Stopwatch s;
			#endif
            pos = LF(pos, c); //! Perform one backward step
			#ifdef MEASURE_TIME
			time_backwardstep += s.seconds();
			#endif
        }
		#ifdef MEASURE_TIME
		cout << "Time backwardsearch: " << time_backwardstep << std::endl;
		cout << "Time lce: " << time_lce << std::endl;
		#endif
        return m;
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += samples_start.serialize(out, child, "samples_start");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;

    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, const std::string& filename)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->samples_last.load(in);
        this->samples_start.load(in);
        
        load_grammar(filename);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

    protected :

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        {
          ulint i = 0;
          while ((c = heads.get()) != EOF)
          {
            size_t length;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
              this->F[c] += length;
            else
            {
              this->F[TERMINATOR] += length;
              this->terminator_position = i;
            }
            i++;
          }
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }
    };

#endif /* end of include guard: _MS_POINTERS_HH */
