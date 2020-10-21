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

    std::vector<size_t> thresholds;
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

    ms_pointers(std::string filename) : 
        ri::r_index<sparse_bv_type, rle_string_t>()
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");
        std::ifstream ifs(bwt_fname);
        this->bwt = rle_string_t(ifs); 
        // std::string istring;
        // read_file(bwt_fname.c_str(), istring);
        // for(size_t i = 0; i < istring.size(); ++i)
        //     if(istring[i]==0)
        //         istring[i] = TERMINATOR;
        // this->bwt = rle_string_t(istring);

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = " , this->r);
        verbose("Rate n/r = " , double(this->bwt.size()) / this->r);
        verbose("log2(r) = " , log2(double(this->r)));
        verbose("log2(n/r) = " , log2(double(this->bwt.size()) / this->r));

        ifs.seekg(0);
        this->build_F(ifs);
        // this->build_F(istring);
        // istring.clear();
        // istring.shrink_to_fit();


        read_samples(filename + ".ssa", this->r, log_n, samples_start);
        read_samples(filename + ".esa", this->r, log_n, this->samples_last);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


        {
            verbose("Load Grammar");
            std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
            ifstream fs(filename + ".slp");
            slp.load(fs);
            std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
            verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        }


        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();

        std::string tmp_filename = filename + std::string(".thr_pos");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invalid file " + tmp_filename);

        size_t length = filestat.st_size / THRBYTES;
        thresholds.resize(length);

        for(size_t i = 0; i < length; ++i )
            if ((fread(&thresholds[i], THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

        fclose(fd);

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    }

    void read_samples(std::string filename, ulint r, int log_n, int_vector<> &samples)
    {

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
            ulint val = (right ? right - 1 : r - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    // Computes the matching statistics pointers for the given pattern
    std::vector<size_t> query(const std::vector<uint8_t>& pattern)
    {
        size_t m = pattern.size();

        std::vector<size_t> ms_pointers(m);
        std::vector<size_t> my_refs(m,0);
        std::vector<size_t> ms_lengths(m,1);
        std::vector<size_t> debug_lengths(m,1);

        // Start with the empty string
        auto pos = this->bwt.select(1, pattern[m-1]);
        {
            const ri::ulint run_of_j = this->bwt.run_of_position(pos);
            ms_pointers[m-1] = samples_start[run_of_j];
            my_refs[m-1] = samples_start[run_of_j];
            assert(slp.charAt(ms_pointers[m-1]) == pattern[m-1]);
            assert(slp.charAt(my_refs[m-1]) == pattern[m-1]);
            ms_lengths[m-1] = 1;
            debug_lengths[m-1] = 1;
        }
        auto sample = ms_pointers[m-1]; 
        pos = LF(pos, pattern[m-1]);

        for (size_t i = 1; i < pattern.size(); ++i)
        {
            auto c = pattern[m - i - 1];

            if (this->bwt.number_of_letter(c) == 0)
            {
                sample = 0;
                ms_lengths[m-i-1] = 0;
                debug_lengths[m-i-1] = 0;
                std::cout << "2 letter " << c  << " not found for " << (m-i-1) << " : " << ms_lengths[m-i-1] << std::endl;
                my_refs[m - i - 1] = 0;
            } 
            else if (pos < this->bwt.size() && this->bwt[pos] == c)
            {
                sample--;
                assert(i != 0);
                ms_lengths[m-i-1] = ms_lengths[m-i]+1;
                debug_lengths[m-i-1] = debug_lengths[m-i]+1;
                std::cout << "0 Len for " << (m-i-1) << " : " << ms_lengths[m-i-1] << std::endl;

                my_refs[m - i - 1] = my_refs[m - i]-1;
                // {//assert
                //     const ri::ulint run = this->bwt.run_of_position(pos);
                //     const size_t textposStart = this->samples_start[run];
                //     const size_t textposLast = this->samples_last[run];
                //     const size_t lenLast = lceToR(slp, textposLast+1, my_refs[m-i]);
                //     const size_t lenStart = lceToR(slp, textposStart+1, my_refs[m-i]);
                //     const size_t lenMax = std::max(lenStart, lenLast);
                //     assert(lenMax >= ms_lengths[m-i-1]);
                // }
            }
            else
            {
                {
                    const ri::ulint rank = this->bwt.rank(pos, c);
                    size_t len0 = 0;
                    size_t len1 = 0;
                    size_t ref0 = 0;
                    size_t ref1 = 0;

                    //TODO: should check that rank > 0, since select(0) should throw an exception!
                    if(rank < this->bwt.number_of_letter(c)) {
                        const ri::ulint sa1 = this->bwt.select(rank, c);
                        const ri::ulint run1 = this->bwt.run_of_position(sa1);
                        const size_t textposStart = this->samples_start[run1];
                        const size_t textposLast = this->samples_last[run1];
                        const size_t lenStart = lceToR(slp, textposStart+1, my_refs[m-i]);
                        const size_t lenLast = lceToR(slp, textposLast+1, my_refs[m-i]);
                        assert(lenStart >= lenLast);
                        ref1 = textposStart;
                        len1 = lenStart;
                        // len1 = lceToR(slp, textpos1+1, my_refs[m-i]);
                        // len1 = std::max(len1, lceToR(slp, textpos2+1, my_refs[m-i]));
                        // ref = lceToR(slp, textpos1+1, my_refs[m-i]) < 
                    }
                    if(rank > 0) {
                        const ri::ulint sa0 = this->bwt.select(rank-1, c);
                        const ri::ulint run0 = this->bwt.run_of_position(sa0);

                        const size_t textposStart = this->samples_start[run0];
                        const size_t textposLast = this->samples_last[run0];
                        const size_t lenStart = lceToR(slp, textposStart+1, my_refs[m-i]);
                        const size_t lenLast = lceToR(slp, textposLast+1, my_refs[m-i]);
                        assert(lenStart <= lenLast);
                        ref0 = textposLast;
                        len0 = lenLast;
                    }
                    ms_lengths[m-i-1] = 1 + std::min(ms_lengths[m-i],std::max(len1,len0));
                    my_refs[m-i-1] = (len0 > len1) ? ref0 : ref1;
                    // assert(ms_lengths[m-i-1] > 0);
                    std::cout << "1 Len for " << (m-i-1) << " : " << ms_lengths[m-i-1] << " with LCE " << std::max(len1,len0) << std::endl;
                }



                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;

                ulint next_pos = pos;

                // if (rnk < (this->F[c] - this->F[c-1]) // I can use F to compute it
                if (rnk < this->bwt.number_of_letter(c))
                {
                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0

                    // Here we should use Phi_inv that is not implemented yet
                    // sample = this->Phi(this->samples_last[run_of_j - 1]) - 1;
                    sample = samples_start[run_of_j];

                    next_pos = j;
                }

                if (pos < thr)
                {

                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);
                    sample = this->samples_last[run_of_j];

                    next_pos = j;
                }
                pos = next_pos;
            }

            {
                debug_lengths[m-i-1] = 1+std::min(lceToR(slp, sample+1, ms_pointers[m-i]), debug_lengths[m-i]);
                    std::cout << "4 Len for " << (m-i-1) << " : " << debug_lengths[m-i-1] << std::endl;
            }
            ms_pointers[m-i-1] = sample;

            // Perform one backward step
            pos = LF(pos, c);
        }

        return ms_pointers;
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

        written_bytes += my_serialize(thresholds, out, child, "thresholds");
        // written_bytes += my_serialize(samples_start, out, child, "samples_start");
        written_bytes += samples_start.serialize(out, child, "samples_start");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {

        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->pred.load(in);
        this->samples_last.load(in);
        this->pred_to_run.load(in);

        my_load(thresholds,in);
        samples_start.load(in);
        // my_load(samples_start,in);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

    protected :

        // // From r-index
        // vector<ulint> build_F(std::ifstream &ifs)
        // {
        //     ifs.clear();
        //     ifs.seekg(0);
        //     F = vector<ulint>(256, 0);
        //     uchar c;
        //     ulint i = 0;
        //     while (ifs >> c)
        //     {
        //         if (c > TERMINATOR)
        //             F[c]++;
        //         else
        //         {
        //             F[TERMINATOR]++;
        //             terminator_position = i;
        //         }
        //         i++;
        //     }
        //     for (ulint i = 255; i > 0; --i)
        //         F[i] = F[i - 1];
        //     F[0] = 0;
        //     for (ulint i = 1; i < 256; ++i)
        //         F[i] += F[i - 1];
        //     return F;
        // }

        // // From r-index
        // vector<pair<ulint, ulint>> &read_run_starts(std::string fname, ulint n, vector<pair<ulint, ulint>> &ssa)
        // {
        //     ssa.clear();
        //     std::ifstream ifs(fname);
        //     uint64_t x = 0;
        //     uint64_t y = 0;
        //     uint64_t i = 0;
        //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
        //     {
        //         ssa.push_back(pair<ulint, ulint>(y ? y - 1 : n - 1, i));
        //         i++;
        //     }
        //     return ssa;
        // }

        // // From r-index
        // vector<ulint> &read_run_ends(std::string fname, ulint n, vector<ulint> &esa)
        // {
        //     esa.clear();
        //     std::ifstream ifs(fname);
        //     uint64_t x = 0;
        //     uint64_t y = 0;
        //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
        //     {
        //         esa.push_back(y ? y - 1 : n - 1);
        //     }
        //     return esa;
        // }
    };

#endif /* end of include guard: _MS_POINTERS_HH */
