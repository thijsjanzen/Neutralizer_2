//
//  simulation.h
//  neutralizer_backbone
//
//  Created by Thijs Janzen on 12/05/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef simulation_h
#define simulation_h

#include <vector>
#include "cell.h"
#include "rand_t.h"
#include <algorithm>
#include <map>
#include <cmath>

class simulation {
private:


  std::vector< cell > world;
  std::vector< species > meta_community;
  std::vector<int> meta_community_octaves;
  std::map<int, int> histogram_local_comm;
  std::vector<int> local_community_octaves;


  std::vector<double> cdf_;

  rnd_t rndgen_;

  const double prob_same;
  const double rel_prob_spec;

  const double dispersal_range;

public:
  size_t t;
    size_t L;
  std::vector<double> rank_abund_curve;

  simulation(size_t one_side,
             double sp,
             double mgr,
             size_t meta_comm_size,
             double disp_range,
             double theta,
             bool init_mono_dom) :
    L(one_side),
    world(one_side * one_side),
    prob_same(1.0 - sp - mgr),
    rel_prob_spec(sp / (sp + mgr)),
    dispersal_range(disp_range),
    t(0)
  {
    rndgen_ = rnd_t();
    rndgen_.set_world_size(one_side * one_side);
    create_meta_community(meta_comm_size, theta);
    size_t cnt = 0;

    auto mono_dom_spec = get_species_from_meta_community();

    for (auto& i : world) {
        if (init_mono_dom) {
            i.set_species(mono_dom_spec);
          } else {
            i.set_species( get_species_from_meta_community() );
          }

        i.set_xy(cnt, L);
        cnt++;
      }
  }


  std::vector < species > create_meta_community(size_t Jm, double theta) {

    std::vector<int> abund(1,0);
    std::size_t nsp = 1;
    abund.push_back(1);

    for(size_t j = 1; j < Jm; ++j) {
        double x = static_cast<double>(rndgen_.uniform());
        double val = theta / (theta + j - 1);
        if(x < val) {
            nsp++;
            if(nsp > (abund.size()-1)) {
                size_t dif = 1 + nsp - abund.size();
                for(size_t k = 0; k < dif; ++k) {
                    abund.push_back(0);
                  }
              }
            abund[nsp] = 1;
          }
        else {
            int translate_to_abund = static_cast<int>((x * j) - 1);
            //now find corresponding species
            std::size_t index = 0;
            while(index < abund.size()) {
                translate_to_abund -= abund[index];
                if(translate_to_abund <= 0) break;

                index++;
              }

            abund[index] = abund[index] + 1;
          }
      }

    meta_community.clear();

    for(std::size_t i = 0; i < abund.size() ;++i) {
      if (abund[i] > 0) {
          meta_community.push_back(species(abund[i], rndgen_));
        }
    }

    double cumsum = 0.0;
    for (const auto& i : meta_community) {
        cumsum += i.count_;
        cdf_.push_back(cumsum);
      }

    update_octave_meta_comm();

    return meta_community;
  }

  species get_species_from_meta_community() {
    double p =  rndgen_.uniform() * cdf_.back(); //rndgen_.random_number(cdf_.back());
    size_t index = std::distance(cdf_.cbegin(), std::lower_bound(cdf_.cbegin(), cdf_.cend(), p));
    return meta_community[index];
  }

  size_t convert_to_pos(size_t x, size_t y) {
    size_t output =  x * L + y;
    if (output >= world.size()) {
        output = world.size() - 1;
    }
    return output;
  }

  size_t get_coordinate(size_t source_x,
                        size_t source_y) {

    static const float Pi = 3.14159265359f;
    int distance = 1 + static_cast<int>(rndgen_.uniform() * dispersal_range);
    float dir = rndgen_.uniform() * 2 * Pi;
    double pY = static_cast<double>(sinf(dir) * distance);
    double pX = static_cast<double>(cosf(dir) * distance);

    int target_x = static_cast<int>(std::round(source_x + pX)); // round to get values >0.5 to be round up (or < -0.5 round down)
    int target_y = static_cast<int>(std::round(source_y + pY));

    int max_val = static_cast<int>(L);
    if (target_x < 0)  target_x += max_val;
    if (target_x >= max_val) target_x -= max_val;
    if (target_y < 0)  target_y += max_val;
    if (target_y >= max_val) target_y -= max_val;

    if (target_x == source_x && target_y == source_y) {
        return get_coordinate(source_x, source_y);
    }

    return convert_to_pos(target_x, target_y);
  }



  species local_reproduction(size_t source_index) {
    auto pos = get_coordinate(world[source_index].x_,
                              world[source_index].y_);

    return world[pos].get_species();
  }


  void update() {
    size_t pos_to_die = rndgen_.random_pos();

    if (rndgen_.bernouilli(prob_same)) {
        // reproduce locally
        world[pos_to_die].set_species( local_reproduction( pos_to_die ) );
      } else {
        if (rndgen_.bernouilli(rel_prob_spec)) {
            // speciation
            world[pos_to_die].set_species( species(1, rndgen_));
          } else {
            // migration
            world[pos_to_die].set_species( get_species_from_meta_community());
          }
      }
    t++;
  }

  std::array<size_t, 3> get_color(size_t pos) const {
    return world[pos].get_species().get_color();
  }

  size_t update_stats() {
    histogram_local_comm.clear();
    for (const auto& i : world) {
        ++histogram_local_comm[i.get_species_id()];
      }

    local_community_octaves.clear();
    local_community_octaves = std::vector<int>(1 + static_cast<int>(log2(world.size())), 0);
    for(const auto& i : histogram_local_comm) {
        auto oct = octave_sort(i.second);
        local_community_octaves[oct]++;
      }

    update_rank_abund_curve();

    return histogram_local_comm.size();
  }

  void update_rank_abund_curve() {
    rank_abund_curve = std::vector<double>(histogram_local_comm.size());
    int cnt = 0;
    double max = -1;
    for (const auto& i : histogram_local_comm) {
        rank_abund_curve[cnt] = i.second;
        if (i.second > max) max = i.second;
        cnt++;
      }
    std::sort(rank_abund_curve.begin(), rank_abund_curve.end(), std::greater<double>());
    double mult = 100.0 / max;
    for(auto& i : rank_abund_curve) {
        i *= mult;
      }
  }

  size_t octave_sort(long ab_in) { //adapted from James
    size_t result;
    if(ab_in <= 0) {
        result = 0;
      } else {
        long min = 1;
        long max = 2;
        result = 0;
        while(!((ab_in < max)&&(ab_in >= min))) {
            min = min*2;
            max = max*2;
            result ++;
          }
      }
    return result;
  }

  void update_octave_meta_comm() {
    /*std::map<int, int> histogram;
    for (const auto& i : meta_community) {
      ++histogram[i.id_];
    }*/

    // should not go beyond 2^100 normally...
    meta_community_octaves = std::vector<int>(100, 0);
    for(const auto& i : meta_community) {
        auto oct = octave_sort(i.count_);
        meta_community_octaves[oct]++;
      }

    // remove trailing zeros:
    while(meta_community_octaves.back() == 0) {
        meta_community_octaves.pop_back();
      }
  }

  std::vector< int > get_meta_octaves() {
    return meta_community_octaves;
  }

  std::vector< int > get_local_octaves() {
    if (local_community_octaves.empty()) {
        update_stats();
      }
    return local_community_octaves;
  }

  int num_species() {
    return histogram_local_comm.size();
  }

  void update_species_area(std::vector< double >& area,
                           std::vector< double >& num_species) {

    area.clear();
    num_species.clear();
    std::vector< species > found_species;

    for(size_t x = 0; x < L; ++x) {
        for(size_t y = 0; y <= x; ++y) {
            auto pos = convert_to_pos(x, y );
            bool match = false;
            auto local = world[pos].get_species();
            for (const auto& i : found_species) {
              if(i == local) {
                    match = true;
                    break;
              }
            }

            if(!match) {
                found_species.push_back(local);
            }
          }
        if(x > 0) {
            area.push_back( static_cast<double>(x * x) );
            num_species.push_back(static_cast<double>(found_species.size()));
        }
      }
    return;
  }
};




#endif /* simulation_h */
