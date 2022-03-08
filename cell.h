//
//  cell.h
//  neutralizer_backbone
//
//  Created by Thijs Janzen on 12/05/2021.
//  Copyright © 2021 Thijs Janzen. All rights reserved.
//

#ifndef cell_h
#define cell_h

#include "rand_t.h"
#include <array>
#include <cassert>

struct species {
  size_t id_;
  int count_;

  species(int count, rnd_t& rndgen) : count_(count) {
    id_ = rndgen.random_number(static_cast<size_t>(1e10));
    for(int i = 0; i < 3; ++i) {
      color_[i] = rndgen.random_number(256); // in range [0, 255]
    }
  }

  species() {
    count_ = 0;
    // fake data, using default rng:
    id_ = std::rand();
    for(int i = 0; i < 3; ++i) {
      color_[i] = std::rand() % 256;
    }
  }

  bool operator==(const species& other) const {
    if (id_ != other.id_) return false;
    if (color_[0] != other.color_[0]) return false;
    if (color_[1] != other.color_[1]) return false;
    if (color_[2] != other.color_[2]) return false;
    return true;
  }

  const std::array<size_t, 3>& get_color() const noexcept{
    return color_;
  }

private:
  std::array<size_t, 3> color_;
};




struct cell {
  species local_species;
  size_t x_;
  size_t y_;

  cell() {
    local_species = species();
  }

  void set_species(const species& target_species) {
    local_species = target_species;
  }

  const species& get_species() const noexcept {
    return local_species;
  }

  int get_species_id() const noexcept {
    return local_species.id_;
  }

  void set_xy(size_t pos, size_t row_size) {
    x_ = pos / row_size;
    y_ = pos % row_size;
  }

  // we assume we have a fixed grid
  cell(cell&&) = delete;
  const cell& operator=(cell&&) = delete;
  cell(const cell&) = delete;
  const cell& operator=(const cell&) = delete;
};


#endif /* cell_h */
