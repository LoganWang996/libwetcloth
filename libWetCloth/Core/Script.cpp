//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Script.h"

#include "DistanceFields.h"
#include "MathUtilities.h"
#include "TwoDScene.h"

using namespace mathutils;

scalar Script::getNextVelocity(const scalar& dt, const scalar& current_time) {
  if (func == Script::COSINE) {
    return cosine_ease_function(current_time + dt, start, end,
                                start + ease_start, end - ease_end, amplitude,
                                frequency);
  } else if (func == Script::CUBIC) {
    return cubic_ease_function(current_time + dt, start, end,
                               start + ease_start, end - ease_end, v(3));
  } else if (func == Script::WENO) {
    const scalar vel = weno_ease_function(current_time + dt, dt, start, end,
                                          base_dt, base_pos, base_vertices);
    base_pos += vel * dt;
    return vel;
  }

  return 0.0;
}

void Script::stepScript(const scalar& dt, const scalar& current_time) {
  if (!m_scene) return;

  if (current_time >= start && current_time + dt <= end) {
    switch (type) {
      case Script::TWIST: {
        m_scene->setStrandForceScript(strand, script_node, v(3));
        break;
      }
      default:
        std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << type << "]!" << std::endl;
        break;
    }
  }
}
