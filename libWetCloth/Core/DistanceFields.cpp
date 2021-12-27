//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifdef WIN32
#define NOMINMAX
#endif

#include "DistanceFields.h"

#include <numeric>

#include "MathUtilities.h"
#include "TwoDScene.h"
#include "Array3.h"
#include "Array3Utils.h"
#include "MakeLevelSet3.h"
#include "Sorter.h"

using namespace mathutils;
DistanceField::DistanceField(DISTANCE_FIELD_TYPE type_,
                             DISTANCE_FIELD_USAGE usage_, int group_,
                             int params_index_, bool sampled_)
    : type(type_),
      usage(usage_),
      parent(NULL),
      group(group_),
      params_index(params_index_),
      sampled(sampled_) {}

void DistanceField::center(Vector3s& cent) const {
  Vector3s low, high;
  local_bounding_box(low, high);
  cent = (low + high) / 2.0;
}

void DistanceField::resample_internal(const std::shared_ptr<TwoDScene>& parent,
                                      const scalar& dx, const VectorXs& exist,
                                      VectorXs& additional) {
  Vector3s bbx_low = Vector3s::Zero();
  Vector3s bbx_high = Vector3s::Zero();

  local_bounding_box(bbx_low, bbx_high);

  bbx_low -= Vector3s(dx, dx, dx);
  bbx_high += Vector3s(dx, dx, dx);

  Vector3s dxyz = (bbx_high - bbx_low) / dx;
  if (dxyz(0) <= 0.0 || dxyz(1) <= 0.0 || dxyz(2) <= 0.0) return;

  Vector3i nxyz =
      Vector3i((int)ceil(dxyz(0)), (int)ceil(dxyz(1)), (int)ceil(dxyz(2)));

  const int num_init = nxyz(0) * nxyz(1) * nxyz(2);
  const int num_exist = exist.size() / 4;
  const int num_total = num_init + num_exist;

  std::vector<int> available(num_total, 0);
  VectorXs pos_init(num_total * 3);

  const int num_nodes = (nxyz(0) + 1) * (nxyz(1) + 1) * (nxyz(2) + 1);

  // init phi
  Array3s phi_init(nxyz(0) + 1, nxyz(1) + 1, nxyz(2) + 1);

  threadutils::for_each(0, num_nodes, [&](int pidx) {
    int iz = pidx / ((nxyz(0) + 1) * (nxyz(1) + 1));
    int iy = (pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1)) / (nxyz(0) + 1);
    int ix = pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1) - iy * (nxyz(0) + 1);
    Vector3s cp = bbx_low + Vector3s(ix, iy, iz) * dx;

    Vector3s vel;

    phi_init(ix, iy, iz) = compute_phi_vel(cp, vel);
  });

  // first pass of selection for new + existing particles
  threadutils::for_each(0, num_total, [&](int pidx) {
    Vector3s cp;

    if (pidx < num_init) {
      int iz = pidx / (nxyz(0) * nxyz(1));
      int iy = (pidx - iz * nxyz(0) * nxyz(1)) / nxyz(0);
      int ix = pidx - iz * nxyz(0) * nxyz(1) - iy * nxyz(0);

      cp = bbx_low + (Vector3s(ix, iy, iz) + Vector3s::Constant(0.5)) * dx +
           Vector3s(mathutils::scalarRand(-dx, dx),
                    mathutils::scalarRand(-dx, dx),
                    mathutils::scalarRand(-dx, dx));
    } else {
      int qidx = pidx - num_init;

      cp = exist.segment<3>(qidx * 4);
    }

    Vector3s vel;

    scalar dist = compute_phi_vel(cp, vel);

    if (dist > dx * sqrt(3.0)) return;

    if (parent->computePhiVel(
            cp, vel, [](const std::shared_ptr<DistanceField>& dfptr) -> bool {
              return dfptr->usage == DFU_SOLID;
            }) < 0.0)
      return;

    Vector3s grad;

    if (pidx < num_init && dist > 0.0) {
      for (int i = 0; i < 5; ++i) {
        dist = compute_phi_vel(cp, vel);

        Vector3s icp = (cp - bbx_low) / dx;

        interpolate_gradient(grad, icp, phi_init);

        if (grad.norm() > 1e-20) grad.normalize();

        cp -= grad * dist;
      }
    }

    pos_init.segment<3>(pidx * 3) = cp;
    available[pidx] = 1;
  });

  // select all particles used to compare
  std::vector<int> mapping(num_total);

  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_parts = mapping[num_total - 1];     // all used particles
  int num_parts_new = mapping[num_init - 1];  // all newly-used particles

  VectorXs pos_selected(num_parts * 3);

  // put all used particles into buffer
  threadutils::for_each(0, num_total, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;
    pos_selected.segment<3>(idx * 3) = pos_init.segment<3>(i * 3);
  });

  available.resize(num_parts);
  for (int i = 0; i < num_parts; ++i) available[i] = 1;

  // sort all used particles in buffer
  Sorter sorter(nxyz(0), nxyz(1), nxyz(2));

  sorter.sort(num_parts, [&](int pidx, int& i, int& j, int& k) {
    i = (int)floor((pos_selected(pidx * 3 + 0) - bbx_low(0)) / dx);
    j = (int)floor((pos_selected(pidx * 3 + 1) - bbx_low(1)) / dx);
    k = (int)floor((pos_selected(pidx * 3 + 2) - bbx_low(2)) / dx);
  });

  // for each new particles, do second pass of selection by dart-picking
  sorter.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        // we ignore existing particles
        if (pidx >= num_parts_new || !available[pidx]) return;

        const Vector3s& pos = pos_selected.segment<3>(pidx * 3);
        sorter.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int np_bucket_idx) -> bool {
              if (pidx == npidx || !available[npidx]) return false;

              const Vector3s& npos = pos_selected.segment<3>(npidx * 3);
              const scalar dist = (npos - pos).norm();

              if (dist < 0.5773502692 * dx) {
                available[pidx] = 0;
                return true;
              }

              return false;
            });
      },
      3);

  // only map new particles
  mapping.resize(num_parts_new);
  std::partial_sum(available.begin(), available.begin() + num_parts_new,
                   mapping.begin());

  int num_results = mapping[num_parts_new - 1];

  additional.resize(num_results * 3);

  // put passed new particles into output
  threadutils::for_each(0, num_parts_new, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;

    additional.segment<3>(idx * 3) = pos_selected.segment<3>(i * 3);
  });
}

void DistanceField::sample(const scalar& dx, VectorXs& result,
                           VectorXs& normals) {
  Vector3s bbx_low = Vector3s::Zero();
  Vector3s bbx_high = Vector3s::Zero();

  local_bounding_box(bbx_low, bbx_high);

  bbx_low -= Vector3s(dx, dx, dx);
  bbx_high += Vector3s(dx, dx, dx);

  Vector3s dxyz = (bbx_high - bbx_low) / dx;
  if (dxyz(0) <= 0.0 || dxyz(1) <= 0.0 || dxyz(2) <= 0.0) return;

  Vector3i nxyz =
      Vector3i((int)ceil(dxyz(0)), (int)ceil(dxyz(1)), (int)ceil(dxyz(2)));

  const int num_init = nxyz(0) * nxyz(1) * nxyz(2);
  std::vector<int> available(num_init, 0);
  VectorXs pos_init(num_init * 3);

  const int num_nodes = (nxyz(0) + 1) * (nxyz(1) + 1) * (nxyz(2) + 1);

  Array3s phi_init(nxyz(0) + 1, nxyz(1) + 1, nxyz(2) + 1);

  threadutils::for_each(0, num_nodes, [&](int pidx) {
    int iz = pidx / ((nxyz(0) + 1) * (nxyz(1) + 1));
    int iy = (pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1)) / (nxyz(0) + 1);
    int ix = pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1) - iy * (nxyz(0) + 1);
    Vector3s cp = bbx_low + Vector3s(ix, iy, iz) * dx;

    Vector3s vel;

    phi_init(ix, iy, iz) = compute_phi_vel(cp, vel);
  });

  threadutils::for_each(0, num_init, [&](int pidx) {
    int iz = pidx / (nxyz(0) * nxyz(1));
    int iy = (pidx - iz * nxyz(0) * nxyz(1)) / nxyz(0);
    int ix = pidx - iz * nxyz(0) * nxyz(1) - iy * nxyz(0);

    Vector3s cp =
        bbx_low + Vector3s(ix, iy, iz) * dx +
        Vector3s(mathutils::scalarRand(0.0, dx), mathutils::scalarRand(0.0, dx),
                 mathutils::scalarRand(0.0, dx));

    Vector3s vel;

    scalar dist = compute_phi_vel(cp, vel);

    if (fabs(dist) > dx * sqrt(3.0)) return;

    Vector3s grad;

    for (int i = 0; i < 5; ++i) {
      dist = compute_phi_vel(cp, vel);

      Vector3s icp = (cp - bbx_low) / dx;

      interpolate_gradient(grad, icp, phi_init);

      if (grad.norm() > 1e-20) grad.normalize();

      cp -= grad * dist;
    }

    pos_init.segment<3>(pidx * 3) = cp;
    available[pidx] = 1;
  });

  std::vector<int> mapping(num_init);

  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_parts = mapping[mapping.size() - 1];

  VectorXs pos_selected(num_parts * 3);

  threadutils::for_each(0, num_init, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;
    pos_selected.segment<3>(idx * 3) = pos_init.segment<3>(i * 3);
  });

  available.resize(num_parts);
  for (int i = 0; i < num_parts; ++i) available[i] = 1;

  Sorter sorter(nxyz(0), nxyz(1), nxyz(2));

  sorter.sort(num_parts, [&](int pidx, int& i, int& j, int& k) {
    i = (int)floor((pos_selected(pidx * 3 + 0) - bbx_low(0)) / dx);
    j = (int)floor((pos_selected(pidx * 3 + 1) - bbx_low(1)) / dx);
    k = (int)floor((pos_selected(pidx * 3 + 2) - bbx_low(2)) / dx);
  });

  sorter.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        if (!available[pidx]) return;

        const Vector3s& pos = pos_selected.segment<3>(pidx * 3);
        sorter.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int np_bucket_idx) -> bool {
              if (pidx == npidx || !available[npidx]) return false;

              const Vector3s& npos = pos_selected.segment<3>(npidx * 3);
              const scalar dist = (npos - pos).norm();

              if (dist < 2.5980762116 * dx) {
                available[pidx] = 0;
                return true;
              }

              return false;
            });
      },
      3);

  mapping.resize(num_parts);
  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_results = mapping[mapping.size() - 1];

  result.resize(num_results * 3);

  normals.resize(num_results * 3);

  threadutils::for_each(0, num_parts, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;

    result.segment<3>(idx * 3) = pos_selected.segment<3>(i * 3);
    Vector3s icp = (result.segment<3>(idx * 3) - bbx_low) / dx;

    Vector3s grad;

    interpolate_gradient(grad, icp, phi_init);

    if (grad.norm() > 1e-20) grad.normalize();

    normals.segment<3>(idx * 3) = grad;
  });
}

bool DistanceFieldObject::check_durations(const scalar& cur_time,
                                          const scalar& cur_vol,
                                          Vector3s& shooting_vel) {
  bool ret = false;
  for (auto& dur : durations) {
    ret =
        (cur_time >= dur.start && cur_time <= dur.end) && cur_vol < dur.maxvol;
    if (ret) {
      shooting_vel = dur.vel;
      break;
    }
  }

  return ret;
}

DistanceFieldObject::DistanceFieldObject(
    const Vector3s& center_, const VectorXs& parameter_,
    DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_, bool inside,
    const Vector3s& raxis, const scalar& rangle, int group_, int params_index_,
    bool sampled_, const std::vector<DF_SOURCE_DURATION>& durations_,
    const std::string& szfn, const std::string& szfn_cache)
    : DistanceField(type_, usage_, group_, params_index_, sampled_),
      center(center_),
      parameter(parameter_),
      sign(inside ? -1.0 : 1.0),
      rot(Eigen::AngleAxis<scalar>(rangle, raxis)),
      durations(durations_) {
  V.setZero();
  omega.setZero();

  switch (type_) {
    default:
      mesh = nullptr;
  }

  future_center = center;
  future_rot = rot;
}

void DistanceFieldObject::resample_mesh(const scalar& dx, VectorXs& result,
                                        VectorXs& normals) {
  const int num_tris = (int)mesh->getIndices().size();
  std::vector<scalar> PDF(num_tris);

  auto& indices = mesh->getIndices();
  auto& verts = mesh->getVertices();

  scalar tot_area = 0.0;

  for (int i = 0; i < num_tris; ++i) {
    const Vector3s& v0 = verts[indices[i](0)];
    const Vector3s& v1 = verts[indices[i](1)];
    const Vector3s& v2 = verts[indices[i](2)];

    const scalar area = (v2 - v0).cross(v1 - v0).norm();

    PDF[i] = area;
    tot_area += area;
  }

  if (tot_area == 0.0) return;

  std::partial_sum(PDF.begin(), PDF.end(), PDF.begin());
  for (int i = 0; i < num_tris; ++i) {
    PDF[i] /= tot_area;
  }

  tot_area *= 0.5;

  const scalar coverage_radius =
      dx * mathutils::defaultRadiusMultiplier() * 0.5;
  const scalar coverage_area = M_PI * coverage_radius * coverage_radius;

  const int num_samples = (int)ceil(tot_area / coverage_area);

  result.resize(num_samples * 3);
  normals.resize(num_samples * 3);

  Eigen::AngleAxis<scalar> rotaa(rot);

  for (int i = 0; i < num_samples; ++i) {
    const scalar seed0 = mathutils::scalarRand(0.0, 1.0);
    std::vector<scalar>::iterator low =
        std::lower_bound(PDF.begin(), PDF.end(), seed0);
    const int tri_idx = (int)(low - PDF.begin());

    if (tri_idx < 0 || tri_idx >= num_tris) continue;

    const scalar s1 = mathutils::scalarRand(0.0, 1.0);
    const scalar s2 = mathutils::scalarRand(0.0, 1.0);

    const Vector3s& x0 = verts[indices[tri_idx](0)];
    const Vector3s& x1 = verts[indices[tri_idx](1)];
    const Vector3s& x2 = verts[indices[tri_idx](2)];

    const scalar ss1 = sqrt(s1);

    const Vector3s p =
        (1.0 - ss1) * x0 + (ss1 * (1.0 - s2)) * x1 + (s2 * ss1) * x2;

    const Vector3s n = (x2 - x0).cross(x1 - x0).normalized();

    result.segment<3>(i * 3) = rotaa * p + center;
    normals.segment<3>(i * 3) = rotaa * n * sign;
  }
}

void DistanceFieldObject::apply_global_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  future_rot = rot * future_rot;
}

void DistanceFieldObject::apply_local_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  future_rot = future_rot * rot;
}

void DistanceFieldObject::apply_translation(const Vector3s& t) {
  future_center += t;
}

scalar DistanceFieldObject::compute_phi(const Vector3s& pos) const {
  scalar phi = 0.0;

  Vector3s dx = pos - center;

  switch (type) {
    case DFT_FILE: {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp_coord =
          ((irot * p0 * irot.inverse()).vec() - volume_origin) / parameter(1);
      phi = sign * interpolate_value(rotp_coord, volume);
      break;
    }
    default:
      break;
  }

  return phi;
}

scalar DistanceFieldObject::compute_phi_vel(const Vector3s& pos,
                                            Vector3s& vel) const {
  scalar phi = compute_phi(pos);

  Vector3s dx = pos - center;

  vel = V + omega.cross(dx);

  return phi;
}

void DistanceFieldObject::advance(const scalar& dt) {
  V = (future_center - center) / dt;
  Eigen::Quaternion<scalar> q = future_rot * rot.conjugate();
  scalar len = q.vec().norm();
  if (len > 0.0) {
    scalar angle = 2.0 * atan2(len, q.w());
    omega = q.vec() / len * angle / dt;
  } else {
    omega.setZero();
  }

  center = future_center;
  rot = future_rot;
}

bool DistanceFieldObject::local_bounding_box(Vector3s& bbx_low,
                                             Vector3s& bbx_high) const {
  switch (type) {
    case DFT_FILE:
      mesh->boundingBox(bbx_low, bbx_high, rot);
      bbx_low += center;
      bbx_high += center;
      break;
    default:
      break;
  }

  return sign < 0.0;
}

void DistanceFieldObject::render(
    const std::function<void(const std::vector<Vector3s>&,
                             const std::vector<Vector3i>&,
                             const Eigen::Quaternion<scalar>&, const Vector3s&,
                             const scalar&)>& func) const {
  if (!mesh) return;

  func(mesh->getVertices(), mesh->getIndices(), rot, center, sign);
}

DistanceFieldOperator::DistanceFieldOperator(DISTANCE_FIELD_TYPE type_,
                                             DISTANCE_FIELD_USAGE usage_,
                                             int group_, int params_index_,
                                             bool sampled_)
    : DistanceField(type_, usage_, group_, params_index_, sampled_) {}

bool DistanceFieldOperator::check_durations(const scalar& cur_time,
                                            const scalar& cur_vol,
                                            Vector3s& shooting_vel) {
  int nb = children.size();
  for (int i = 0; i < nb; ++i) {
    if (children[i]->check_durations(cur_time, cur_vol, shooting_vel))
      return true;
  };

  return false;
}

void DistanceFieldOperator::apply_global_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  int nb = children.size();
  threadutils::for_each(
      0, nb, [&](int i) { children[i]->apply_global_rotation(rot); });
}

void DistanceFieldOperator::resample_mesh(const scalar& dx, VectorXs& result,
                                          VectorXs& normals) {
  int nb = children.size();
  threadutils::for_each(0, nb, [&](int i) {
    VectorXs tmp_result;
    VectorXs tmp_norm;
    children[i]->resample_mesh(dx, tmp_result, tmp_norm);
    const int old_size = result.size();
    result.conservativeResize(result.size() + tmp_result.size());
    normals.conservativeResize(normals.size() + tmp_norm.size());

    result.segment(old_size, tmp_result.size()) = tmp_result;
    normals.segment(old_size, tmp_norm.size()) = tmp_norm;
  });
}

void DistanceFieldOperator::apply_local_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  int nb = children.size();
  threadutils::for_each(0, nb,
                        [&](int i) { children[i]->apply_local_rotation(rot); });
}

void DistanceFieldOperator::apply_translation(const Vector3s& t) {
  int nb = children.size();
  threadutils::for_each(0, nb,
                        [&](int i) { children[i]->apply_translation(t); });
}

void DistanceFieldOperator::advance(const scalar& dt) {
  int nb = children.size();
  threadutils::for_each(0, nb, [&](int i) { children[i]->advance(dt); });
}

void DistanceFieldOperator::render(
    const std::function<void(const std::vector<Vector3s>&,
                             const std::vector<Vector3i>&,
                             const Eigen::Quaternion<scalar>&, const Vector3s&,
                             const scalar&)>& func) const {
  int nb = children.size();
  for (int i = 0; i < nb; ++i) {
    children[i]->render(func);
  }
}

scalar DistanceFieldOperator::compute_phi(const Vector3s& pos) const {
  switch (DistanceField::type) {
    case DFT_UNION: {
      scalar min_phi = 1e+20;

      for (auto& child : children) {
        scalar phi = child->compute_phi(pos);
        if (phi < min_phi) {
          min_phi = phi;
        }
      }

      return min_phi;
    }
    case DFT_INTERSECT: {
      scalar max_phi = -1e+20;

      for (auto& child : children) {
        scalar phi = child->compute_phi(pos);
        if (phi > max_phi) {
          max_phi = phi;
        }
      }

      return max_phi;
    }
    default:
      return 1e+20;
  }
}

scalar DistanceFieldOperator::compute_phi_vel(const Vector3s& pos,
                                              Vector3s& vel) const {
  switch (DistanceField::type) {
    case DFT_UNION: {
      scalar min_phi = 1e+20;

      for (auto& child : children) {
        Vector3s sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if (phi < min_phi) {
          min_phi = phi;
          vel = sub_vel;
        }
      }

      return min_phi;
    }
    case DFT_INTERSECT: {
      scalar max_phi = -1e+20;

      for (auto& child : children) {
        Vector3s sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if (phi > max_phi) {
          max_phi = phi;
          vel = sub_vel;
        }
      }

      return max_phi;
    }
    default:
      vel = Vector3s::Zero();
      return 1e+20;
  }
}

int DistanceFieldOperator::vote_param_indices() {
  int m = params_index;
  int i = 0;
  for (auto& child : children) {
    int cpi = child->vote_param_indices();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  params_index = m;

  return m;
}

DISTANCE_FIELD_USAGE DistanceFieldOperator::vote_usage() {
  DISTANCE_FIELD_USAGE m = usage;
  int i = 0;
  for (auto& child : children) {
    DISTANCE_FIELD_USAGE cpi = child->vote_usage();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  usage = m;

  return m;
}

bool DistanceFieldOperator::vote_sampled() {
  bool m = sampled;
  int i = 0;
  for (auto& child : children) {
    bool cpi = child->vote_sampled();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  sampled = m;

  return m;
}

bool DistanceFieldOperator::local_bounding_box(Vector3s& bbx_low,
                                               Vector3s& bbx_high) const {
  bool inside = false;

  switch (type) {
    case DFT_UNION:
      bbx_low = Vector3s::Constant(1e+20);
      bbx_high = Vector3s::Constant(-1e+20);

      for (auto& child : children) {
        Vector3s bbx_l, bbx_h;
        inside = inside || child->local_bounding_box(bbx_l, bbx_h);

        for (int r = 0; r < 3; ++r) {
          bbx_low(r) = std::min(bbx_low(r), bbx_l(r));
          bbx_high(r) = std::max(bbx_high(r), bbx_h(r));
        }
      }
      break;
    case DFT_INTERSECT:
      bbx_low = Vector3s::Constant(-1e+20);
      bbx_high = Vector3s::Constant(1e+20);

      for (auto& child : children) {
        Vector3s bbx_l, bbx_h;
        inside = inside || child->local_bounding_box(bbx_l, bbx_h);

        for (int r = 0; r < 3; ++r) {
          bbx_low(r) = std::max(bbx_low(r), bbx_l(r));
          bbx_high(r) = std::min(bbx_high(r), bbx_h(r));
        }
      }
      break;

    default:
      break;
  }

  return inside;
}
